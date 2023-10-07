# File: QI2report.py
# Author: jiang tao
# Time: 2023/2/8 12:28
# Desc: 中药成分鉴定报告
import copy
import datetime
import functools
import logging
import os
import re
import sys
import tempfile
import uuid

import requests
import yaml
import shutil
import zipfile
from typing import Tuple, Union

from scipy.signal import find_peaks
from sqlalchemy.engine import Engine
from tqdm import tqdm

import numpy as np
import pandas as pd
import subprocess as sp
import matplotlib as mpl
from pathlib import Path
from PIL import Image

from .exception import SaveOnlineRegisterTableError, RawDataNotFoundError, ShellRunTimeError
from .Constants import _engine, BASE_URL, TOKEN, MS_CONVERT

from matplotlib import pyplot as plt, ticker
from matplotlib.ticker import MultipleLocator
from matplotlib import font_manager
from matplotlib.colors import ListedColormap
from pyopenms import MSExperiment, MzMLFile, MSSpectrum, SpectrumAlignmentScore
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from sqlalchemy import create_engine, text
from matplotlib.axes import Axes

from adjustText import adjust_text
from oebio.report import Report as oebioReport

# 5 / 8  248/100
FIG_ASPECT = (16*2.7) / 100
logger = logging.getLogger("TCMreport")


def ax_extra_settings(bottom: int):
    """This decorator is deprecated."""

    def out_wrapper(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            func(*args, **kwargs)
            ax = args[1]
            ax.yaxis.set_major_locator(ticker.AutoLocator())
            ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax.set_axisbelow(True)
            y_ticks = ax.get_yticks()
            ax.set_yticks(y_ticks[(y_ticks <= 100) & (y_ticks >= bottom)])
            ax.yaxis.set_major_formatter(
                ticker.FuncFormatter(lambda x, pos: f"{abs(int(x))}")
            )
            ax.tick_params(axis="both", which="both", labelsize="small")

        return wrapper

    return out_wrapper


class MS2spec:
    def __init__(self, mz: list, intensity: list, metadata: dict):
        self.mz = np.asarray(mz)
        self.intensity = np.asarray(intensity)
        self.metadata = metadata

    def __add__(self, other):
        return [self, other]

    def __len__(self):
        return 1

    def __repr__(self):
        return self.metadata["Comment"] + self.metadata.get("Precursor_type", "None")


class Report:
    """
    从原始数据到oebio_report的完整流程工具类
    Attention: make sure your mzml path do not contains "Chinese char"!
    """
    ENGINE = create_engine(_engine)
    group = ("Control", "Treatment", "TCM")

    # MASS * nM / charge[1] + deltaM = mz
    QI_ADDUCT = {
        "POS": {
            "M+H": 1.007276,
            "M+Na": 22.989218,
            "M+NH4": 18.033823,
            "M+K": 38.963158,
            "M+H-H2O": 17.0027
                        },
        "NEG": {
            "M-H": -1.007276,
            "M+FA-H": 44.998199,
            "M-H2O-H": -19.018391,
            "2M-H": -1.007276
        }
    }
    database_priority = ('TCM', 'ANIMAL', 'HERB', 'PRODUCT')

    sql1 = """SELECT
                    InChIKey,
                    SMILES
                FROM
                    herb_network
                WHERE
                    network_id = "%s"
            """
    sql2 = """SELECT
                    b.inchikey,
                    b.smiles,
                    c.`name`,
                    c.hmdbid,
                    c.metlinid,
                    c.lipidmapsid,
                    c.keggid,
                    c.pubchemid,
                    c.casNumber
                FROM
                    compound_structure AS b
                    INNER JOIN
                    compound_identifier AS c
                    ON
                        b.cid = c.cid
                WHERE
                    b.inchikey = "%s"
            """
    sql3 = """SELECT
                    a.cn_name,
                    a.inchikey,
                    a.cn_class1,
                    a.cn_class2,
                    a.factory,
                    a.item_number,
                    a.purity,
                    a.class1,
                    a.class2
                FROM
                    herb_annotation AS a
                WHERE
                    a.inchikey = "%s"
            """
    sql4 = """SELECT
                    a.inchikey, 
                    a.smiles, 
                    b.hmdbid, 
                    b.metlinid, 
                    b.lipidmapsid, 
                    b.keggid, 
                    d.chebiid,
                    b.pubchemid, 
                    b.casNumber, 
                    c.superclass, 
                    c.class, 
                    c.subclass,
                    a.cid
                FROM
                    compound_structure AS a
                    INNER JOIN
                    compound_identifier AS b
                    ON 
                        a.cid = b.cid
                    LEFT JOIN
                    compound_classification AS c
                    ON 
                        b.cid = c.cid
                    LEFT JOIN
                    extra_hmdb AS d
                    ON 
                        c.cid = d.cid
                WHERE
                    a.inchikey = "%s"
    """

    def __init__(self, data, project: str = "TCM", form: dict = None, local=False):
        self.data = data
        assert project in ("TCM", "Blood"), "Value Error: 项目类型必须为 TCM 或 Blood。"
        self._project = project
        self.form = form
        self.midware = Path(self.data) / "midware"
        self.local = local

    @staticmethod
    def sci_big(string):
        main, sci = string.split("e")
        _ = str(np.ceil(float(main)))
        return float(_ + 'e' + sci)

    @staticmethod
    def get_online_register_form(analysis_id: str) -> dict:
        resp = requests.post(BASE_URL, data=dict(analysis_id=analysis_id, token=TOKEN), timeout=10)
        if resp.status_code != 200:
            raise ValueError(f"不存在的分析编号：{analysis_id}")
        r: dict = resp.json()
        resp.close()
        # todo: parse valid field and replace sample info
        key = r['分析基本信息']['key'].values()
        values = r['分析基本信息']['value'].values()
        parse_data = dict(zip(key, values))
        return parse_data

    @staticmethod
    def estimate_noise(arr: np.ndarray):
        """
         calc noise
        :param arr: peak_heights
        :return: relative noise 1-100
        """
        if arr.size <= 30:
            return 0.1
        counts, bins = np.histogram(arr, bins=12)
        max_count_index = np.argmax(counts)
        densest_range_start = bins[max_count_index]
        densest_range_end = bins[max_count_index + 1]
        return arr[np.where((densest_range_start <= arr) & (arr <= densest_range_end))].mean()

    @staticmethod
    def subset_ms1(rts, its, expect_rt, interval=0.2):
        x = rts[np.where((rts >= expect_rt - interval) & (rts <= expect_rt + interval))]
        y = its[np.where((rts >= expect_rt - interval) & (rts <= expect_rt + interval))]
        return x, y

    @staticmethod
    def is_ms1_complex(its: np.ndarray):
        return len(find_peaks(its, height=5, width=1)[0]) >= 24

    def check_ms1(self, rts: np.ndarray, its: np.ndarray, expect_rt: float) -> bool:
        if self.is_ms1_complex(rts):
            rts, its = self.subset_ms1(rts, its, expect_rt, interval=0.5)
            peak_idx, peak_info = find_peaks(rts, height=0.1, width=(1.4, 15+5))
            peak_noise = peak_info["peak_heights"] / (peak_info["peak_heights"] - peak_info["prominences"])
            r = (its[peak_idx][peak_noise >= 3] - expect_rt) <= 0.2
            return r.any()
        else:
            # ndarray, dict
            peak_idx, peak_info = find_peaks(its, height=0.05)
            noise = self.estimate_noise(peak_info['peak_heights'])
            rts, its = self.subset_ms1(rts, its, expect_rt)
            # 0.05-0.5
            peak_idx, peak_info = find_peaks(its, height=noise, width=(1.4, 15+5))
            # print(rts[peak_idx])
            # print(peak_info)
            # print(noise)
            snr = peak_info['peak_heights'] / noise
            # print(snr)
            return np.any(snr >= 3)

    def rm_specific_cpd(self):
        """deprecated"""
        logger.info("正在从表格中删除指定的代谢物信息...")
        project_name, project_info = self.parse_project()

        data_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "数据矩阵.xlsx")
        del_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "delCpds.xlsx")
        pic_path = str(Path(self.data) / project_name / "数据分析" / "3.质谱图" / "png")
        pdf_path = str(Path(self.data) / project_name / "数据分析" / "3.质谱图" / "pdf")

        delCpds = pd.read_excel(del_path, sheet_name=0)["Metabolites"].tolist()
        data = pd.read_excel(data_path, sheet_name=0).set_index(keys=["Metabolites", ], drop=False)
        fail_list = []
        for cpd in delCpds:
            # 删除图片和 pdf
            filename = self.filename_handler(cpd, escape=True)
            logger.info(f"正在从删除代谢物{cpd}...")
            try:
                os.unlink(pic_path+'/'+filename+'.png')
                os.unlink(pdf_path+'/'+filename+'.pdf')
            except Exception as e:
                print(str(e))
                fail_list.append(filename)
        print(fail_list)
        #data = data.query("Metabolites not in @delCpds").reset_index(drop=True)
        #data.to_excel(data_path, index=False)

    @classmethod
    def find_best_adduct(cls, mz, smi, mode):
        # MASS * nM / charge[1] + deltaM = mz
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        exact_mass = Descriptors.ExactMolWt(mol)
        for k, v in cls.QI_ADDUCT[mode.upper()].items():
            deltaM = (mz - exact_mass * 2) if k == "2M-H" else (mz - exact_mass)
            _ = deltaM - v
            if abs(_) <= 0.2:
                return k
        else:
            return None

    @classmethod
    def transPNG(cls, im: Image) -> Image:
        """去除图片背景"""
        img = im.convert("RGBA")
        datas = img.getdata()
        newData = list()
        for item in datas:
            if item[0] > 220 and item[1] > 220 and item[2] > 220:
                newData.append((255, 255, 255, 0))
            else:
                newData.append(item)
        img.putdata(newData)
        # 抗锯齿
        img = img.resize(img.size, Image.Resampling.LANCZOS)
        return img

    @staticmethod
    def cut_bins(arr: np.ndarray, sep: int) -> Tuple[list, list]:
        """ 确保索引sorted 按照固定差值分bin"""
        length = len(arr)
        arr_index = list(range(length))
        index = []
        result = []
        # arr length is 1
        start = arr[0]
        if length == 1:
            result.append(arr[:])
        else:
            i = 0
            j = 0
            # 如果索引未到倒数第一个元素，就继续前进
            while i < length - 1:
                if abs(arr[i + 1] - start) <= sep:
                    i += 1
                else:
                    result.append(arr[j:i + 1])
                    index.append(arr_index[j:i + 1])
                    j = i + 1
                    i += 1
                    start = arr[j]
            # i为最后一个元素
            result.append(arr[j:i + 1])
            index.append(arr_index[j:i + 1])
        return index, result

    @staticmethod
    def center_point_generator(min_mz, max_mz, img_width, interval_width, size) -> np.ndarray:
        center = (min_mz + max_mz) / 2
        _ = img_width + interval_width
        x = []
        n = int(np.floor(size / 2))
        if size % 2 == 0:
            for i in range(0, n):
                distance = (0.5 + i) * _
                x.append(center + distance)
                x.append(center - distance)
        else:
            x.append(center)
            for i in range(1, n + 1):
                distance = i * _
                x.append(center + distance)
                x.append(center - distance)
        _x = np.asarray(x)
        return _x[np.argsort(_x)]

    @staticmethod
    def filename_handler(filename: str, escape: bool = True, full: bool = False) -> str:
        """
        :param filename: 文件名
        :param escape: 是否转义特殊字符
        :param full: 不缩减长度
        :return: new name
        """
        _name = "".join(re.findall(r'[^\*"/:?\\|<>]', filename, re.S)) if escape else filename
        if full:
            if len(_name) >= 200: # max 255
                return _name[:75] + 3*"." + _name[81:91] + 3*"." + _name[-75:]
            return _name
        else:
            # title 用缩写
            if _name.__contains__("_M"):
                p = "_M" + _name.split("_M", maxsplit=1)[1]
                _name = _name.split("_M", maxsplit=1)[0]
                if len(_name) > 40:
                    # 13 ~ _Product_M2-1
                    _name = _name[:30] + ".."
                _name = _name + p
            else:
                if len(_name) > 40:
                    _name = _name[:37] + "..."
            return _name


    @staticmethod
    def make_stems(mzs: np.ndarray, intensities: np.ndarray):
        """calculate where the stems of the spectrum peaks are going to be"""
        x = np.zeros([2, mzs.size], dtype="float")
        y = np.zeros(x.shape)
        x[:, :] = np.tile(mzs, (2, 1))
        y[1, :] = intensities
        return x, y

    @staticmethod
    def selectIndex(df: pd.DataFrame, col: str) -> list:
        """
        select best index
        :param df: sorted & reset_index
        :param col: col name
        :return: index selected
        """
        index_arr = []

        tmp_index = 0
        key = df.loc[0, col]
        _score, _f_score = df.loc[0, "Score"], df.loc[0, "Fragmentation Score"]

        for i in range(1, df.shape[0]):
            compound = df.loc[i, col]
            if key == compound:
                score, f_score = df.loc[i, "Score"], df.loc[i, "Fragmentation Score"]
                if score > _score or (score == _score and f_score > _f_score):
                    tmp_index = i
                    _score = score
                    _f_score = f_score
            else:
                index_arr.append(tmp_index)
                tmp_index = i
                _score, _f_score = df.loc[i, "Score"], df.loc[i, "Fragmentation Score"]
                key = compound
        index_arr.append(tmp_index)
        return index_arr

    @staticmethod
    def selectIndexOfBestName(df: pd.DataFrame, col: str) -> list:
        """
        根据名称去重，保留最佳数据库来源
        ['TCM','ANIMAL','HERB','PRODUCT']
        select best index
        :param df: sorted[按照 Description] & reset_index
        :param col: col name = Description
        :return: index selected
        """
        index_arr = []

        tmp_index = 0
        key = df.loc[0, col]
        _db = df.loc[0, "DB"]
        for i in range(1, df.shape[0]):
            compound = df.loc[i, col]
            if key == compound:
                db = df.loc[i, "DB"]
                if Report.database_priority.index(db) < Report.database_priority.index(_db):
                    tmp_index = i
                    _db = db
            else:
                index_arr.append(tmp_index)
                tmp_index = i
                _db = df.loc[i, "DB"]
                key = compound
        index_arr.append(tmp_index)
        return index_arr


    @staticmethod
    def msp_parser(msp: str) -> dict:
        with open(msp, "r", encoding="utf8") as f:
            metadata = {}
            mz = []
            intensity = []
            for line in f:
                line = line.rstrip()
                if len(line) == 0:
                    yield MS2spec(mz, intensity, metadata)
                    metadata = {}
                    mz = []
                    intensity = []
                    continue
                if line.__contains__(": "):
                    k, v = line.split(": ")
                    metadata[k] = v
                else:
                    m, _ = line.split(" ", maxsplit=1)
                    mz.append(float(m))
                    intensity.append(float(_))

    @staticmethod
    def spec_normalize(spec: MSSpectrum) -> MSSpectrum:
        "sort by mz ASC & normalize intensities to max 100"
        _spec = copy.deepcopy(spec)
        mzs, intensities = _spec.get_peaks()
        _mzs, _intensities = mzs[np.argsort(mzs)], intensities[np.argsort(mzs)]
        _intensities = _intensities * (100 / _intensities.max())
        _spec.set_peaks((_mzs, _intensities))
        return _spec

    def is_cpd_valid(self, inchikey, is_product=False):
        with self.ENGINE.connect() as conn:
            # 从标准品库查
            sql = f"""select * from herb_network where network_id={repr(inchikey)}""" if is_product else f"""select * from herb_annotation where inchikey={repr(inchikey)}"""
            fetch_data = pd.read_sql(text(sql), conn)
            if fetch_data.shape[0] == 0:
                # 从大库查
                _sql = f"""SELECT
                                a.exactMass AS exact_mass, 
                                b.smiles AS smiles, 
                                b.inchi AS inchi
                            FROM
                                compound_property AS a
                                INNER JOIN
                                compound_structure AS b
                                ON 
                                    a.cid = b.cid
                            WHERE
                                b.inchikey = {repr(inchikey)}"""
                fetch_data = pd.read_sql(text(_sql), conn)
            # 临时搜库项目
            # if fetch_data.shape[0] == 0:
            #     local_db = "/data/hstore1/database/database/tcm/localdb/local_DZLM2023050885.xlsx" if sys.platform == "linux" \
            #         else r"D:\pycharm_project\data_processing_pipeline\CM_DB\external\local_DZLM2023050885.xlsx"
            #     db = pd.read_excel(local_db, sheet_name=0)
            #     fetch_data = db.loc[db["inchikey"] == inchikey, :].reset_index(drop=True)
            if fetch_data.shape[0] == 0:
                return False
            return True

    @staticmethod
    def get_smiles(spec: MSSpectrum, engine: Engine, formula: str, inchikey: str, cpd_name: str,
                   ion_mode: str, adducts: list, is_product=False) -> Union[None, list]:
        """
        get fragment smiles using MetFrag
        :param spec: Sort by mz & normalize MSSpectrum
        :param engine: from sqlalchemy.engine import Engine
        :param formula: formula
        :param inchikey: inchikey
        :param cpd_name: cpd_name
        :param ion_mode: ion_mode
        :param adducts: adducts
        :return: None or [smiles],parent_smiles
        """
        #src_path = Path(__file__).parent / "data"
        src_path = Path("/data/hstore1/database/database/tcm")
        ADDUCT = {"M+H": 1,
                  "M+Na": 23,
                  "M+NH4": 18,
                  "M+K": 39,
                  "M-H": -1,
                  "M+FA-H": 59
                  }
        _mz, _intensity = spec.get_peaks()
        conn = engine.connect()
        for adduct in adducts:
            index = ADDUCT.get(adduct, 1) if ion_mode == "POS" else ADDUCT.get(adduct, -1)
            # prepare metFrag input
            with tempfile.TemporaryDirectory() as folder:
                task_id = uuid.uuid4()
                tmp_db_file_path = Path(folder) / f"{task_id}_localdb.csv"
                tmp_peak_file_path = Path(folder) / f"{task_id}_MS2.txt"
                tmp_config_file_path = Path(folder) / f"{task_id}_parameter_file.txt"
                tmp_result_file_path = Path(folder)
                tmp_result_file_name = f"{task_id}_FragmentSmilesExpl"
                tmp_result_file = tmp_result_file_path / f"{task_id}_FragmentSmilesExpl.psv"

                with open(str(tmp_db_file_path), "w", encoding="utf-8") as d:
                    d.write(
                        f""""Identifier","MonoisotopicMass","MolecularFormula","SMILES","InChI","InChIKey1","InChIKey2","InChIKey3","Name","InChIKey"\n""")
                    # 从标准品库查
                    # sql = f"""select * from herb_product where `NAME`={repr(inchikey)}""" if is_product else f"""select * from herb_msms where inchikey={repr(inchikey)}"""
                    sql = f"""SELECT MonoisotopicMass AS exact_mass, SMILES AS smiles, InChIKey AS inchikey, InChI AS inchi FROM herb_network where network_id={repr(inchikey)}""" if is_product else f"""select * from herb_annotation where inchikey={repr(inchikey)}"""
                    res = pd.read_sql(text(sql), conn)
                    # Change: 从大库查，标准品库有限
                    _sql = f"""SELECT
                                        a.computedExactMass AS exact_mass, 
                                        b.smiles AS smiles, 
                                        b.inchi AS inchi,
                                        b.inchikey AS inchikey
                                    FROM
                                        compound_computed_property AS a
                                        INNER JOIN
                                        compound_structure AS b
                                        ON 
                                            a.cid = b.cid
                                    WHERE
                                        b.inchikey = {repr(inchikey)}"""
                    _res = pd.read_sql(text(_sql), conn)
                    # 临时搜库项目
                    # local_db = "/data/hstore1/database/database/tcm/localdb/local_DZLM2023050885.xlsx"
                    # db = pd.read_excel(local_db, sheet_name=0)
                    # _res1 = db.loc[db["inchikey"] == inchikey, :].reset_index(drop=True)
                    #
                    fetch_data = res if res.shape[0] else _res  # if _res.shape[0] else _res1
                    if fetch_data.shape[0] == 0:
                        conn.close()
                        return None
                    smiles = fetch_data.loc[0, "smiles"]
                    inchi = fetch_data.loc[0, "inchi"]
                    exact_mass = fetch_data.loc[0, "exact_mass"]
                    # NEW FOR PRODUCT
                    inchikey = fetch_data.loc[0, "inchikey"]
                    inchikey1, inchikey2, inchikey3 = inchikey.split("-")
                    d.write(
                        f""""{inchikey}","{exact_mass}","{formula}","{smiles}","{inchi}","{inchikey1}","{inchikey2}","{inchikey3}","{cpd_name}","{inchikey}"\n""")
                with open(str(tmp_peak_file_path), "w", encoding="utf-8") as f:
                    for mz, ints in zip(_mz, _intensity):
                        f.write(f"{mz}\t{ints}\n")
                with open(str(src_path / "metfrag" / "parameter_template.txt"), "r", encoding="utf-8") as p:
                    content = p.read().replace("NeutralPrecursorMolecularFormula = C10H10O3",
                                               f"NeutralPrecursorMolecularFormula = {formula}"). \
                        replace("PeakListPath = MS2.txt", f"PeakListPath = {str(tmp_peak_file_path)}"). \
                        replace("LocalDatabasePath = localdb.csv",
                                f"LocalDatabasePath = {str(tmp_db_file_path)}") \
                        .replace("ResultsPath = .", f"ResultsPath = {str(tmp_result_file_path)}") \
                        .replace("SampleName = FragmentSmilesExpl", f"SampleName = {str(tmp_result_file_name)}")
                    if ion_mode == "NEG":
                        content = content.replace("IsPositiveIonMode = True", f"IsPositiveIonMode = False").replace(
                            "PrecursorIonMode = 1", f"PrecursorIonMode = {index}")
                    else:
                        content = content.replace("PrecursorIonMode = 1", f"PrecursorIonMode = {index}")
                    with open(str(tmp_config_file_path), "w", encoding="utf-8") as ff:
                        ff.write(content)
                # use subprocess instead
                r = sp.run(['java', '-jar', '/root/MetFrag2.4.5-CL.jar',
                            f"{str(tmp_config_file_path)}"],
                           stdout=sp.PIPE, stderr=sp.STDOUT)
                if r.returncode == 0:
                    with open(str(tmp_result_file), "r", encoding="utf-8") as file:
                        content = file.read().strip()
                    if content:
                        df = pd.read_csv(str(tmp_result_file), sep="|", header=0,
                                         index_col=False)
                        smilesOfExplPeaks = df.loc[0, "SmilesOfExplPeaks"]
                        smilesOfExplPeaks_dict = dict()
                        if isinstance(smilesOfExplPeaks, str) and smilesOfExplPeaks:
                            for kv in smilesOfExplPeaks.split(";"):
                                m, smiles = kv.split(":")
                                mz = float(m)
                                smilesOfExplPeaks_dict[mz] = smiles
                            _smiles = []
                            for z in _mz:
                                smi = smilesOfExplPeaks_dict.get(z, "")
                                _smiles.append(smi)
                            conn.close()
                            return _smiles
                # else:
                #     err = r.stderr.decode()
                #     logger.info("METFRAG INTERNAL ERROR".format(err))
        conn.close()
        return None

        #     with open(str(src_path / "metfrag" / "localdb.csv"), "w", encoding="utf-8") as d:
        #         d.write(
        #             f""""Identifier","MonoisotopicMass","MolecularFormula","SMILES","InChI","InChIKey1","InChIKey2","InChIKey3","Name","InChIKey"\n""")
        #         # 从标准品库查
        #         # sql = f"""select * from herb_product where `NAME`={repr(inchikey)}""" if is_product else f"""select * from herb_msms where inchikey={repr(inchikey)}"""
        #         sql = f"""SELECT MonoisotopicMass AS exact_mass, SMILES AS smiles, InChIKey AS inchikey, InChI AS inchi FROM herb_network where network_id={repr(inchikey)}""" if is_product else f"""select * from herb_annotation where inchikey={repr(inchikey)}"""
        #         res = pd.read_sql(text(sql), conn)
        #         # Change: 从大库查，标准品库有限
        #         _sql = f"""SELECT
        #                         a.computedExactMass AS exact_mass,
        #                         b.smiles AS smiles,
        #                         b.inchi AS inchi,
        #                         b.inchikey AS inchikey
        #                     FROM
        #                         compound_computed_property AS a
        #                         INNER JOIN
        #                         compound_structure AS b
        #                         ON
        #                             a.cid = b.cid
        #                     WHERE
        #                         b.inchikey = {repr(inchikey)}"""
        #         _res = pd.read_sql(text(_sql), conn)
        #         # 临时搜库项目
        #         # local_db = "/data/hstore1/database/database/tcm/localdb/local_DZLM2023050885.xlsx"
        #         # db = pd.read_excel(local_db, sheet_name=0)
        #         # _res1 = db.loc[db["inchikey"] == inchikey, :].reset_index(drop=True)
        #         #
        #         fetch_data = res if res.shape[0] else _res # if _res.shape[0] else _res1
        #         if fetch_data.shape[0] == 0:
        #             conn.close()
        #             return None
        #         smiles = fetch_data.loc[0, "smiles"]
        #         inchi = fetch_data.loc[0, "inchi"]
        #         exact_mass = fetch_data.loc[0, "exact_mass"]
        #         # NEW FOR PRODUCT
        #         inchikey = fetch_data.loc[0, "inchikey"]
        #         inchikey1, inchikey2, inchikey3 = inchikey.split("-")
        #         d.write(
        #             f""""{inchikey}","{exact_mass}","{formula}","{smiles}","{inchi}","{inchikey1}","{inchikey2}","{inchikey3}","{cpd_name}","{inchikey}"\n""")
        #     with open(str(src_path / "metfrag" / "MS2.txt"), "w", encoding="utf-8") as f:
        #         for mz, ints in zip(_mz, _intensity):
        #             f.write(f"{mz}\t{ints}\n")
        #     with open(str(src_path / "metfrag" / "parameter_template.txt"), "r", encoding="utf-8") as p:
        #         content = p.read().replace("NeutralPrecursorMolecularFormula = C10H10O3",
        #                                    f"NeutralPrecursorMolecularFormula = {formula}"). \
        #             replace("PeakListPath = MS2.txt", f"PeakListPath = {str(src_path / 'metfrag' / 'MS2.txt')}"). \
        #             replace("LocalDatabasePath = localdb.csv",
        #                     f"LocalDatabasePath = {str(src_path / 'metfrag' / 'localdb.csv')}") \
        #             .replace("ResultsPath = .", f"ResultsPath = {str(src_path / 'metfrag')}")
        #         if ion_mode == "NEG":
        #             content = content.replace("IsPositiveIonMode = True", f"IsPositiveIonMode = False").replace(
        #                 "PrecursorIonMode = 1", f"PrecursorIonMode = {index}")
        #         else:
        #             content = content.replace("PrecursorIonMode = 1", f"PrecursorIonMode = {index}")
        #         with open(str(src_path / "metfrag" / "parameter_file.txt"), "w", encoding="utf-8") as ff:
        #             ff.write(content)
        #
        #     # run MetFrag in cmd using java >= 8
        #     # status = os.system(
        #     #     f"java -jar /root/MetFrag2.4.5-CL.jar {str(src_path / 'metfrag' / 'parameter_file.txt')}")
        #     # use subprocess instead
        #     r = sp.run(['java', '-jar', '/root/MetFrag2.4.5-CL.jar', f"{str(src_path / 'metfrag' / 'parameter_file.txt')}"],
        #                     stdout=sp.PIPE, stderr=sp.STDOUT)
        #     if r.returncode == 0:
        #         with open(str(src_path / "metfrag" / "FragmentSmilesExpl.psv"), "r", encoding="utf-8") as file:
        #             content = file.read().strip()
        #         if content:
        #             df = pd.read_csv(str(src_path / "metfrag" / "FragmentSmilesExpl.psv"), sep="|", header=0,
        #                              index_col=False)
        #             smilesOfExplPeaks = df.loc[0, "SmilesOfExplPeaks"]
        #             smilesOfExplPeaks_dict = dict()
        #             if isinstance(smilesOfExplPeaks, str) and smilesOfExplPeaks:
        #                 for kv in smilesOfExplPeaks.split(";"):
        #                     m, smiles = kv.split(":")
        #                     mz = float(m)
        #                     smilesOfExplPeaks_dict[mz] = smiles
        #                 _smiles = []
        #                 for z in _mz:
        #                     smi = smilesOfExplPeaks_dict.get(z, "")
        #                     _smiles.append(smi)
        #                 conn.close()
        #                 return _smiles
        # conn.close()
        # return None

    def get_parent_smiles(cls, inchikey: str, is_product=False) -> Union[None, str]:
        engine = cls.ENGINE
        conn = engine.connect()
        # 从标准品库查
        sql = f"""select * from herb_network where network_id={repr(inchikey)}""" if is_product else f"""select * from herb_annotation where inchikey={repr(inchikey)}"""
        fetch_data = pd.read_sql(text(sql), conn)
        if fetch_data.shape[0] == 0:
            # 从大库查
            _sql = f"""SELECT
                            a.exactMass AS exact_mass, 
                            b.smiles AS smiles, 
                            b.inchi AS inchi
                        FROM
                            compound_property AS a
                            INNER JOIN
                            compound_structure AS b
                            ON 
                                a.cid = b.cid
                        WHERE
                            b.inchikey = {repr(inchikey)}"""
            fetch_data = pd.read_sql(text(_sql), conn)
        # if fetch_data.shape[0] == 0:
        #     # 临时搜库项目
        #     local_db = "/data/hstore1/database/database/tcm/localdb/local_DZLM2023050885.xlsx" if sys.platform == "linux" \
        #         else r"D:\pycharm_project\data_processing_pipeline\CM_DB\external\local_DZLM2023050885.xlsx"
        #     db = pd.read_excel(local_db, sheet_name=0)
        #     fetch_data = db.loc[db["inchikey"] == inchikey, :].reset_index(drop=True)
        if fetch_data.shape[0] == 0:
            return None
        if "SMILES" in fetch_data.columns.values:
            fetch_data.rename(columns={"SMILES": "smiles"}, inplace=True)
        smiles = fetch_data.loc[0, "smiles"]
        conn.close()
        return smiles

    @staticmethod
    def get_fragment_ions(exp_spec: MSSpectrum):
        mzs, intensities = exp_spec.get_peaks()
        # subset
        mzs = mzs[intensities >= 5]
        intensities = intensities[intensities >= 5]

        # sort
        _ = mzs[np.argsort(intensities)[::-1]]
        _ = _[np.argsort(_)]
        if _.size <= 10:
            return ", ".join([str(m) for m in _])
        return ", ".join([f"{m:.4f}" for m in _[-10:]])

    @staticmethod
    def plot_defect(ax: Axes, parent_smiles: Union[str, None], parentMz: float, **kwargs):
        img_size = 300
        if kwargs.get("img_size"):
            img_size = int(kwargs["img_size"])

        min_mz = max(0, np.floor(parentMz / 100 - 1) * 100)
        max_mz = np.floor(parentMz / 100 + 2) * 100
        anchor1 = min_mz + 75
        anchor2 = max_mz - 75
        ax.set_ylim([0, 100])
        ax.set_xlim([min_mz, max_mz])
        # ratio_xy = (max_mz-min_mz) / 100
        # fig_aspect = 16 / 100
        # aspect = ratio_xy * fig_aspect

        # plot smiles
        if parent_smiles:
            mol = Chem.MolFromSmiles(parent_smiles)
            formula = CalcMolFormula(mol)
            img = Draw.MolToImage(mol, size=(int(img_size * 8 / 5), img_size))
            img = Report.transPNG(img)
            ax.imshow(img, extent=(anchor1 - 25, anchor1 + 25, 40, 90), origin="upper", aspect=5 / 8, zorder=999)
            ax.text(anchor1, 40, f"Parent Compound: {formula}", fontsize=8, color="red", zorder=99,
                    horizontalalignment="center", verticalalignment='top')

        ax.text(anchor2, 55, "No MS/MS spectrum detected!", fontsize=8, color="red", bbox=dict(alpha=0.2, color="red"),
                horizontalalignment="center", verticalalignment='top')
        # axis
        ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.yaxis.set_major_locator(MultipleLocator(20))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        y_ticks = ax.get_yticks()
        ax.set_yticks(y_ticks[(y_ticks <= 100) & (y_ticks >= 0)])

        # global config
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlabel("m/z", style="italic")
        ax.set_ylabel("Relative Intensity")

    def plot_spectrum(self, ax: Axes, exp_spec: MSSpectrum, parent_smiles: str, parentMz: float, **kwargs):
        img_size = 300
        if kwargs.get("img_size"):
            img_size = int(kwargs["img_size"])
        # load data, coordinates
        mzs_exp, intensities_exp = exp_spec.get_peaks()
        x1, y1 = self.make_stems(mzs_exp, intensities_exp)

        # plot lines
        lines = ax.plot(x1, y1, color="teal", linewidth=1.0, zorder=5)

        # x limit
        mz_min = mzs_exp.min()
        mz_max = mzs_exp.max()
        min_mz = max(0, np.floor(mz_min / 100 - 1) * 100)
        max_mz = np.floor(mz_max / 100 + 1) * 100
        xlim = max_mz - min_mz
        ax.set_xlim(min_mz, max_mz)

        # y limit 150 + 100 = 250
        ax.set_ylim([0, 160])
        ylim = 160
        # aspect
        ratio_xy = xlim / ylim
        fig_aspect = 16*1.6/100
        aspect = ratio_xy * fig_aspect

        # imshow params
        img_width = int(np.floor(xlim / 7))
        img_height = 40
        ratio = 40 / img_width

        # mz范围分bin
        bin_width = xlim * 0.1
        bin_idx1, bin_grp1 = self.cut_bins(x1[1, :], bin_width)
        # 保留每个Group中最大强度的peak in exp
        select_index1 = []
        for index_arr in bin_idx1:
            intensities = y1[1, :][index_arr]
            idx_max = index_arr[np.argmax(intensities)]
            select_index1.append(idx_max)

        # plot texts
        adjust_y = 1
        text_plotted = []
        for i in zip(select_index1):
            if y1[1, :][i] >= 10:
                text_plotted.append(x1[1, :][i])
                ax.text(x1[1, :][i], y1[1, :][i] + adjust_y, f"{x1[1, :][i]:.4f}", fontsize=6, color="red", zorder=99,
                        horizontalalignment="center")
        # plot smiles
        # find parentMZ in 5 ppm
        _ = np.abs(mzs_exp - parentMz).argmin()
        parent_mz = mzs_exp[_]
        parent_intensity = intensities_exp[_]
        # calc ppm
        ppm = abs(parent_mz - parentMz) * 1e6 / parentMz
        # parentMZ was found
        if ppm <= 5:
            ax.text(parent_mz, parent_intensity + adjust_y, f"{parent_mz:.4f}", fontsize=6, color="red", zorder=99,
                    horizontalalignment="center")
            # # iter data for plot
            # ax_y = 116
            # # plot Parent
            # center = (min_mz + max_mz) / 2
            # mol = Chem.MolFromSmiles(parent_smiles, sanitize=False)
            # img = Draw.MolToImage(mol, size=(int(img_size / aspect / ratio), img_size))
            # img = Report.transPNG(img)
            # ax.imshow(img, extent=(center - img_width / 2, center + img_width / 2, ax_y, ax_y + img_height),
            #           origin="upper", aspect=aspect, zorder=999)
            # ax.arrow(parent_mz, parent_intensity, center - parent_mz, ax_y - parent_intensity, linewidth=0.2, alpha=0.1,
            #          color="grey")
            # if parent_mz not in text_plotted:
            #     ax.text(center, ax_y, f"{parent_mz:.4f}", fontsize=6, color="red", zorder=99,
            #             horizontalalignment="center", verticalalignment='top')
        # axis
        if xlim >= 1000:
            ax.xaxis.set_major_locator(MultipleLocator(200))
            ax.xaxis.set_minor_locator(MultipleLocator(50))
        else:
            ax.xaxis.set_major_locator(MultipleLocator(50))
            ax.xaxis.set_minor_locator(MultipleLocator(10))

        ax.yaxis.set_major_locator(MultipleLocator(20))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        y_ticks = ax.get_yticks()
        ax.set_yticks(y_ticks[(y_ticks <= 100) & (y_ticks >= 0)])
        # global config
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlabel("m/z", style="italic")
        ax.set_ylabel("Relative Intensity")
        ax.legend([lines[0], ], ['expSpec', ], fontsize=6, loc="best")
        ax.set_aspect(aspect)
        return self.get_fragment_ions(exp_spec)

    def plot_spectrum_with_smiles(self, ax: Axes, exp_spec: MSSpectrum, _smiles: list, parent_smiles: str,
                                  parentMz: float, **kwargs):
        # fragment_ions_str = None
        img_size = 300
        if kwargs.get("img_size"):
            img_size = int(kwargs["img_size"])
        # load data, coordinates
        mzs_exp, intensities_exp = exp_spec.get_peaks()
        x1, y1 = self.make_stems(mzs_exp, intensities_exp)

        # plot lines
        lines = ax.plot(x1, y1, color="teal", linewidth=1.0, zorder=5)

        # x limit
        mz_min = mzs_exp.min()
        mz_max = mzs_exp.max()
        min_mz = max(0, np.floor(mz_min / 100 - 1) * 100)
        max_mz = np.floor(mz_max / 100 + 1) * 100
        xlim = max_mz - min_mz
        ax.set_xlim(min_mz, max_mz)

        # y limit 150 + 100 = 250
        ax.set_ylim([0, 160])
        ylim = 160

        # aspect
        ratio_xy = xlim / ylim
        fig_aspect = 16*1.6/100
        aspect = ratio_xy * fig_aspect

        # imshow params
        # imshow_aspect<W/H> * ylim / xlim  = axes_aspect<H/W>
        img_width = int(np.floor(xlim / 7))
        img_height = 40
        interval_width = img_width / 5
        ratio = 40 / img_width
        # Draw.IMG 尺寸 [300 * 8 / 5, 300]

        # mz范围分bin
        bin_width = xlim * 0.1
        bin_idx1, bin_grp1 = self.cut_bins(x1[1, :], bin_width)
        # 保留每个Group中最大强度的peak in exp
        select_index1 = []
        for index_arr in bin_idx1:
            intensities = y1[1, :][index_arr]
            idx_max = index_arr[np.argmax(intensities)]
            select_index1.append(idx_max)

        # plot texts
        adjust_y = 1
        text_plotted = []
        for i in zip(select_index1):
            if y1[1, :][i] >= 10:
                text_plotted.append(x1[1, :][i])
                ax.text(x1[1, :][i], y1[1, :][i] + adjust_y, f"{x1[1, :][i]:.4f}", fontsize=6, color="red", zorder=99,
                        horizontalalignment="center")

        if _smiles:
            # plot smiles
            smiles = np.asarray(_smiles)
            # find parentMZ in 5 ppm
            _ = np.abs(mzs_exp - parentMz).argmin()
            parent_mz = mzs_exp[_]
            parent_intensity = intensities_exp[_]
            # calc ppm
            ppm = abs(parent_mz - parentMz) * 1e6 / parentMz
            # filter None and mz < 60
            smiles_arr = smiles[np.where((smiles != '') & (mzs_exp >= 60))]
            mz_arr = mzs_exp[np.where((smiles != '') & (mzs_exp >= 60))]
            intensity_arr = intensities_exp[np.where((smiles != '') & (mzs_exp >= 60))]
            #
            # filter smiles max 5
            if smiles_arr.size > 5:
                smiles_arr = smiles_arr[-5:]
                mz_arr = mz_arr[-5:]
                intensity_arr = intensity_arr[-5:]
            # parentMZ was found
            if ppm <= 5:
                ax.text(parent_mz, parent_intensity + adjust_y, f"{parent_mz:.4f}", fontsize=6, color="red", zorder=99,
                        horizontalalignment="center")
                # 不再绘制 `parent_smiles`
                # if not np.isin(parent_mz, mz_arr):
                #     mz_arr = np.append(mz_arr, parent_mz)
                #     intensity_arr = np.append(intensity_arr, parent_intensity)
                #     smiles_arr = np.append(smiles_arr, parent_smiles)
            # sort by mz
            mz_arr = mz_arr[np.argsort(mz_arr)]
            # fragment_ions
            # fragment_ions_str = "; ".join([f"{m:.4f}" for m in mz_arr])
            intensity_arr = intensity_arr[np.argsort(mz_arr)]
            smiles_arr = smiles_arr[np.argsort(mz_arr)]
            # iter data for plot
            count_img = smiles_arr.size
            ax_y = 116
            center_points = self.center_point_generator(min_mz, max_mz, img_width, interval_width, count_img)
            for _x, _y, smi, p in zip(mz_arr, intensity_arr, smiles_arr, center_points):
                # plot Fragments
                mol = Chem.MolFromSmiles(smi, sanitize=False)
                img = Draw.MolToImage(mol, size=(int(img_size / aspect / ratio), img_size))
                img = Report.transPNG(img)
                ax.imshow(img, extent=(p - img_width / 2, p + img_width / 2, ax_y, ax_y + img_height), origin="upper",
                          aspect=aspect, zorder=999)
                ax.arrow(_x, _y, p - _x, ax_y - _y, linewidth=0.2, alpha=0.1, color="grey")
                if _x not in text_plotted:
                    ax.text(p, ax_y, f"{_x:.4f}", fontsize=6, color="red", zorder=99,
                            horizontalalignment="center", verticalalignment='top')
        # axis
        if xlim >= 1000:
            ax.xaxis.set_major_locator(MultipleLocator(200))
            ax.xaxis.set_minor_locator(MultipleLocator(50))
        else:
            ax.xaxis.set_major_locator(MultipleLocator(50))
            ax.xaxis.set_minor_locator(MultipleLocator(10))

        ax.yaxis.set_major_locator(MultipleLocator(20))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        y_ticks = ax.get_yticks()
        ax.set_yticks(y_ticks[(y_ticks <= 100) & (y_ticks >= 0)])
        # global config
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlabel("m/z", style="italic")
        ax.set_ylabel("Relative Intensity")
        ax.legend([lines[0], ], ['expSpec', ], fontsize=6, loc="best")
        ax.set_aspect(aspect)
        return self.get_fragment_ions(exp_spec)

    def plot_mirror_spectrum_with_smiles_deprecated(self, ax: Axes, exp_spec: MSSpectrum, ref_spec: MSSpectrum,
                                                    _smiles: list, parent_smiles: str, parentMz: float, **kwargs):
        """
        This method is deprecated! plot Fragments in top!
        plot mirror spectrum with fragments
        :param ax: Axes object
        :param exp_spec: exp_spec
        :param ref_spec: ref_spec
        :param _smiles: smiles
        :param parent_smiles: parent_smiles
        :param parentMz: parentMz
        """
        img_size = 300
        if kwargs.get("img_size"):
            img_size = int(kwargs["img_size"])
        # load data, coordinates
        mzs_exp, intensities_exp = exp_spec.get_peaks()
        x1, y1 = self.make_stems(mzs_exp, intensities_exp)

        mzs_ref, intensities_ref = ref_spec.get_peaks()
        x2, y2 = self.make_stems(mzs_ref, intensities_ref)

        # plot lines
        ax.plot(x1, y1, color="teal", linewidth=1.0, zorder=5)
        ax.plot(x2, -y2, color="red", linewidth=1.0, zorder=5)

        # x limit
        mz_min = np.hstack((mzs_ref, mzs_exp)).min()
        mz_max = np.hstack((mzs_ref, mzs_exp)).max()
        min_mz = max(0, np.floor(mz_min / 100 - 1) * 100)
        max_mz = np.floor(mz_max / 100 + 1) * 100
        xlim = max_mz - min_mz
        ax.set_xlim(min_mz, max_mz)

        # y limit
        ax.set_ylim([-110, 155])
        ylim = 265

        # aspect
        ratio_xy = xlim / ylim
        fig_aspect = FIG_ASPECT
        aspect = ratio_xy * fig_aspect

        # imshow params
        # imshow_aspect<W/H> * ylim / xlim  = axes_aspect<H/W>
        img_width = int(np.floor(xlim / 7))
        interval_width = img_width / 5
        img_height = 40
        ratio = 40 / img_width
        # Draw.IMG 尺寸 [300 * 8 / 5, 300]

        # mz范围分bin
        bin_width = xlim * 0.1
        bin_idx1, bin_grp1 = self.cut_bins(x1[1, :], bin_width)
        # 保留每个Group中最大强度的peak in exp
        select_index1 = []
        for index_arr in bin_idx1:
            intensities = y1[1, :][index_arr]
            idx_max = index_arr[np.argmax(intensities)]
            select_index1.append(idx_max)

        bin_idx2, bin_grp2 = self.cut_bins(x2[1, :], bin_width)
        # 保留每个Group中最大强度的peak in ref
        select_index2 = []
        for index_arr in bin_idx2:
            intensities = y2[1, :][index_arr]
            idx_max = index_arr[np.argmax(intensities)]
            select_index2.append(idx_max)
        # plot texts
        adjust_y = 1
        text_plotted = []
        for i, j in zip(select_index1, select_index2):
            if y1[1, :][i] >= 10:
                text_plotted.append(x1[1, :][i])
                ax.text(x1[1, :][i], y1[1, :][i] + adjust_y, f"{x1[1, :][i]:.4f}", fontsize=6, color="red",
                        horizontalalignment="center")
            if y2[1, :][j] >= 10:
                ax.text(x2[1, :][j], (-y2)[1, :][j] - adjust_y, f"{x2[1, :][j]:.4f}", fontsize=6, color="red",
                        horizontalalignment="center", verticalalignment='top')

        # plot smiles
        smiles = np.asarray(_smiles)
        # find parentMZ in 5 ppm
        _ = np.abs(mzs_exp - parentMz).argmin()
        parent_mz = mzs_exp[_]
        parent_intensity = intensities_exp[_]
        # calc ppm
        ppm = abs(parent_mz - parentMz) * 1e6 / parentMz
        # filter None and mz < 60
        smiles_arr = smiles[np.where((smiles != '') & (mzs_exp >= 60))]
        mz_arr = mzs_exp[np.where((smiles != '') & (mzs_exp >= 60))]
        intensity_arr = intensities_exp[np.where((smiles != '') & (mzs_exp >= 60))]
        #
        # filter smiles max 5
        if smiles_arr.size > 5:
            smiles_arr = smiles_arr[-5:]
            mz_arr = mz_arr[-5:]
            intensity_arr = intensity_arr[-5:]
        # parentMZ was found
        if ppm <= 5:
            ax.text(parent_mz, parent_intensity + adjust_y, f"{parent_mz:.4f}", fontsize=6, color="red",
                    horizontalalignment="center")
            # 不再绘制 `parent_smiles`
            # if not np.isin(parent_mz, mz_arr):
            #     mz_arr = np.append(mz_arr, parent_mz)
            #     intensity_arr = np.append(intensity_arr, parent_intensity)
            #     smiles_arr = np.append(smiles_arr, parent_smiles)
        # sort by mz
        mz_arr = mz_arr[np.argsort(mz_arr)]
        intensity_arr = intensity_arr[np.argsort(mz_arr)]
        smiles_arr = smiles_arr[np.argsort(mz_arr)]
        # iter data for plot
        count_img = smiles_arr.size
        ax_y = 113
        center_points = self.center_point_generator(min_mz, max_mz, img_width, interval_width, count_img)
        for _x, _y, smi, p in zip(mz_arr, intensity_arr, smiles_arr, center_points):
            # plot Fragments
            mol = Chem.MolFromSmiles(smi, sanitize=False)
            img = Draw.MolToImage(mol, size=(int(img_size / aspect / ratio), img_size))
            img = Report.transPNG(img)
            ax.imshow(img, extent=(p - img_width / 2, p + img_width / 2, ax_y, ax_y + img_height), origin="upper",
                      aspect=aspect, zorder=999)
            ax.arrow(_x, _y, p - _x, ax_y - _y, linewidth=0.3, alpha=0.1)
            if _x not in text_plotted:
                ax.text(p, ax_y, f"{_x:.4f}", fontsize=6, color="red",
                        horizontalalignment="center", verticalalignment='top')

        # axis
        ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.yaxis.set_major_locator(MultipleLocator(20))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        y_ticks = ax.get_yticks()
        ax.set_yticks(y_ticks[(y_ticks <= 100) & (y_ticks >= -100)])
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{abs(int(x))}"))

        # global config
        ax.axhline(0, color="#9E9E9E", zorder=10, lw=0.8)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlabel("m/z", style="italic")
        ax.set_ylabel("Relative Intensity")

    def plot_mirror_spectrum_with_smiles(self, ax: Axes, exp_spec: MSSpectrum, ref_spec: MSSpectrum,
                                         _smiles: list, parent_smiles: str, parentMz: float, **kwargs):
        """
        New features:
        ____________
        plot mirror spectrum with fragments in bottom.
        :param ax: Axes object
        :param exp_spec: exp_spec
        :param ref_spec: ref_spec
        :param _smiles: smiles
        :param parent_smiles: parent_smiles
        :param parentMz: parentMz
        """
        # fragment_ions_str = None

        img_size = 300
        if kwargs.get("img_size"):
            img_size = int(kwargs["img_size"])

        # load data, coordinates
        mzs_exp, intensities_exp = exp_spec.get_peaks()
        x1, y1 = self.make_stems(mzs_exp, intensities_exp)

        mzs_ref, intensities_ref = ref_spec.get_peaks()
        x2, y2 = self.make_stems(mzs_ref, intensities_ref)

        # plot lines
        tops = ax.plot(x1, y1, color="teal", linewidth=1.0, zorder=5)
        bots = ax.plot(x2, -y2, color="red", linewidth=1.0, zorder=5)

        # x limit
        mz_min = np.hstack((mzs_ref, mzs_exp)).min()
        mz_max = np.hstack((mzs_ref, mzs_exp)).max()
        min_mz = max(0, np.floor(mz_min / 100 - 1) * 100)
        max_mz = np.floor(mz_max / 100 + 1) * 100
        xlim = max_mz - min_mz
        ax.set_xlim(min_mz, max_mz)

        # y limit
        ax.set_ylim([-160, 110])
        ylim = 270

        # aspect
        ratio_xy = xlim / ylim
        fig_aspect = FIG_ASPECT
        aspect = ratio_xy * fig_aspect

        # imshow params
        # imshow_aspect<W/H> * ylim / xlim  = axes_aspect<H/W>
        img_width = int(np.floor(xlim / 7))
        interval_width = img_width / 5
        img_height = 40
        ratio = 40 / img_width
        # Draw.IMG 尺寸 [300 * 8 / 5, 300]

        # mz范围分bin
        bin_width = xlim * 0.1
        bin_idx1, bin_grp1 = self.cut_bins(x1[1, :], bin_width)
        # 保留每个Group中最大强度的peak in exp
        select_index1 = []
        for index_arr in bin_idx1:
            intensities = y1[1, :][index_arr]
            idx_max = index_arr[np.argmax(intensities)]
            select_index1.append(idx_max)

        bin_idx2, bin_grp2 = self.cut_bins(x2[1, :], bin_width)
        # 保留每个Group中最大强度的peak in ref
        select_index2 = []
        for index_arr in bin_idx2:
            intensities = y2[1, :][index_arr]
            idx_max = index_arr[np.argmax(intensities)]
            select_index2.append(idx_max)
        # plot texts
        adjust_y = 1
        text_plotted = []
        for i, j in zip(select_index1, select_index2):
            if y1[1, :][i] >= 10:
                ax.text(x1[1, :][i], y1[1, :][i] + adjust_y, f"{x1[1, :][i]:.4f}", fontsize=6, color="red", zorder=99,
                        horizontalalignment="center")
            if y2[1, :][j] >= 10:
                text_plotted.append(x2[1, :][j])
                ax.text(x2[1, :][j], (-y2)[1, :][j] - adjust_y, f"{x2[1, :][j]:.4f}", fontsize=6, color="red",
                        zorder=99,
                        horizontalalignment="center", verticalalignment='top')
        # new
        _ref = np.abs(mzs_ref - parentMz).argmin()
        parent_mz_ref = mzs_ref[_ref]
        parent_intensity_ref = intensities_ref[_ref]
        ppm_ref = abs(parent_mz_ref - parentMz) * 1e6 / parentMz

        # find parentMZ in 5 ppm EXP
        _ = np.abs(mzs_exp - parentMz).argmin()
        parent_mz = mzs_exp[_]
        parent_intensity = intensities_exp[_]
        # calc ppm
        ppm = abs(parent_mz - parentMz) * 1e6 / parentMz
        # parentMZ was found
        if ppm <= 5:
            ax.text(parent_mz, parent_intensity + adjust_y, f"{parent_mz:.4f}", fontsize=6, color="red", zorder=99,
                    horizontalalignment="center")
        # plot smiles
        if _smiles:
            smiles = np.asarray(_smiles)
            # filter None and mz < 60
            smiles_arr = smiles[np.where((smiles != '') & (mzs_ref >= 60))]
            mz_arr = mzs_ref[np.where((smiles != '') & (mzs_ref >= 60))]
            intensity_arr = intensities_ref[np.where((smiles != '') & (mzs_ref >= 60))]

            # filter smiles max 5
            if smiles_arr.size > 5:
                smiles_arr = smiles_arr[-5:]
                mz_arr = mz_arr[-5:]
                intensity_arr = intensity_arr[-5:]
            # 不再绘制 `parent_smiles`
            # if ppm_ref <= 5:
            #     if not np.isin(parent_mz_ref, mz_arr):
            #         mz_arr = np.append(mz_arr, parent_mz_ref)
            #         intensity_arr = np.append(intensity_arr, parent_intensity_ref)
            #         smiles_arr = np.append(smiles_arr, parent_smiles)
            # sort by mz
            mz_arr = mz_arr[np.argsort(mz_arr)]
            # fragment_ions
            # fragment_ions_str = "; ".join([f"{m:.4f}" for m in mz_arr])

            intensity_arr = intensity_arr[np.argsort(mz_arr)]
            smiles_arr = smiles_arr[np.argsort(mz_arr)]
            # iter data for plot
            count_img = smiles_arr.size
            ax_y = -116
            center_points = self.center_point_generator(min_mz, max_mz, img_width, interval_width, count_img)
            for _x, _y, smi, p in zip(mz_arr, -intensity_arr, smiles_arr, center_points):
                # plot Fragments
                mol = Chem.MolFromSmiles(smi, sanitize=False)
                img = Draw.MolToImage(mol, size=(int(img_size / aspect / ratio), img_size))
                img = Report.transPNG(img)
                ax.imshow(img, extent=(p - img_width / 2, p + img_width / 2, ax_y, ax_y - img_height), origin="lower",
                          aspect=aspect, zorder=999)
                ax.arrow(_x, _y, p - _x, ax_y - _y, linewidth=0.2, alpha=0.1, color="grey")
                if _x not in text_plotted:
                    ax.text(p, ax_y, f"{_x:.4f}", fontsize=6, color="red", zorder=99,
                            horizontalalignment="center", verticalalignment='bottom')
        # else:
        #     if ppm_ref <= 5:
        #         fragment_ions_str = str(parent_mz_ref)
        #         mol = Chem.MolFromSmiles(parent_smiles, sanitize=False)
        #         img = Draw.MolToImage(mol, size=(int(img_size / aspect / ratio), img_size))
        #         img = Report.transPNG(img)
        #         ax_y = -116
        #         p = (max_mz + min_mz) / 2
        #         ax.imshow(img, extent=(p - img_width / 2, p + img_width / 2, ax_y, ax_y - img_height), origin="lower",
        #                   aspect=aspect, zorder=999)
        #         ax.arrow(parent_mz_ref, -parent_intensity_ref, p - parent_mz_ref, ax_y + parent_intensity_ref,
        #                  linewidth=0.2, alpha=0.1, color="grey")
        #         if parent_mz_ref not in text_plotted:
        #             ax.text(p, ax_y, f"{parent_mz_ref:.4f}", fontsize=6, color="red", zorder=99,
        #                     horizontalalignment="center", verticalalignment='bottom')
        # axis
        if xlim >= 1000:
            ax.xaxis.set_major_locator(MultipleLocator(200))
            ax.xaxis.set_minor_locator(MultipleLocator(50))
        else:
            ax.xaxis.set_major_locator(MultipleLocator(50))
            ax.xaxis.set_minor_locator(MultipleLocator(10))

        ax.yaxis.set_major_locator(MultipleLocator(20))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        y_ticks = ax.get_yticks()
        ax.set_yticks(y_ticks[(y_ticks <= 100) & (y_ticks >= -100)])
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{abs(int(x))}"))

        # global config
        ax.axhline(0, color="#9E9E9E", zorder=10, lw=0.8)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlabel("m/z", style="italic")
        ax.set_ylabel("Relative Intensity")
        # new
        ax.legend([tops[0], bots[0]], ['expSpec', 'refSpec'], fontsize=6, loc="best")
        return self.get_fragment_ions(exp_spec)

    @functools.lru_cache(maxsize=1000)
    def parse_meta(self) -> Tuple[list, pd.DataFrame]:
        """parse sample info"""
        mode_used = []

        SAMPLE_INFO = Path(self.data) / "sampleinfo.xlsx"
        sample_info = pd.read_excel(SAMPLE_INFO, sheet_name="样品信息", index_col=0)

        is_null = sample_info.isnull().any(axis=1)
        if not is_null["POS"]:
            mode_used.append("POS")
        if not is_null["NEG"]:
            mode_used.append("NEG")

        return mode_used, sample_info

    @functools.lru_cache(maxsize=1000)
    def parse_project(self) -> Tuple[str, pd.DataFrame]:
        """parse project info"""
        PROJECT_INFO = Path(self.data) / "sampleinfo.xlsx"
        project_info = pd.read_excel(PROJECT_INFO, sheet_name="项目信息").applymap(lambda x: x.strip() if isinstance(x, str) else x)
        # 使用线上登记单
        if not self.local:
            if self.form is None:
                aid: str
                obj = re.search('-a[a-z]', project_info.loc[0, "项目编号"])
                if isinstance(obj, re.Match):
                    aid = project_info.loc[0, "项目编号"]
                else:
                    aid = project_info.loc[0, "项目编号"] + "-aa"
                self.form = self.get_online_register_form(aid)
                project_path = Path(self.data).cwd().name + "-" + self.form.get("客户名称", "***") + "-" + project_info.loc[0, "项目名称"] + "结题报告"
                r = sp.run(["tool_GetAnalystInfo", "-ad", aid, "-sp", str(Path(self.data) / project_path / "分析确认单"), "-fn", "分析确认单.xlsx", "--overwrite"])
                if r.returncode != 0:
                    raise SaveOnlineRegisterTableError("保存线上登记单出错!")
        else:
            self.form = dict()
        project_info.loc[0, "客户名称"] = self.form.get("客户名称", project_info.loc[0, "客户名称"])
        project_info.loc[0, "联系人名称"] = self.form.get("联系人", project_info.loc[0, "联系人名称"])
        project_info.loc[0, "样本"] = self.form.get("样本类型", project_info.loc[0, "样本"])
        project_path = Path(self.data).cwd().name + "-" + project_info.loc[0, "客户名称"] + "-" + project_info.loc[0, "项目名称"] + "结题报告"
        return project_path, project_info

    @functools.lru_cache(maxsize=1000)
    def parse_exp_col(self):
        """
        返回样本的实际列名
        """
        mode_used, sample_info = self.parse_meta()
        sample_info = sample_info.T
        mode = mode_used[0]
        tcm_cols = sample_info[mode].drop(index=["QI定性结果", "QI定量结果", "QI二级MSP"])\
            .apply(lambda x: x.rsplit(".", 1)[0]).tolist()
        return tcm_cols

    def fill_qi_df(self, data_processed, dedup=True):
        # concat res in two ion modes
        res = pd.concat(data_processed).sort_values("Compound").reset_index(drop=True)
        # add new columns
        res["InChIKey"] = ""
        res["SMILES"] = ""
        res["HMDB"] = ""
        res["METLIN"] = ""
        res["LipidMaps"] = ""
        res["KEGG"] = ""
        res["PubChem"] = ""
        res["CAS"] = ""
        res["中文名"] = ""
        res["中文大类"] = ""
        res["中文子类"] = ""
        res["英文大类"] = ""
        res["英文子类"] = ""
        res["厂家"] = ""
        res["货号"] = ""
        res["纯度"] = ""

        # assign values to new columns from database
        conn = self.ENGINE.connect()
        for i in range(res.shape[0]):
            # 1.identifiers
            cid = res.loc[i, "Compound ID"]
            # if isinstance(hmdbid, str) and hmdbid.startswith("HMDB"):
            #     sql = self.sql1 % hmdbid
            # else:
            sql = self.sql2 % cid
            res_proxy = conn.execute(text(sql)).fetchall()

            if len(res_proxy) != 0:
                inchikey, smiles, _, hmdbid, metlin, lipidmaps, kegg, pubchem, cas = res_proxy[0][0], res_proxy[0][1], \
                                                                                     res_proxy[0][2], \
                                                                                     res_proxy[0][3], res_proxy[0][4], \
                                                                                     res_proxy[0][5], \
                                                                                     res_proxy[0][6], res_proxy[0][7], \
                                                                                     res_proxy[0][8]
                res.loc[i, "InChIKey"] = inchikey
                res.loc[i, "SMILES"] = smiles
                res.loc[i, "HMDB"] = hmdbid
                res.loc[i, "METLIN"] = metlin
                res.loc[i, "LipidMaps"] = lipidmaps
                res.loc[i, "KEGG"] = kegg
                res.loc[i, "PubChem"] = pubchem
                res.loc[i, "CAS"] = cas
            # 2.cn class
            key = res.loc[i, "Compound ID"]
            sql = self.sql3 % key
            res_proxy = conn.execute(text(sql)).fetchall()
            # TODO: HERB库字段完善后添加来源和这些信息
            if len(res_proxy):
                res.loc[i, "中文名"] = res_proxy[0][0]
                res.loc[i, "中文大类"] = res_proxy[0][2]
                res.loc[i, "中文子类"] = res_proxy[0][3]
                res.loc[i, "厂家"] = res_proxy[0][4]
                res.loc[i, "货号"] = res_proxy[0][5]
                res.loc[i, "纯度"] = res_proxy[0][6]
                res.loc[i, "英文大类"] = res_proxy[0][7]
                res.loc[i, "英文子类"] = res_proxy[0][8]
        conn.close()

        # del duplicated cpd and select best
        if dedup:
            res = res.sort_values("Compound ID").reset_index(drop=True)
            best_index = self.selectIndex(res, "Compound ID")
            res = res.loc[best_index, :].sort_values("Compound").reset_index(drop=True)
        return res

    def raw2mzml(self):
        logger.info("convert rawdata to mzmldata, time-consuming step")
        cwd = Path(self.data).cwd()
        rawdata = cwd / "rawdata"
        if len(list(rawdata.rglob('*.raw'))) == 0:
            raise RawDataNotFoundError('请检查rawdata文件夹！')
        mzmldata = cwd / "mzmldata"
        if len(list(mzmldata.rglob('*.mzML'))) != 0:
            logger.info("detected mzML file in folder mzmldata, skip this step!")
            return
        for mode in ("POS", "NEG"):
            (mzmldata / mode).mkdir(parents=True, exist_ok=True)
            args = ["sh", MS_CONVERT, str(rawdata / mode), str(mzmldata / mode)]
            r = sp.run(args)
            if r.returncode != 0:
                raise ShellRunTimeError("mzmldata转换失败！")
        logger.info("convert rawdata to mzmldata successfully!")

    def raw2obs(self):
        if not self.local:
            logger.info("upload rawdata to HuaWei OBS server")
            cwd = Path(self.data).cwd()
            if (cwd / 'log.txt').exists():
                logger.info("log.txt found in path, skip upload!")
                return
            project_path, project_info = self.parse_project()
            server_path = "代谢实验部/2023"
            local_path = str(cwd / "rawdata")
            filename = project_path.replace("项目报告", "原始数据")
            args = ["tool_UpDataToObs", "-s", server_path, "-m", "-p", local_path, "-f", filename]
            r = sp.run(args, cwd=str(cwd))
            if r.returncode != 0:
                raise ShellRunTimeError("update rawdata to HuaWei OBS server failed！")
            logger.info("upload rawdata to HuaWei OBS server successfully!")
        else:
            logger.info("Local project, skip upload rawdata to HuaWei OBS server!")

    def preprocess_qi(self, using_herb=False, treshold_tcm=1000, interested_cpd=None, peak_id=None):
        MIDWARE_PATH = Path(self.data) / "midware"
        MIDWARE_PATH.mkdir(exist_ok=True)
        """ QI result process"""
        logger.info("正在预处理QI结果，添加代谢物注释...")
        mode_used, sample_info = self.parse_meta()

        tcm_cols = self.parse_exp_col()
        n_tcm_cols = len(tcm_cols)

        data_processed = []
        interested_records = []
        for mode in mode_used:
            src = Path(self.data)
            id_path = src / "QIdata" / mode / sample_info.loc[mode, "QI定性结果"]
            m_path = src / "QIdata" / mode / sample_info.loc[mode, "QI定量结果"]
            # read data
            qualitative = pd.read_csv(id_path, encoding='utf-8')
            quantitative = pd.read_csv(m_path, encoding='utf-8', header=[2, ])
            # drop columns
            qualitative = qualitative.drop(
                columns=["Accepted?", "Link", "Neutral mass (Da)", "m/z", "Charge", "Retention time (min)",
                         "Isotope Similarity", "Theoretical Isotope Distribution"])
            quantitative = quantitative.drop(
                # columns[13:16] rm Raw Abundance
                # TODO: 13:13+n
                columns=[*quantitative.columns[-10:], *quantitative.columns[13:13+n_tcm_cols], "Neutral mass (Da)",
                         "Chromatographic peak width (min)",
                         "Identifications", "Isotope Distribution", "Maximum Abundance", "Minimum CV%"])
            # merge and del rows with na
            qual_quant = pd.merge(qualitative, quantitative, on="Compound", how="inner").dropna(axis=0, how="any")
            # filter with different rules when dealing with various db
            qual_quant["DB"] = qual_quant["Description"].apply(lambda x: x.split(":", 1)[0])
            qual_quant["Description"] = qual_quant["Description"].apply(lambda x: x.split(":", 1)[1])
            # add annotation level
            qual_quant["Level"] = qual_quant.apply(
                lambda x: "level3" if x["DB"] == "HERB" else "level1" if x["Fragmentation Score"] >= 50 else "level2",
                axis=1)
            # add theory mz
            qual_quant["theoretical m/z"] = (1e6 * qual_quant["m/z"]) / (1e6 + qual_quant["Mass Error (ppm)"])
            # TODO: retain interested compounds
            if interested_cpd:
                _tmp = qual_quant.query("Description in @interested_cpd")
                if peak_id:
                    _nmap = dict(zip(interested_cpd, peak_id))
                    for i, row in _tmp.iterrows():
                        if _nmap.get(row["Description"]) != row["Compound"]:
                            _tmp.drop(index=i, inplace=True)
                if _tmp.shape[0]:
                    _tmp["Ion mode"] = mode
                interested_records.append(_tmp)
            # TCM
            qual_quant_tcm = qual_quant.query("DB == 'TCM'")
            # qual_quant_tcm = qual_quant_tcm.query("Score >= 50 | (Score >= 40 & `Fragmentation Score` >= 50)")
            qual_quant_tcm = qual_quant_tcm.query("Score >= 40")
            # Animal
            qual_quant_animal = qual_quant.query("DB == 'ANIMAL'")
            qual_quant_animal = qual_quant_animal.\
                query("(Score >= 50 & `Fragmentation Score` >= 50) | (Score >= 40 & `Fragmentation Score` >= 50)")
            # filter with rules and Sort
            if using_herb:
                # Herb
                qual_quant_herb = qual_quant.query("DB == 'HERB'")
                qual_quant_herb = qual_quant_herb.query("Score >= 50 | (Score >= 40 & `Fragmentation Score` >= 50)")

                qual_quant = pd.concat([qual_quant_tcm, qual_quant_animal, qual_quant_herb]).sort_values(
                "Compound").reset_index(drop=True)
            else:
                qual_quant = pd.concat([qual_quant_tcm, qual_quant_animal]).sort_values(
                "Compound").reset_index(drop=True)
            # filter with rules and Sort
            # qual_quant = qual_quant.query("Score >= 50 | (Score >= 40 & `Fragmentation Score` >= 50)").sort_values(
            #     "Compound").reset_index(drop=True)
            # select best cpd
            best_index = self.selectIndex(qual_quant, "Compound")
            qual_quant = qual_quant.loc[best_index, :]
            qual_quant["Ion mode"] = mode
            data_processed.append(qual_quant)

        res = self.fill_qi_df(data_processed)
        # abundance filter
        res = res[res[tcm_cols].mean(axis=1) > treshold_tcm].reset_index(drop=True)
        if interested_records:
            if not pd.concat(interested_records).empty:
                interested_records = self.fill_qi_df(interested_records, dedup=False)
                tcm_percent_cols = []
                for col in tcm_cols:
                    interested_records[f"{col}（峰面积的比值%）"] = 100 * interested_records[col] / interested_records[col].sum()
                    tcm_percent_cols.append(f"{col}（峰面积的比值%）")
                # add M mean
                interested_records["M均值（峰面积的比值%）"] = interested_records[tcm_percent_cols].mean(axis=1)
                tcm_percent_cols.append("M均值（峰面积的比值%）")
                # add compound No.
                interested_records["No."] = pd.Series(data=[f"cpdInterested{x:0>3d}" for x in range(1, interested_records.shape[0] + 1)])
                # arrange columns
                interested_records = interested_records[["No.", "Compound", "Compound ID", "Adducts", "Formula", "Score", "Fragmentation Score",
                           "Mass Error (ppm)",
                           "Description", "theoretical m/z", "m/z", "Retention time (min)", "Ion mode", *tcm_cols,
                           *tcm_percent_cols, "InChIKey", "SMILES", "HMDB", "METLIN",
                           "LipidMaps", "KEGG", "PubChem", "CAS", "中文名", "中文大类", "中文子类", "英文大类", "英文子类", "厂家", "货号", "纯度",
                           "DB", "Level"]]

                interested_records["lower_name"] = interested_records['Description'].apply(lambda x: x.lower())
                interested_records = interested_records.sort_values(by="lower_name").reset_index(drop=True)
                indexOfBestName = self.selectIndexOfBestName(interested_records, 'lower_name')
                interested_records = interested_records.loc[indexOfBestName, :].reset_index(drop=True).applymap(
                    lambda x: x.strip() if isinstance(x, str) else x)
                interested_records.sort_values(by="No.").to_excel(str(self.midware / "interested_cpd.xlsx"), index=False)
        # add relative M <占位>
        tcm_percent_cols = []
        for col in tcm_cols:
            res[f"{col}（峰面积的比值%）"] = 100 * res[col] / res[col].sum()
            tcm_percent_cols.append(f"{col}（峰面积的比值%）")
        # add M mean
        res["M均值（峰面积的比值%）"] = res[tcm_percent_cols].mean(axis=1)
        tcm_percent_cols.append("M均值（峰面积的比值%）")

        # 根据名称去重
        res["lower_name"] = res['Description'].apply(lambda x: x.lower())
        res = res.sort_values(by="lower_name").reset_index(drop=True)
        indexOfBestName = self.selectIndexOfBestName(res, 'lower_name')
        res = res.loc[indexOfBestName, :].reset_index(drop=True).applymap(lambda x:x.strip() if isinstance(x, str) else x)
        # add compound No.
        res["No."] = pd.Series(data=[f"compound{x:0>4d}" for x in range(1, res.shape[0] + 1)])
        # arrange columns
        res = res[["No.", "Compound", "Compound ID", "Adducts", "Formula", "Score", "Fragmentation Score", "Mass Error (ppm)",
                   "Description", "theoretical m/z", "m/z", "Retention time (min)", "Ion mode", *tcm_cols, *tcm_percent_cols, "InChIKey", "SMILES", "HMDB", "METLIN",
                   "LipidMaps", "KEGG", "PubChem", "CAS", "中文名", "中文大类", "中文子类", "英文大类", "英文子类", "厂家", "货号", "纯度", "DB", "Level"]]
        # export
        dirname, project_info = self.parse_project()

        is_del_animal = False if project_info.loc[0, "复方中是否包含动物"] == '是' else True
        if is_del_animal:
            with self.ENGINE.connect() as tx:
                df = pd.read_sql('herb_annotation', tx)
                dels = df[(df.is_animal_source.astype(int) + df.is_lipid.astype(int)) >= 1]["inchikey"].to_list()
            res = res.query('`Compound ID` not in @dels').reset_index(drop=True)
        res = res.sort_values(by="Level").reset_index(drop=True)
        SAVE_PATH = Path(self.data) / dirname
        # if not SAVE_PATH.exists():
        SAVE_PATH.mkdir(exist_ok=True)
        (SAVE_PATH / "实验内容").mkdir(exist_ok=True)
        # (SAVE_PATH / "鹿明信息").mkdir(exist_ok=True)
        (SAVE_PATH / "数据分析").mkdir(exist_ok=True)
        (SAVE_PATH / "数据分析" / "1.紫外吸收图").mkdir(exist_ok=True)
        (SAVE_PATH / "数据分析" / "2.色谱图").mkdir(exist_ok=True)
        (SAVE_PATH / "数据分析" / "3.质谱图").mkdir(exist_ok=True)
        # (SAVE_PATH / "数据分析" / "3.质谱图" / "png").mkdir(exist_ok=True)
        # (SAVE_PATH / "数据分析" / "3.质谱图" / "pdf").mkdir(exist_ok=True)
        (SAVE_PATH / "数据分析" / "3.质谱图" / "中药原方成分" / "LuMet-TCM" / "png").mkdir(exist_ok=True, parents=True)
        (SAVE_PATH / "数据分析" / "3.质谱图" / "中药原方成分" / "LuMet-TCM" / "pdf").mkdir(exist_ok=True, parents=True)
        if using_herb:
            (SAVE_PATH / "数据分析" / "3.质谱图" / "中药原方成分" / "HerbDB" / "png").mkdir(exist_ok=True, parents=True)
            (SAVE_PATH / "数据分析" / "3.质谱图" / "中药原方成分" / "HerbDB" / "pdf").mkdir(exist_ok=True, parents=True)
        (SAVE_PATH / "数据分析" / "4.定性定量结果").mkdir(exist_ok=True)
        (SAVE_PATH / "数据分析" / "5.中药成分分类").mkdir(exist_ok=True)
        res.to_excel(str(SAVE_PATH / "数据分析" / "4.定性定量结果" / "数据矩阵.xlsx"), sheet_name="data", encoding="utf-8",
                     index=False)
        res.to_excel(str(MIDWARE_PATH / "preprocessed.xlsx"), sheet_name="data", encoding="utf-8",
                     index=False)

    def loadMS2FromMsp(self, path: Path, mode: list) -> dict:
        """
        load MS2 from MSP files
        :param path: pathlib.Path("../../MSP")
        :param mode: ["POS",["NEG"]
        :return:
        """
        MS2 = dict()
        for m in mode:
            # MS2spec container
            spectra = dict()
            # msp path
            msp = str(list(path.rglob(f"*{m}*"))[0])
            for spec in self.msp_parser(msp):
                _id = spec.metadata["Comment"]
                if spectra.get(_id):
                    if isinstance(spectra.get(_id), MS2spec):
                        spectra[_id] = spectra.get(_id) + spec
                    elif isinstance(spectra.get(_id), list):
                        spectra[_id] = spectra.get(_id) + [spec]
                else:
                    spectra[_id] = spec
            MS2[m] = spectra
        return MS2

    def readMS2FromDB(self, inchikey: str, ion_mode: str, adducts: list,
                      exp_spec: Union[None, MSSpectrum] = None, is_product=False) -> Union[MSSpectrum, None, bool]:
        conn = self.ENGINE.connect()
        for adduct in adducts:
            sql = f"""select * from herb_product where `NAME`={repr(inchikey)} and `ion_mode`={repr(ion_mode)} and adduct={repr(adduct)}""" if is_product \
                else f"""select * from herb_msms where inchikey={repr(inchikey)} and `mode`={repr(ion_mode)} and adduct={repr(adduct)}"""
            res = pd.read_sql(text(sql), conn)
            if res.shape[0] == 0:
                continue
            mzs = [float(x) for x in res.loc[0, "mz"].split(",")]
            intensities = [float(x) for x in res.loc[0, "intensity"].split(",")]
            spectrum = MSSpectrum()
            spectrum.setMSLevel(2)
            spectrum.set_peaks((mzs, intensities))
            conn.close()
            return spectrum
        # return None
        else:
            sql = f"""select * from herb_product where `NAME`={repr(inchikey)} and `ion_mode`={repr(ion_mode)}""" if is_product \
                else f"""select * from herb_msms where inchikey={repr(inchikey)} and `mode`={repr(ion_mode)}"""
            res = pd.read_sql(text(sql), conn)
            conn.close()
            if res.shape[0] == 0:
                return None
            else:
                if exp_spec is None:
                    return True
                else:
                    scores, spectra = [], []
                    for i, row in res.iterrows():
                        mzs = [float(x) for x in row["mz"].split(",")]
                        intensities = [float(x) for x in row["intensity"].split(",")]
                        spectrum = MSSpectrum()
                        spectrum.setMSLevel(2)
                        spectrum.set_peaks((mzs, intensities))
                        ref_spec = self.spec_normalize(spectrum)
                        # 计算相似性得分
                        scorer = SpectrumAlignmentScore()
                        param = scorer.getParameters()
                        param.setValue("tolerance", 5.0)
                        param.setValue("is_relative_tolerance", "true")
                        scorer.setParameters(param)
                        # exp_spec must be normalized
                        score = scorer(ref_spec, exp_spec)
                        scores.append(score)
                        spectra.append(ref_spec)
                    best = np.asarray(scores).argmax()
                    return spectra[best]

    def plot_tic(self, dpi: int = 200, fmt: str = "png",
                 absolute_intensity: bool = True):
        """
        plot bpc
        :param dpi: dpi
        :param fmt: picture format
        :param absolute_intensity: y axis
        :return: None
        """
        logger.info("正在绘制TIC图...")
        _type = "中药复方" if self._project == "TCM" else "中药入血成分"
        mode_used, sample_info = self.parse_meta()
        project, project_info = self.parse_project()

        for mode in mode_used:
            # prepare files
            if self._project == 'TCM':
                dir_path = Path(self.data) / "mzmldata" / mode
                files = [str(path) for path in dir_path.rglob("*.mzML")]
            else:
                control = Path(self.data) / "mzmldata" / mode / sample_info["Control"].loc[mode, "质谱数据1"]
                blood = Path(self.data) / "mzmldata" / mode / sample_info["Blood"].loc[mode, "质谱数据1"]
                tcm = Path(self.data) / "mzmldata" / mode / sample_info["TCM"].loc[mode, "质谱数据1"]
                files = [control, blood, tcm]
                files = [str(f) for f in files]
            # plot
            n = len(files)
            fig, ax = plt.subplots(n, 1, figsize=(8, 2 * n), dpi=dpi, gridspec_kw={'hspace': 0.4})
            ax = [ax] if n == 1 else ax
            plt.style.use(style="default")
            colormap = plt.get_cmap('Set1')(range(n))

            for i, file in enumerate(files):
                exp: MSExperiment = MSExperiment()
                MzMLFile().load(file, exp)
                tic = exp.calculateTIC()
                retention_times, intensities = tic.get_peaks()
                intensities = intensities if absolute_intensity else intensities * (100 / intensities.max())
                retention_times = retention_times / 60

                axis = ax[i]
                # 设置 轴刻度
                if absolute_intensity:
                    axis.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
                    axis.set_ylabel('Abundance')
                    axis.set_xlim([0, 16])
                else:
                    axis.axis([0, 16, 0, 100])
                    axis.yaxis.set_major_locator(MultipleLocator(20))
                    axis.yaxis.set_minor_locator(MultipleLocator(5))
                    axis.set_ylabel('Relative Abundance')

                axis.xaxis.set_major_locator(MultipleLocator(4))
                axis.xaxis.set_minor_locator(MultipleLocator(1))

                axis.spines['right'].set_visible(False)
                axis.spines['top'].set_visible(False)


                if self._project == 'TCM':
                    axis.plot(retention_times, intensities, color=colormap[i], lw=0.8)
                else:
                    axis.plot(retention_times, intensities, color=colormap[i], lw=0.8, label=self.group[i])
                    axis.legend(fontsize=6)
                if i == n-1:
                    axis.set_xlabel('Time (min)')
                if absolute_intensity:
                    yticks = axis.get_yticks()
                    axis.set_ylim([0, yticks[-1]])

            ax[0].set_title(f"TIC ({mode})")

            save_path = str(Path(self.data) / project / "数据分析" / "2.色谱图" / f"{_type}TIC图-{mode}")
            plt.savefig(f"{save_path}.{fmt}")
            plt.savefig(f"{save_path}.pdf")
            plt.close()

    def plot_tic_split(self, dpi: int = 200, fmt: str = "png",
                       absolute_intensity: bool = True, add_peak_idx: bool = False):
        """
        plot bpc
        :param dpi: dpi
        :param fmt: picture format
        :return: None
        """
        logger.info("正在绘制单张TIC图...")
        _type = "中药复方" if self._project == "TCM" else "中药入血成分"
        mode_used, sample_info = self.parse_meta()
        project, project_info = self.parse_project()

        for mode in mode_used:
            # prepare files
            if self._project == 'TCM':
                dir_path = Path(self.data) / "mzmldata" / mode
                files = [str(path) for path in dir_path.rglob("*.mzML")]
            else:
                control = Path(self.data) / "mzmldata" / mode / sample_info["Control"].loc[mode, "质谱数据1"]
                blood = Path(self.data) / "mzmldata" / mode / sample_info["Blood"].loc[mode, "质谱数据1"]
                tcm = Path(self.data) / "mzmldata" / mode / sample_info["TCM"].loc[mode, "质谱数据1"]
                files = [control, blood, tcm]
                files = [str(f) for f in files]
            # plot
            n = len(files)
            plt.style.use(style="default")
            colormap = plt.get_cmap('Set1')(range(n))

            for i, file in enumerate(files):
                exp: MSExperiment = MSExperiment()
                MzMLFile().load(file, exp)
                tic = exp.calculateTIC()

                retention_times, intensities = tic.get_peaks()
                retention_times = retention_times / 60
                intensities = intensities if absolute_intensity else intensities * (100 / intensities.max())

                fig, axis = plt.subplots(1, 1, figsize=(12, 8), dpi=dpi)

                # 设置 轴刻度
                if absolute_intensity:
                    axis.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
                    axis.set_ylabel('Abundance')
                    axis.set_xlim([0, 16])
                else:
                    axis.axis([0, 16, 0, 100])
                    axis.yaxis.set_major_locator(MultipleLocator(20))
                    axis.yaxis.set_minor_locator(MultipleLocator(5))
                    axis.set_ylabel('Relative Abundance')

                axis.xaxis.set_major_locator(MultipleLocator(4))
                axis.xaxis.set_minor_locator(MultipleLocator(1))

                axis.spines['right'].set_visible(False)
                axis.spines['top'].set_visible(False)
                if self._project == 'TCM':
                    axis.plot(retention_times, intensities, color=colormap[i], lw=0.8)
                else:
                    axis.plot(retention_times, intensities, color=colormap[i], lw=0.8, label=self.group[i])
                    axis.legend(fontsize=6)
                axis.set_xlabel('Time (min)')

                if absolute_intensity:
                    yticks = axis.get_yticks()
                    axis.set_ylim([0, yticks[-1]])

                    # y_min, y_max = axis.get_ylim()
                    # y_max = np.format_float_scientific(y_max, precision=2, exp_digits=1)
                    # y_max = self.sci_big(y_max)
                    # axis.set_ylim([0, y_max])

                if add_peak_idx:
                    x, y = retention_times, intensities
                    _, info = find_peaks(y, width=(3,), prominence=0.02 * 1e10)
                    idx = 1
                    for a, b in zip(x[_], y[_]):
                        axis.text(a, b, str(idx), fontsize=8, ha="center", va="bottom")
                        idx = idx + 1
                    sp = str(Path(self.data) / project / "数据分析" / "2.色谱图" / f"{_type}TIC图-Peak_index_rt_single-{mode}-{i + 1}.xlsx")
                    pd.DataFrame(data={"No.Peak": range(1, x[_].size + 1), "rt (min)": x[_]}).to_excel(sp, index=False)
                axis.set_title(f"TIC ({mode})")

                save_path = str(Path(self.data) / project / "数据分析" / "2.色谱图" / f"{_type}TIC图-single-{mode}-{i+1}")
                plt.savefig(f"{save_path}.{fmt}")
                plt.savefig(f"{save_path}.pdf")
                plt.close()

    def plot_bpc(self, dpi: int = 200, fmt: str = "png", absolute_intensity: bool = True, water_mask: bool = False):
        """
        plot bpc
        :param dpi: dpi
        :param fmt: picture format
        :return: None
        """
        logger.info("正在绘制BPC图...")
        _type = "中药复方" if self._project == "TCM" else "中药入血成分"
        mode_used, sample_info = self.parse_meta()
        project, project_info = self.parse_project()

        for mode in mode_used:
            # prepare files
            if self._project == 'TCM':
                dir_path = Path(self.data) / "mzmldata" / mode
                files = [str(path) for path in dir_path.rglob("*.mzML")]
            else:
                control = Path(self.data) / "mzmldata" / mode / sample_info["Control"].loc[mode, "质谱数据1"]
                blood = Path(self.data) / "mzmldata" / mode / sample_info["Blood"].loc[mode, "质谱数据1"]
                tcm = Path(self.data) / "mzmldata" / mode / sample_info["TCM"].loc[mode, "质谱数据1"]
                files = [control, blood, tcm]
                files = [str(f) for f in files]
            # plot
            n = len(files)
            fig, ax = plt.subplots(n, 1, figsize=(8, 2 * n), dpi=dpi, gridspec_kw={'hspace': 0.4})
            ax = [ax] if n == 1 else ax
            plt.style.use(style="default")
            colormap = plt.get_cmap('Set1')(range(n))

            for i, file in enumerate(files):
                exp: MSExperiment = MSExperiment()
                MzMLFile().load(file, exp)
                retention_times = []
                intensities = []
                for spec in exp:
                    if spec.getMSLevel() == 1:
                        retention_times.append(spec.getRT())
                        intensities.append(max(spec.get_peaks()[1]))
                intensities = np.asarray(intensities)
                intensities = intensities if absolute_intensity else intensities * (100 / intensities.max())
                retention_times = np.asarray(retention_times) / 60

                axis = ax[i]
                # 设置 轴刻度
                if absolute_intensity:
                    axis.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
                    axis.set_ylabel('Abundance')
                    axis.set_xlim([0, 16])
                else:
                    axis.axis([0, 16, 0, 100])
                    axis.yaxis.set_major_locator(MultipleLocator(20))
                    axis.yaxis.set_minor_locator(MultipleLocator(5))
                    axis.set_ylabel('Relative Abundance')

                axis.xaxis.set_major_locator(MultipleLocator(4))
                axis.xaxis.set_minor_locator(MultipleLocator(1))

                axis.spines['right'].set_visible(False)
                axis.spines['top'].set_visible(False)

                if self._project == 'TCM':
                    axis.plot(retention_times, intensities, color=colormap[i], lw=0.8)
                else:
                    axis.plot(retention_times, intensities, color=colormap[i], lw=0.8, label=self.group[i])
                    axis.legend(fontsize=6)
                if i == n-1:
                    axis.set_xlabel('Time (min)')
                if absolute_intensity:
                    yticks = axis.get_yticks()
                    axis.set_ylim([0, yticks[-1]])

            ax[0].set_title(f"BPC ({mode})")

            if water_mask:
                fig.text(0.5, 0.5, 'created by lumingbio',
                         fontsize=40, color='gray', alpha=0.5,
                         ha='center', va='center', rotation=30)

            save_path = str(Path(self.data) / project / "数据分析" / "2.色谱图" / f"{_type}BPC图-{mode}")
            plt.savefig(f"{save_path}.{fmt}")
            plt.savefig(f"{save_path}.pdf")
            plt.close()

    def plot_bpc_split(self, dpi: int = 200, fmt: str = "png", absolute_intensity=True, water_mask: bool = False):
        """
        plot bpc
        :param dpi: dpi
        :param fmt: picture format
        :return: None
        """
        logger.info("正在绘制单张BPC图...")
        _type = "中药复方" if self._project == "TCM" else "中药入血成分"
        mode_used, sample_info = self.parse_meta()
        project, project_info = self.parse_project()

        for mode in mode_used:
            # prepare files
            if self._project == 'TCM':
                dir_path = Path(self.data) / "mzmldata" / mode
                files = [str(path) for path in dir_path.rglob("*.mzML")]
            else:
                control = Path(self.data) / "mzmldata" / mode / sample_info["Control"].loc[mode, "质谱数据1"]
                blood = Path(self.data) / "mzmldata" / mode / sample_info["Blood"].loc[mode, "质谱数据1"]
                tcm = Path(self.data) / "mzmldata" / mode / sample_info["TCM"].loc[mode, "质谱数据1"]
                files = [control, blood, tcm]
                files = [str(f) for f in files]
            # plot
            n = len(files)
            plt.style.use(style="default")
            colormap = plt.get_cmap('Set1')(range(n))

            for i, file in enumerate(files):
                exp: MSExperiment = MSExperiment()
                MzMLFile().load(file, exp)
                retention_times = []
                intensities = []
                for spec in exp:
                    if spec.getMSLevel() == 1:
                        retention_times.append(spec.getRT())
                        intensities.append(max(spec.get_peaks()[1]))
                intensities = np.asarray(intensities)
                intensities = intensities if absolute_intensity else intensities * (100 / intensities.max())
                retention_times = np.asarray(retention_times) / 60
                fig, axis = plt.subplots(1, 1, figsize=(12, 8), dpi=dpi)

                # 设置 轴刻度
                if absolute_intensity:
                    axis.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
                    axis.set_ylabel('Abundance')
                    axis.set_xlim([0, 16])
                else:
                    axis.axis([0, 16, 0, 100])
                    axis.yaxis.set_major_locator(MultipleLocator(20))
                    axis.yaxis.set_minor_locator(MultipleLocator(5))
                    axis.set_ylabel('Relative Abundance')

                axis.xaxis.set_major_locator(MultipleLocator(4))
                axis.xaxis.set_minor_locator(MultipleLocator(1))

                axis.spines['right'].set_visible(False)
                axis.spines['top'].set_visible(False)
                if self._project == 'TCM':
                    axis.plot(retention_times, intensities, color=colormap[i], lw=0.8)
                else:
                    axis.plot(retention_times, intensities, color=colormap[i], lw=0.8, label=self.group[i])
                    axis.legend(fontsize=6)
                axis.set_xlabel('Time (min)')

                axis.set_title(f"BPC ({mode})")
                if absolute_intensity:
                    yticks = axis.get_yticks()
                    axis.set_ylim([0, yticks[-1]])
                if water_mask:
                    fig.text(0.5, 0.5, 'created by lumingbio',
                             fontsize=40, color='gray', alpha=0.5,
                             ha='center', va='center', rotation=30)
                save_path = str(Path(self.data) / project / "数据分析" / "2.色谱图" / f"{_type}BPC图-single-{mode}-{i+1}")
                plt.savefig(f"{save_path}.{fmt}")
                plt.savefig(f"{save_path}.pdf")
                plt.close()


    def plot_chromatograms(self, dpi: int = 200, fmt: str = "png", water_mask: bool = False, absolute_intensity=True):
        self.plot_bpc(dpi=dpi, fmt=fmt, water_mask=water_mask, absolute_intensity=absolute_intensity)
        self.plot_bpc_split(dpi=dpi, fmt=fmt, water_mask=water_mask, absolute_intensity=absolute_intensity)
        self.plot_tic(dpi=dpi, fmt=fmt, absolute_intensity=absolute_intensity)
        self.plot_tic_split(dpi=dpi, fmt=fmt, absolute_intensity=absolute_intensity)

    def plot_xic(self, ppm: int = 5, fmt: str = "png",  water_mask: bool = False, interested_cpd=None, **kwargs):
        dpi = 600
        img_size = 250
        if kwargs.get("dpi"):
            dpi = int(kwargs.get("dpi"))
        if kwargs.get("img_size"):
            img_size = int(kwargs.get("img_size"))

        # 获取项目信息
        project, project_info = self.parse_project()
        # 1.MSExp container for ex xic
        MSExps = dict()
        # mzML path
        mode_used, sample_info = self.parse_meta()
        for mode in mode_used:
            # TODO: 确保 mzML PATH 中没有中文, 否则会 error.
            mzml_path = Path(self.data) / "mzmldata" / mode
            # use the first to extract EIC chr
            mzml = [str(path) for path in mzml_path.rglob("*.mzML")][0]
            # load mzml
            exp = MSExperiment()
            MzMLFile().load(mzml, exp)
            MSExps[mode] = exp

        # 2.MS2 container for ref MS/MS
        msp_path = Path(self.data) / "QIdata" / "MSP"
        refMS2: dict = self.loadMS2FromMsp(msp_path, mode_used)

        # 3.iter DataFrame of metabolites
        project_name, project_info = self.parse_project()
        data_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "数据矩阵.xlsx") if not interested_cpd else str(self.midware / "interested_cpd.xlsx")
        # load peak table
        peaks = pd.read_excel(data_path, sheet_name=0)
        # plot each line
        fragment_ions = []
        # TODO: 速度太慢改用多线程
        # add progress bar
        with tqdm(zip(peaks["No."], peaks["Compound"], peaks["Compound ID"],peaks["Adducts"],peaks["Formula"], peaks["Description"],
                      peaks["m/z"],peaks["Retention time (min)"],peaks["Ion mode"],peaks["DB"]), total=peaks.shape[0],
                  desc="INFO - Plot MS2") as pb:
            for cpdn, _id, inchikey, adduct, formula, name, _mz, _rt, ion_mode, db in pb:
                pb.set_postfix(compound=name)
                # logger.info("正在绘制代谢物%s的MS/MS图...", name)
                fragment_ions_str = None

                # 3.1 extract data from mzML & plot XIC
                exp = MSExps[ion_mode]
                mass_tolerance = _mz * ppm / 1e6
                _rt, _mz = float(_rt), float(_mz)
                # EIC
                retention_times = []
                intensities = []
                for spec in exp:
                    if spec.getMSLevel() == 1:
                        filtered_int = []
                        for mz, i in zip(*spec.get_peaks()):
                            if _mz - mass_tolerance < mz < _mz + mass_tolerance:
                                filtered_int.append(i)
                        rt = spec.getRT()
                        ints = sum(filtered_int)
                        retention_times.append(rt)
                        intensities.append(ints)

                intensities = np.asarray(intensities)
                retention_times = np.asarray(retention_times) / 60
                # pyplot top plot
                fig, axes = plt.subplots(2, 1, figsize=(8, 7.4), dpi=dpi, gridspec_kw={'height_ratios': [1, 2.7]})
                # 设置白底
                fig.patch.set_alpha(1.)
                if water_mask:
                    fig.text(0.5, 0.5, 'created by lumingbio',
                             fontsize=40, color='gray', alpha=0.5,
                             ha='center', va='center', rotation=30)
                ax = axes[0]
                # ratio_xy = 16 / 100
                # fig_aspect = 16 / (100*682/755)
                # aspect = ratio_xy * fig_aspect
                # ax.set_aspect(aspect)

                ax1 = axes[1]
                # 设置 轴刻度
                # ax.yaxis.set_major_locator(MultipleLocator(20))
                # ax.yaxis.set_minor_locator(MultipleLocator(5))
                # ax.axis([0, 16, 0, 100])
                ax.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
                ax.set_xlim([0, 16])

                ax.xaxis.set_major_locator(MultipleLocator(4))
                ax.xaxis.set_minor_locator(MultipleLocator(1))
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.plot(retention_times, intensities, lw=0.8, color=plt.get_cmap('Set1')(range(3))[1])
                # ax.plot(retention_times, intensities * (100 / intensities.max()), lw=0.8, color=plt.get_cmap('Set1')(range(3))[1])
                yticks = ax.get_yticks()
                ax.set_ylim([0, yticks[-1]])
                ax.set_ylabel('Abundance')
                # ax.set_ylabel('Relative Abundance')
                title = self.filename_handler(name, escape=False)
                ax.set_title(f"{title}: Mode_{ion_mode}  RT_{_rt:.2f} min  MZ_{_mz:.4f}")

                # new: check ms1 20230623
                rts, its = retention_times, intensities * (100 / intensities.max())
                is_ms1_valid = self.check_ms1(rts, its, _rt) if not interested_cpd else True
                if not is_ms1_valid:
                    # pic_name = self.filename_handler(name, escape=True, full=True)
                    # save_path_png = Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "delete"
                    # if not save_path_png.exists():
                    #     save_path_png.mkdir(parents=True, exist_ok=True)
                    # ax.axvline(_rt, color='green')
                    # plt.savefig(f"{str(save_path_png)}/{pic_name}.{fmt}")
                    fragment_ions_str = 'delete'
                else:
                    # 3.2 query data & plot mirror or MS2
                    adducts = adduct.split(", ")
                    expMS2 = refMS2[ion_mode].get(_id, "")

                    precursor = self.get_parent_smiles(inchikey)

                    if self.is_cpd_valid(inchikey):
                        if len(expMS2) > 0:
                            exp_MS2 = expMS2[0] if isinstance(expMS2, list) else expMS2
                            exp_spec = MSSpectrum()
                            exp_spec.set_peaks((exp_MS2.mz, exp_MS2.intensity))
                            exp_spec = self.spec_normalize(exp_spec)

                            ref_spec = self.readMS2FromDB(inchikey, ion_mode.lower(), adducts, exp_spec=exp_spec)
                            if ref_spec:
                                # TODO: 在HERB MS2库出来之前临时添加，后续需要删除
                                if db == "HERB":
                                    # herb找到std二级，说明RT不对，需要删除
                                    fragment_ions_str = "delete"
                                else:
                                    ref_spec = self.spec_normalize(ref_spec)
                                    _smiles: Union[None, list] = self.get_smiles(ref_spec, self.ENGINE, formula, inchikey,
                                                                                 name, ion_mode, adducts)
                                    if _smiles:
                                        fragment_ions_str = self.plot_mirror_spectrum_with_smiles(ax1, exp_spec, ref_spec, _smiles, precursor, _mz,
                                                                                  img_size=img_size)
                                    else:
                                        fragment_ions_str = self.plot_mirror_spectrum_with_smiles(ax1, exp_spec, ref_spec, [], precursor, _mz,
                                                                              img_size=img_size)
                                # if not fragment_ions_str:
                                #     mzs, intensities = ref_spec.get_peaks()
                                #     _mzs, _intensities = mzs[np.argsort(intensities)], intensities[np.argsort(intensities)]
                                #     _mzs = _mzs if len(_mzs) <= 5 else _mzs[:6]
                                #     _mzs = _mzs[np.argsort(_mzs)]
                                #     fragment_ions_str = "; ".join([f"{m:.4f}" for m in _mzs])
                            else:
                                # std无MS2，实测有 MS2
                                if db in ["TCM", "ANIMAL"]:
                                    fragment_ions_str = "delete"
                                    if interested_cpd:
                                        # 非标品库绘制
                                        _smiles: Union[None, list] = self.get_smiles(exp_spec, self.ENGINE, formula,
                                                                                     inchikey, name,
                                                                                     ion_mode, adducts)
                                        if _smiles:
                                            fragment_ions_str = self.plot_spectrum_with_smiles(ax1, exp_spec, _smiles,
                                                                                               precursor,
                                                                                               _mz,
                                                                                               img_size=img_size)
                                        else:
                                            fragment_ions_str = self.plot_spectrum(ax1, exp_spec, precursor, _mz,
                                                                                   img_size=img_size)
                                else:
                                    # 非标品库绘制
                                    _smiles: Union[None, list] = self.get_smiles(exp_spec, self.ENGINE, formula, inchikey, name,
                                                                                 ion_mode, adducts)
                                    if _smiles:
                                        fragment_ions_str = self.plot_spectrum_with_smiles(ax1, exp_spec, _smiles, precursor,
                                                                                           _mz,
                                                                                           img_size=img_size)
                                    else:
                                        fragment_ions_str = self.plot_spectrum(ax1, exp_spec, precursor, _mz, img_size=img_size)
                                # ADD_TO_MS2_DIR
                                # best_adduct = Report.find_best_adduct(_mz, precursor, ion_mode)
                                # if best_adduct:
                                #     tmp_dir = Path(self.MS2_DIR) / project
                                #     if not tmp_dir.exists():
                                #         tmp_dir.mkdir(exist_ok=True)
                                #     fn = inchikey + "_" + ion_mode + "_" + best_adduct
                                #     file = tmp_dir / fn
                                #     mzs = ",".join([str(m) for m in exp_spec.get_peaks()[0]])
                                #     ints = ",".join([str(i) for i in exp_spec.get_peaks()[1]])
                                #     with open(str(file), "w", encoding="utf8") as f:
                                #         line = f"{inchikey}\t{formula}\t{name}\t{_mz}\t{_rt}\t{ion_mode}\t{precursor}\t{mzs}\t{ints}\n"
                                #         f.write(line)
                                #         for m, i in zip(*exp_spec.get_peaks()):
                                #             f.write(f"{m}\t{i}\n")

                                # _smiles: Union[None, list] = self.get_smiles(exp_spec, self.ENGINE, formula,
                                #                                              inchikey, name,
                                #                                              ion_mode, adducts)
                                # if _smiles:
                                #     # fragment_ions_str =
                                #     self.plot_spectrum_with_smiles(ax1, exp_spec, _smiles, precursor, _mz, img_size=img_size)
                                # else:
                                #     self.plot_spectrum(ax1, exp_spec, precursor, _mz, img_size=img_size)
                                # if not fragment_ions_str:
                                #     mzs, intensities = exp_spec.get_peaks()
                                #     _mzs, _intensities = mzs[np.argsort(intensities)], intensities[np.argsort(intensities)]
                                #     _mzs = _mzs if len(_mzs) <= 5 else _mzs[:6]
                                #     _mzs = _mzs[np.argsort(_mzs)]
                                #     fragment_ions_str = "; ".join([f"{m:.4f}" for m in _mzs])
                        else:
                            self.plot_defect(ax1, precursor, _mz, img_size=img_size)
                            fragment_ions_str = ""
                            #     ref_spec = self.spec_normalize(ref_spec)
                            #     mzs, intensities = ref_spec.get_peaks()
                            #     _mzs, _intensities = mzs[np.argsort(intensities)], intensities[np.argsort(intensities)]
                            #     _mzs = _mzs if len(_mzs) <= 5 else _mzs[:6]
                            #     _mzs = _mzs[np.argsort(_mzs)]
                            #     fragment_ions_str = "; ".join([f"{m:.4f}" for m in _mzs])
                    else:
                        fragment_ions_str = 'delete'

                fragment_ions.append(fragment_ions_str if fragment_ions_str else '')
                # save fig
                if fragment_ions_str != 'delete':
                    pic_name = cpdn
                    save_path_png = str(
                        Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "LuMet-TCM" / "png") if db in ["TCM","ANIMAL"] \
                        else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "HerbDB" / "png")
                    save_path_pdf = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "LuMet-TCM" / "pdf") if db in ["TCM","ANIMAL"] \
                        else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "HerbDB" / "pdf")
                    plt.savefig(f"{save_path_png}/{pic_name}.{fmt}")
                    plt.savefig(f"{save_path_pdf}/{pic_name}.pdf")
                plt.close()
        # excel 列名修改\drop & overwrite  20230413 del 厂家
        peaks["Fragment Ions"] = pd.Series(fragment_ions)
        peaks = peaks.query("`Fragment Ions` != 'delete'").reset_index(drop=True)
        peaks = peaks[[*peaks.columns[:9], "Fragment Ions", *peaks.columns[9:-1]]]
        peaks.drop(columns=["Compound ID", "厂家", "DB"], inplace=True)
        peaks.rename(columns={"Compound": "ID", "Description": "Metabolites"}, inplace=True)

        # 重新计算峰面积比值
        tcm_cols = self.parse_exp_col()
        tcm_percent_cols = []
        for col in tcm_cols:
            peaks[f"{col}（峰面积的比值%）"] = 100 * peaks[col] / peaks[col].sum()
            tcm_percent_cols.append(f"{col}（峰面积的比值%）")
        # add M mean
        peaks["M均值（峰面积的比值%）"] = peaks[tcm_percent_cols].mean(axis=1)
        # float handle
        peaks["m/z"] = peaks["m/z"].apply(lambda x: f"{x:.4f}")
        peaks["Retention time (min)"] = peaks["Retention time (min)"].apply(lambda x: f"{x:.2f}")
        peaks["Mass Error (ppm)"] = peaks["Mass Error (ppm)"].apply(lambda x: f"{x:.2f}")
        peaks.to_excel(data_path, sheet_name="data", encoding="utf-8", index=False)

    def plot_pie(self, type: str, used_column: str = "M均值（峰面积的比值%）",
                 fmt: str = "png", dpi: int = 200, water_mask: bool = False):
        """
            plot pie of gradients
        :param type: 饼图类型 [count, content, top10]
        :param used_column: 用于定量的列名
        :param fmt: png
        :param dpi: 200
        """
        # 自定义colormap
        from lmbio.basic.colors import SelectColors
        my_map = ListedColormap(colors=list(SelectColors("blindless")))
        logger.info("正在绘制中药代谢物%s的分布饼图...", type)
        """绘制分类饼图"""
        if sys.platform == "linux":
            # Path(__file__).parent / "data" / "fonts" / "Microsoft YaHei.ttf"
            font_path = Path("/usr/share/fonts") / "Microsoft YaHei.ttf"
            font_manager.fontManager.addfont(str(font_path))
        # 读取数据表
        project_name, project_info = self.parse_project()
        data_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "数据矩阵.xlsx") if self._project != "Blood" else str(self.midware / "数据矩阵.xlsx")

        df = pd.read_excel(data_path, sheet_name=0).applymap(lambda x: x.strip() if isinstance(x, str) else x)

        # 兼容入血项目
        used_column = "TCM均值（峰面积的比值%）" if self._project == "Blood" else used_column

        # 删除 Na row
        df.dropna(axis=0, how="any", subset=["中文大类"], inplace=True)
        # plt rcParams
        plt.style.use(style="default")
        plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
        fig, ax = plt.subplots(figsize=(10, 6), dpi=dpi, constrained_layout=True)

        # save_path
        save_path = str(Path(self.data) / project_name / "数据分析" / "5.中药成分分类")
        # variables external
        data = None
        title = None
        if type == "top10":
            save_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果")
            # top10
            df = pd.read_excel(data_path, sheet_name=0)
            # 名称缩短
            df['Metabolites'] = df['Metabolites'].apply(self.filename_handler)
            data = df.sort_values(by=used_column, ascending=False).head(10).set_index(["Metabolites"], drop=True)
            title = "中药成分含量TOP10"
        elif type == "count":
            # 数量占比
            data = df.groupby("中文大类").agg("count")["ID"]
            title = "中药成分分类数量分布图"
        elif type == "count_en":
            # 数量占比 en
            data = df.groupby("英文大类").agg("count")["ID"]
            title = "中药成分分类数量分布图-英文版"
        elif type == "content_en":
            # 含量占比 en
            data = df.groupby("英文大类").agg("sum")[used_column]
            title = "中药成分分类含量分布图-英文版"
        elif type == "content":
            # 含量占比
            data = df.groupby("中文大类").agg("sum")[used_column]
            title = "中药成分分类含量分布图"

        x = data[used_column].to_list() if type == "top10" else data.to_list()
        labels = data.index.to_list()
        patches, ltexts, ntexts = ax.pie(x, labels=labels,
                                         # colors=mpl.cm.get_cmap("tab20")(range(len(x))),
                                         colors=my_map(range(len(x))),
                                         autopct="%1.2f%%", pctdistance=0.9, labeldistance=1.1,
                                         # textprops=dict(horizontalalignment="center")
                                         )

        for l, n in zip(ltexts, ntexts):
            l.set_size(6)
            n.set_size(6)
            if type != "top10":
                # ignore < 1%
                digit = float(n.get_text().replace("%", ''))
                if digit < 1:
                    l.set(text="")
                    n.set(text="")
        # from oebio.plot.adjust_text import adjust_text
        # adjust_text(ltexts, only_move={'text': 'y'},
        #             lim=3,
        #             precision=0.5,
        #             force_text=(0.001, 0.002)
        #             )
        # adjust_text(ltexts)
        ax.set_title(title, fontsize=16)
        if type != "top10":
            ax.legend(loc="lower right",
                      # x,y,w,h
                      bbox_to_anchor=(0.95, 0, 0.4, 1),
                      fontsize=5)
        if water_mask:
            fig.text(0.5, 0.5, 'created by lumingbio', transform=ax.transAxes,
                     fontsize=40, color='gray', alpha=0.5,
                     ha='center', va='center', rotation=30)

        plt.savefig(f"{save_path}/{title}.{fmt}")
        plt.savefig(f"{save_path}/{title}.pdf")
        plt.close()

    def plot_pies(self, water_mask: bool = False):
        self.plot_pie("count", water_mask=water_mask)
        self.plot_pie("count_en", water_mask=water_mask)
        self.plot_pie("content", water_mask=water_mask)
        self.plot_pie("content_en", water_mask=water_mask)
        self.plot_pie("top10", water_mask=water_mask)

    def add_src(self):
        logger.info("正在添加中药成分归属信息...")
        # 读取数据表
        project_name, project_info = self.parse_project()
        data_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "数据矩阵.xlsx") if self._project != "Blood" else str(self.midware / "数据矩阵.xlsx")
        df = pd.read_excel(data_path, sheet_name=0)

        # load db
        with self.ENGINE.connect() as conn:
            hs = pd.read_sql('herb_src', conn)[["Ingredient_INCHIKEY", "Herb_cn_name"]]
        m1 = df.merge(hs, how='left', left_on='InChIKey',
                      right_on="Ingredient_INCHIKEY").drop_duplicates()
        m1.Herb_cn_name = m1.Herb_cn_name.fillna('')

        # link data
        df = df.set_index(keys="InChIKey", drop=False)
        df["来源"] = m1.groupby("InChIKey").apply(lambda x: ", ".join(list(set(filter(None, x["Herb_cn_name"])))))
        # export
        df.to_excel(data_path, index=False)

    def add_src_in_prescript(self, prescript=["麝香", "蟾酥", "牛黄", "冰片", "红参", "三七", "琥珀", "丹参", "苏合香"]):
        logger.info("正在添加中药成分归属信息...")
        # 读取数据表
        project_name, project_info = self.parse_project()
        data_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "数据矩阵.xlsx") if self._project != "Blood" else str(self.midware / "数据矩阵.xlsx")
        df = pd.read_excel(data_path, sheet_name=0)

        # load db
        with self.ENGINE.connect() as conn:
            hs = pd.read_sql('herb_src', conn)[["Ingredient_INCHIKEY", "Herb_cn_name"]]
        m1 = df.merge(hs, how='left', left_on='InChIKey',
                      right_on="Ingredient_INCHIKEY").drop_duplicates()
        m1.Herb_cn_name = m1.Herb_cn_name.fillna('')

        # link data
        df = df.set_index(keys="InChIKey", drop=False)
        df["来源"] = m1.groupby("InChIKey").apply(lambda x: ", ".join([yc for yc in list(set(filter(None, x["Herb_cn_name"]))) if yc in prescript]))
        # export
        df.to_excel(data_path, index=False)

    def oebio_report_deprecated(self, version="/data/hstore1/database/report/2023-03-04~2023/src", company="LM", using_herb=False, cloud=False):

        from lmbio.basic.reportinfo import getreportinfo
        os.environ['OEBIO'] = version
        reportinfo = getreportinfo(sys_logo=company)
        if cloud:
            os.environ['CLOUD'] = "True"
            os.environ['CLOUDUSER'] = "lmbioinformatics@lumingbio.com"
            os.environ['PASSWORD'] = "Lmbio@123"

        logger.info("正在生成中药复方成分鉴定报告...")

        # this part can delete when no bug
        # project_name, project_info = self.parse_project()
        # data_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "数据矩阵.xlsx")
        # peaks = pd.read_excel(data_path, sheet_name=0)
        # # 重新计算峰面积比值
        # tcm_cols = self.parse_exp_col()
        # tcm_percent_cols = []
        # for col in tcm_cols:
        #     peaks[f"{col}（峰面积的比值%）"] = 100 * peaks[col] / peaks[col].sum()
        #     tcm_percent_cols.append(f"{col}（峰面积的比值%）")
        # # add M mean
        # peaks["M均值（峰面积的比值%）"] = peaks[tcm_percent_cols].mean(axis=1)
        # peaks.to_excel(data_path, index=False)
        # this part can delete when no bug -----end

        # src_path = Path(__file__).parent / "data"
        src_path = Path("/data/hstore1/database/database/tcm/")
        project_name, project_info = self.parse_project()

        # 批量重命名
        data_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "数据矩阵.xlsx")
        df = pd.read_excel(data_path, sheet_name=0).applymap(lambda x: x.strip() if isinstance(x, str) else x)
        new_names = pd.Series(data=[f"compound{x:0>5d}" for x in range(1, df.shape[0] + 1)])
        name_map = dict(zip(df["No."], new_names))
        fig_path = Path(self.data) / project_name / "数据分析" / "3.质谱图"
        formats = ("pdf", "png")
        for fm in formats:
            for p in fig_path.rglob(f"*.{fm}"):
                _ = name_map.get(p.stem, None)
                if _:
                    dst = Path(p.parent) / (_ + f".{fm}")
                    os.rename(p, dst)
                else:
                    logger.info(f"文件名{p.stem} 映射出错!")
        df["No."] = new_names
        df.to_excel(data_path, index=False)

        # 加载配置文件
        config = yaml.load(open(str(src_path / "tcm" / "config.yaml"), encoding="utf8"),
                           Loader=yaml.FullLoader)
        # 拷贝文件到实验内容
        report_xls = [path for path in (src_path / "tables").rglob("*.xlsx")] + [src_path / "tcm" / "数据矩阵字段说明.xlsx"]
        for xls in report_xls:
            shutil.copy(str(xls), f"{str(Path(self.data) / project_name / '实验内容')}/{xls.name}")

        # 拷贝UV data
        uv_path = [path for path in (Path(self.data) / "UVdata").rglob("*紫外吸收图*")]
        for uv in uv_path:
            shutil.copy(str(uv), f"{str(Path(self.data) / project_name / '数据分析' / '1.紫外吸收图')}/{uv.name}")

        # copy pic data
        # path = src_path / "luming"
        # shutil.copy(str(path / "LM1.jpg"), str(Path(self.data) / project_name / '鹿明信息' / 'LM1.jpg'))
        # shutil.copy(str(path / "LM2.png"), str(Path(self.data) / project_name / '鹿明信息' / 'LM2.png'))
        # shutil.copy(str(path / "LM3.png"), str(Path(self.data) / project_name / '鹿明信息' / 'LM3.png'))

        # 设置项目信息
        config["header_info"]["项目名称"] = '中药成分鉴定分析'
        # config["header_info"]["客户单位"] = project_info.loc[0, "客户单位"].values
        config["header_info"]["任务单号"] = project_info.loc[0, "项目编号"] if "-b" in project_info.loc[0, "项目编号"] else project_info.loc[0, "项目编号"] + "-b1"
        config["header_info"]["客户名称"] = project_info.loc[0, "客户名称"]
        config["header_info"]["联系人名称"] = project_info.loc[0, "联系人名称"]
        config["header_info"]["项目编号"] = project_info.loc[0, "项目编号"]
        config["header_info"]["样本"] = project_info.loc[0, "样本"]
        # config["header_info"]["完成时间"] = datetime.datetime.now().strftime("%Y-%m")
        # 默认 jiang tao
        # config["header_info"]["执行编码"] = 'LM0460'
        header_info = config.get("header_info")

        # OUT PATH
        report_dir = str(Path(self.data) / project_name)

        # 实例化报告对象
        report = oebioReport('中药成分鉴定分析', title='中药成分鉴定分析', header_info=header_info, oe_welcome = reportinfo.oe_weclome)

        # 从配置文件加载 main section 内容
        # with open(str(src_path / 'herb_template.yaml'), "r", encoding="utf8") as f1:
        #     content = f1.read()
        #     content = content.replace("html_header: 'header.jpg'",
        #                               f"html_header: {repr(str(src_path / 'luming' /'header.jpg'))}")
        #     with open(str(src_path / 'herb.yaml'), "w", encoding="utf8") as f2:
        #         f2.write(content)

        report.add_yaml_config(str(src_path / 'herb.yaml'))

        # 前言
        intro = report.add_section('前言')

        # 实验内容
        exp = report.add_section('实验内容')

        ## 材料
        sub1 = exp.add_section('材料', description=project_info.loc[0, "材料描述信息"]) \
            if isinstance(project_info.loc[0, "材料描述信息"], str) \
            else exp.add_section('材料', description=config['实验内容']['材料']['description'])
        son1 = sub1.add_section('试剂')
        son1.add_table(f"{report_dir}/{config['实验内容']['材料']['试剂']['tablepath']}",
                       caption=config['实验内容']['材料']['试剂']['tablename'])
        son2 = sub1.add_section('仪器')
        son2.add_table(f"{report_dir}/{config['实验内容']['材料']['仪器']['tablepath']}",
                       caption=config['实验内容']['材料']['仪器']['tablename'])

        ## 方法
        sub2 = exp.add_section('方法')
        sub2.add_section('前处理', description=project_info.loc[0, "前处理描述信息"]) \
            if isinstance(project_info.loc[0, "前处理描述信息"], str) \
            else sub2.add_section('前处理', description=config['实验内容']['方法']['前处理']['description'])
        sub2.add_section('液相色谱-质谱条件', description=config['实验内容']['方法']['液相色谱-质谱条件']['description'])
        sub2.add_table(f"{report_dir}/{config['实验内容']['方法']['液相色谱-质谱条件']['tablepath1']}",
                       caption=config['实验内容']['方法']['液相色谱-质谱条件']['tablename1'])
        sub2.add_comment(description=config['实验内容']['方法']['液相色谱-质谱条件']['comment'])
        sub2.add_table(f"{report_dir}/{config['实验内容']['方法']['液相色谱-质谱条件']['tablepath2']}",
                       caption=config['实验内容']['方法']['液相色谱-质谱条件']['tablename2'])

        # 数据分析
        analysis = report.add_section('数据分析')
        a = analysis.add_section('紫外吸收图', description=config["数据分析"]["紫外吸收图"]["description"])
        a.add_fig(f"{report_dir}/{config['数据分析']['紫外吸收图']['imgpath1']}", caption=config['数据分析']['紫外吸收图']['imgname1'])
        a.add_fig(f"{report_dir}/{config['数据分析']['紫外吸收图']['imgpath2']}", caption=config['数据分析']['紫外吸收图']['imgname2'])

        a1 = analysis.add_section('基峰图', description=config["数据分析"]['基峰图']['description'])
        a1.add_fig(f"{report_dir}/{config['数据分析']['基峰图']['imgpath1']}", caption=config['数据分析']['基峰图']['imgname1'])
        a1.add_fig(f"{report_dir}/{config['数据分析']['基峰图']['imgpath2']}", caption=config['数据分析']['基峰图']['imgname2'])

        a2 = analysis.add_section('数据预处理')
        a2.add_section("Progenesis QI v3.0的定性分析",
                       description=config["数据分析"]['数据预处理']['Progenesis QI v3.0的定性分析']['description'])
        son3 = a2.add_section("定性定量结果", description=config["数据分析"]['数据预处理']['定性定量结果']['description'])
        son3.add_table(f"{report_dir}/{config['数据分析']['数据预处理']['定性定量结果']['tablepath1']}",
                       caption=config['数据分析']['数据预处理']['定性定量结果']['tablename1'])
        son3.add_table(f"{report_dir}/{config['数据分析']['数据预处理']['定性定量结果']['tablepath2']}",
                       caption=config['数据分析']['数据预处理']['定性定量结果']['tablename2'], show_search=True)
        son3.add_fig(f"{report_dir}/{config['数据分析']['数据预处理']['中药成分分类']['imgpath0']}",
                     caption=config['数据分析']['数据预处理']['中药成分分类']['imgname0'])
        # MS/MS
        son4 = a2.add_section("中药主要化学成分的EIC图及其与标准品比对情况",
                              description=config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['description'])

        son4.add_fig(f"{config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['imgpath']}",
                     caption=config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['imgname'],
                     path=report_dir)
        # LuMet-TCM
        lumet_tcm = son4.add_section("中药原方成分的EIC图及其与标准品库LuMet-TCM的比对情况")
        lumet_tcm.add_fig(f"{config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与标准品库LuMet-TCM的比对情况']['imgpath']}",
                     caption=config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与标准品库LuMet-TCM的比对情况']['imgname'],
                     path=report_dir)
        if using_herb:
            # HerbDB
            herb_db = son4.add_section("中药原方成分的EIC图及其与公共库HerbDB的比对情况")
            herb_db.add_fig(f"{config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与公共库HerbDB的比对情况']['imgpath']}",
                         caption=config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与公共库HerbDB的比对情况']['imgname'],
                         path=report_dir)
        # pie
        son6 = a2.add_section("中药成分分类", description=config['数据分析']['数据预处理']['中药成分分类']['description'])
        son6.add_fig(f"{report_dir}/{config['数据分析']['数据预处理']['中药成分分类']['imgpath1']}",
                     caption=config['数据分析']['数据预处理']['中药成分分类']['imgname1'], description="为了显示美观，自动隐藏了占比<1%的数据标签")
        son6.add_fig(f"{report_dir}/{config['数据分析']['数据预处理']['中药成分分类']['imgpath2']}",
                     caption=config['数据分析']['数据预处理']['中药成分分类']['imgname2'], description="为了显示美观，自动隐藏了占比<1%的数据标签")

        # 公司简介
        reportinfo.companyinfo(report=report,
                               yamlpath=os.environ['OEBIO'])
        # add extra info
        # company = report.add_section('公司简介')
        # company.add_fig(f"{report_dir}/鹿明信息/LM1.jpg", caption='上海鹿明生物科技有限公司资质')
        # company_instrument = company.add_section('公司仪器')
        # company_instrument.add_fig(f"{report_dir}/鹿明信息/LM2.png", caption='上海鹿明生物科技有限公司仪器平台')
        # company_vision = company.add_section('公司服务')
        # company_vision.add_fig(f"{report_dir}/鹿明信息/LM3.png", caption='上海鹿明生物科技有限公司服务')
        # declare = report.add_section('申明')

        # export
        report.write_to(f"{report_dir}/项目报告.html", zip_report_name=report_dir+f"-OECloud.zip")

        # 添加与报告到 /public/项目报告检查空间/预报告文件夹/
        pnum = Path(self.data).name # project_info.loc[0, "项目编号"]
        pname = project_info.loc[0, "项目名称"]
        # shutil.copy(f"{report_dir}/项目报告.html", f"/public/项目报告检查空间/预报告文件夹/{pnum}-{pname}-LCMS-1.0.0-分析预报告.html")
        shutil.copy(f"{report_dir}/项目报告.html", f"/data/nas/177/代谢/项目报告检查空间/预报告文件夹/{pnum}-{pname}-LCMS-1.0.0-分析预报告.html")
        # del picture
        # os.unlink(str(Path(self.data) / project_name / '数据分析' / 'LM1.jpg'))
        # os.unlink(str(Path(self.data) / project_name / '数据分析' / 'LM2.png'))
        # os.unlink(str(Path(self.data) / project_name / '数据分析' / 'LM3.png'))

    def oebio_report(self, version="/data/hstore1/database/report/2023-03-04~2023/src", company="LM",
                                using_herb=False, cloud=False):
        script_path = str(Path(__file__).parent / "oecloud.py")
        local = "0" if self.local else "1"
        args = ["python", script_path, "-d", self.data, "-t", "tcm", "-l", local]
        if using_herb:
            args.append("--herb")
        if cloud:
            args.append("--cloud")
        r = sp.run(args, stderr=sp.PIPE, stdout=sp.PIPE, text=True, encoding='utf8')
        if r.returncode != 0:
            logger.info(r.stderr)
            logger.info("oebio_report TCM Failed!")
        else:
            logger.info("oebio_report TCM SUCCESS!")

    def zip_report(self):
        project_name, info = self.parse_project()
        zip_in = self.data + os.sep + project_name
        zip_out = zip_in + "_local.zip"
        logger.info("正在生成压缩文件%s...", zip_out)

        zip_file = zipfile.ZipFile(zip_out, 'w', zipfile.ZIP_DEFLATED)
        for dir_path, dir_names, file_names in os.walk(zip_in):
            save_path = dir_path.replace(zip_in, '')
            # nas virus dir
            if Path(dir_path).name not in ["@eaDir", "delete"]:
                for filename in file_names:
                    # nas virus file
                    if filename not in ["Thumbs.db", "delCpds.xlsx", "preprocessed.xlsx"] and \
                            not filename.startswith('~$') and \
                            not filename.endswith(".tmp"):
                        zip_file.write(os.path.join(dir_path, filename), os.path.join(save_path, filename))
        zip_file.close()
        logger.info("压缩文件完成%s...", zip_out)
        # rm zip_in
        # shutil.rmtree(zip_in, ignore_errors=True)

    def run_pipeline(self, using_herb=False, water_mask=False, cloud=True):
        logger.info("正在执行中药成分鉴定分析的完整pipeline...")
        """ whole procedure of this report!"""
        self.raw2mzml()
        self.raw2obs()
        self.preprocess_qi(using_herb=using_herb)
        self.plot_chromatograms(water_mask=water_mask)
        self.plot_xic(dpi=600, img_size=250, water_mask=water_mask)
        self.plot_pies(water_mask=water_mask)
        self.add_src()
        self.oebio_report(using_herb=using_herb, cloud=cloud)
        self.zip_report()
        logger.info("中药成分鉴定分析成功完成...")


    def feedback_all_detected_MS2(self):
        mode_used, sample_info = self.parse_meta()
        project_name, info = self.parse_project()
        data_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "数据矩阵.xlsx") if self._project == "TCM" else \
            str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "入血成分分析数据矩阵.xlsx")
        df = pd.read_excel(data_path, sheet_name=0)

        msp_path = Path(self.data) / "QIdata" / "MSP"
        refMS2: dict = self.loadMS2FromMsp(msp_path, mode_used)

        res = []
        for i, row in df.iterrows():
            _id, ion_mode = row["ID"], row["Ion mode"]
            expMS2 = refMS2[ion_mode].get(_id, "")
            if len(expMS2) > 0:
                exp_MS2 = expMS2[0] if isinstance(expMS2, list) else expMS2
                mz_str = ", ".join([str(e) for e in exp_MS2.mz])
                res.append(mz_str)
            else:
                res.append('')
        df["MS2"] = pd.Series(res)
        df.to_excel(data_path, index=False)

    def test(self):
        from oebio.report2 import thumbnail
        p1 = Path("/data/hstore4/lumingos/project/DZLM2023050885-1/DZLM2023050885-薛娟-中药入血成分分析(项目报告)/数据分析/3.质谱图/中药原方成分/LuMet-TCM/png")
        p2 = Path("/data/hstore4/lumingos/project/DZLM2023050885-1/DZLM2023050885-薛娟-中药入血成分分析(项目报告)/数据分析/3.质谱图/中药原方成分/HerbDB/png")

        #thumbnail(p1)
        for f in p1.rglob('*.png'):
            #f = str(file)
            try:
                thumbnail(f)
            except Exception as e:
                print(e)
                print("p1", f)
        for f in p2.rglob('*.png'):
            # f = str(file)
            try:
                thumbnail(f)
            except Exception as e:
                print(e)
                print("p2", f)

if __name__ == '__main__':
    # report = Report(data=r"D:\pycharm_project\pyopenms\report\cm_gredient_report\data2")
    # report.plot_tic()
    # print(Report.get_smiles(MSSpectrum(),Report.ENGINE,"eds","sd","dsa","ds",[]))
    ...
