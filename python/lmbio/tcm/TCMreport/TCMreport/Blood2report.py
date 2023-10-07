# File: Blood2report
# Author: jiang tao
# Time: 2023/2/17 10:57
import functools
import logging
import os

import shutil
import zipfile
from glob import glob
from pathlib import Path
from typing import Tuple, List, Union

import xlsxwriter
from tqdm import tqdm

import numpy as np
import pandas as pd
import subprocess as sp
import yaml

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.ticker import MultipleLocator
from pyopenms import MSExperiment, MzMLFile, MSSpectrum
from sqlalchemy import create_engine, text


from .QI2report import Report
from .network import plot_networkx
from .Constants import _engine
from oebio.report import Report as oebioReport

logger = logging.getLogger("TCMreport")


class BloodReport(Report):
    """入血成分报告类"""

    @staticmethod
    def write_sheet(df, worksheet, _format, header_format):
        """
        :param df: data matrix
        :param worksheet:
        :param _format:
        :param header_format:
        :return:
        """
        # write header
        for ncol, value in enumerate(df.columns.values):
            worksheet.write(0, ncol, value, header_format)
        # write body
        for idx, row in enumerate(df.iterrows()):
            for ncol in range(df.shape[1]):
                worksheet.write(idx + 1, ncol, row[1][ncol], _format)

    @staticmethod
    def prepare_table(df, file: str, ncomp: tuple):
        # workbook
        workbook = xlsxwriter.Workbook(file, {'nan_inf_to_errors': True})
        _format = workbook.add_format({'font_name': 'calibri', 'font_size': 11})
        header_format = workbook.add_format({'font_name': 'calibri', 'font_size': 11, 'bold': True})
        # worksheet
        # add sheet dataframe
        worksheet1 = workbook.add_worksheet("数据矩阵")
        BloodReport.write_sheet(df, worksheet1, _format, header_format)
        # add sheet group
        worksheet2 = workbook.add_worksheet("分组")
        sample = pd.Series(data=[f"Treatment-{i}" for i in range(1, 1+ncomp[0])] + [f"Control-{i}" for i in range(1, 1+ncomp[0])], name="Sample")
        exp = pd.Series(data=["Treatment"] * ncomp[0] + [""] * ncomp[1], name="Treatment")
        ctl = pd.Series(data=[""] * ncomp[0] + ["Control"] * ncomp[1], name="Control")
        group = pd.concat([sample, exp, ctl], axis=1)
        BloodReport.write_sheet(group, worksheet2, _format, header_format)
        # add sheet compare
        worksheet3 = workbook.add_worksheet("比较")
        compare = pd.DataFrame(data=[("Treatment/Control", "FALSE")], columns=["比较组", "配对"])
        BloodReport.write_sheet(compare, worksheet3, _format, header_format)
        workbook.close()

    def fill_qi_df(self, data_processed: List[pd.DataFrame], dedup=True) -> pd.DataFrame:
        """post_process steps in preprocess_qi"""
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
            db = res.loc[i, "DB"]
            cid = res.loc[i, "Compound ID"]
            if db == "PRODUCT":
                sql = self.sql1 % cid
                res_proxy = conn.execute(text(sql)).fetchall()
                if len(res_proxy) != 0:
                    inchikey, smiles = res_proxy[0][0], res_proxy[0][1]
                    res.loc[i, "InChIKey"] = inchikey
                    res.loc[i, "SMILES"] = smiles
                    continue
            # 1.identifiers
            sql = self.sql2 % cid
            res_proxy = conn.execute(text(sql)).fetchall()

            if len(res_proxy) != 0:
                inchikey, smiles, _,  hmdbid, metlin, lipidmaps, kegg, pubchem, cas = res_proxy[0][0], res_proxy[0][1], res_proxy[0][2], \
                                                                                      res_proxy[0][3], res_proxy[0][4], res_proxy[0][5], \
                                                                                      res_proxy[0][6], res_proxy[0][7], res_proxy[0][8]
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
        if dedup:
            # del duplicated cpd and select best
            res = res.sort_values("Compound ID").reset_index(drop=True)
            best_index = self.selectIndex(res, "Compound ID")
            res = res.loc[best_index, :].sort_values("Compound").reset_index(drop=True)
        return res

    def fill_qi_df1(self, data_processed: List[pd.DataFrame], dedup=True) -> pd.DataFrame:
        """only for lc-ms"""
        # concat res in two ion modes
        res = pd.concat(data_processed).sort_values("Compound").reset_index(drop=True)

        # add new columns
        res["InChIKey"] = ""
        res["SMILES"] = ""
        res["HMDB"] = ""
        res["METLIN"] = ""
        res["LipidMaps"] = ""
        res["KEGG"] = ""
        res["ChEBI"] = ""
        res["PubChem"] = ""
        res["CAS"] = ""
        res["Super Class"] = ""
        res["Class"] = ""
        res["Sub Class"] = ""
        res["cid"] = ""

        # assign values to new columns from database
        conn = self.ENGINE.connect()
        for i in range(res.shape[0]):
            ikey = res.loc[i, "Compound ID"]
            # 1.identifiers
            sql = self.sql4 % ikey
            res_proxy = conn.execute(text(sql)).fetchall()

            if len(res_proxy) != 0:
                inchikey, smiles, hmdbid, metlin, lipidmaps, kegg, chebi, pubchem, cas,\
                    supercls, cls, subcls, cid = res_proxy[0][0], res_proxy[0][1], res_proxy[0][2], \
                                            res_proxy[0][3], res_proxy[0][4], res_proxy[0][5], \
                                            res_proxy[0][6], res_proxy[0][7], res_proxy[0][8], \
                                            res_proxy[0][9], res_proxy[0][10], res_proxy[0][11], res_proxy[0][12]
                res.loc[i, "InChIKey"] = ikey
                res.loc[i, "SMILES"] = smiles
                res.loc[i, "HMDB"] = hmdbid
                res.loc[i, "METLIN"] = metlin
                res.loc[i, "LipidMaps"] = lipidmaps
                res.loc[i, "KEGG"] = kegg
                res.loc[i, "ChEBI"] = chebi
                res.loc[i, "PubChem"] = pubchem
                res.loc[i, "CAS"] = cas
                res.loc[i, "Super Class"] = supercls
                res.loc[i, "Class"] = cls
                res.loc[i, "Sub Class"] = subcls
                res.loc[i, "cid"] = cid
        conn.close()
        if dedup:
            # del duplicated cpd and select best
            res = res.sort_values("Compound ID").reset_index(drop=True)
            best_index = self.selectIndex(res, "Compound ID")
            res = res.loc[best_index, :].sort_values("Compound").reset_index(drop=True)
        return res

    @functools.lru_cache(maxsize=1000)
    def parse_meta(self) -> Tuple[list, pd.DataFrame]:
        """override parse sample info"""
        mode_used = []

        xls_path = Path(self.data) / "sampleinfo.xlsx"
        sample_info = pd.read_excel(str(xls_path), sheet_name="样品信息", index_col=[0, 1]).T

        is_null = sample_info.isnull().any(axis=1)
        if not is_null["POS"]:
            mode_used.append("POS")
        if not is_null["NEG"]:
            mode_used.append("NEG")

        return mode_used, sample_info

    @functools.lru_cache(maxsize=1000)
    def parse_exp_col(self):
        """
        返回样本的实际列名
        """
        mode_used, sample_info = self.parse_meta()
        mode = mode_used[0]
        blood_cols = sample_info["Blood"].drop(columns=["QI定性结果", "QI定量结果", "QI二级MSP"]).loc[mode, :].apply(
            lambda x: x.rsplit(".", 1)[0]).tolist()
        control_cols = sample_info["Control"].loc[mode, :].apply(lambda x: x.rsplit(".", 1)[0]).tolist()
        tcm_cols = sample_info["TCM"].loc[mode, :].apply(lambda x: x.rsplit(".", 1)[0]).tolist()

        return blood_cols, control_cols, tcm_cols

    def preprocess_qi(self, FC=10, treshold_ctl=5000, min_treshold_exp=1000, treshold_tcm=2500, score_exp=48,
                      peak_width: tuple = (0.075, 0.40),
                      export_unknown=False, using_herb=False, using_product=False, interested_cpd=None, peak_id=None):
        """override preprocess_qi func
            对TCM和Blood组进行QI result合并过滤，取交集
        """
        unknown_bloods = None
        logger.info("正在预处理QI结果，添加代谢物注释...")
        logger.info(f"差异倍数： {FC}")
        mode_used, sample_info = self.parse_meta()

        MIDWARE_PATH = Path(self.data) / "midware"
        MIDWARE_PATH.mkdir(exist_ok=True)

        blood_cols, control_cols, tcm_cols = self.parse_exp_col()
        ncomp = len(blood_cols), len(control_cols)
        # calc n
        n = len(blood_cols + control_cols + tcm_cols) + 1

        data_processed_blood = []
        data_processed_auto = []

        if export_unknown:
            logger.info("导出unknown入血成分...")
            data_unknown_blood = []
            for mode in mode_used:
                src = Path(self.data)
                # BLOOD GROUP
                m_path_blood = src / "QIdata" / mode / sample_info["Blood"].loc[mode, "QI定量结果"]
                quantitative_blood = pd.read_csv(m_path_blood, encoding='utf-8', header=[2, ])
                quantitative_blood = quantitative_blood.drop(columns=[*quantitative_blood.columns[-(10 + n):]])
                # unknown
                quantitative_blood = quantitative_blood.query("Identifications == 0")
                # filter
                # 排除 0值
                # M_TCM = quantitative_blood.columns[-1]
                # res_blood = quantitative_blood.query(f"{M_TCM} > 0").reset_index(drop=True)

                # 不管 M_TCM 有没有，全部导出
                res_blood = quantitative_blood.reset_index(drop=True)

                # 扣除空白 FC>=4
                ctl_blood = res_blood[control_cols].mean(axis=1)
                exp_blood = res_blood[blood_cols].mean(axis=1)
                unknown = res_blood.loc[(ctl_blood <= treshold_ctl) & (exp_blood >= min_treshold_exp) &
                                        (exp_blood / ctl_blood >= FC), :].reset_index(drop=True)
                unknown["Ion mode"] = mode
                data_unknown_blood.append(unknown)
            unknown_bloods = pd.concat(data_unknown_blood, ignore_index=True)
        interested_records = []
        for mode in mode_used:
            src = Path(self.data)
            # BLOOD GROUP
            id_path_blood = src / "QIdata" / mode / sample_info["Blood"].loc[mode, "QI定性结果"]
            m_path_blood = src / "QIdata" / mode / sample_info["Blood"].loc[mode, "QI定量结果"]

            # ISO-8859-1
            qualitative_blood = pd.read_csv(id_path_blood, encoding='utf-8')
            quantitative_blood = pd.read_csv(m_path_blood, encoding='utf-8', header=[2, ])

            qualitative_blood = qualitative_blood.drop(
                columns=["Accepted?", "Link", "Neutral mass (Da)", "m/z", "Charge", "Retention time (min)",
                         "Isotope Similarity", "Theoretical Isotope Distribution"])
            # TODO：QC列写死了
            quantitative_blood = quantitative_blood.drop(
                columns=[*quantitative_blood.columns[-(10 + n):], "QC", "Neutral mass (Da)",
                          "Chromatographic peak width (min)",
                         "Identifications", "Isotope Distribution", "Maximum Abundance", "Minimum CV%"])

            qual_quant_blood = pd.merge(qualitative_blood, quantitative_blood, on="Compound", how="inner")

            # 峰宽过滤
            # qual_quant_blood = qual_quant_blood.query(
            #     f"`Chromatographic peak width (min)` >= {peak_width[0]} & `Chromatographic peak width (min)` <= {peak_width[1]}")
            # qual_quant_blood.drop(columns=["Chromatographic peak width (min)"], inplace=True)

            # filter with different rules when dealing with various db
            # Description & add DB_ID
            qual_quant_blood["DB"] = qual_quant_blood["Description"].apply(lambda x: x.split(":", 1)[0])
            qual_quant_blood["Description"] = qual_quant_blood["Description"].apply(lambda x: x.split(":", 1)[1])
            # todo: 临时保留代谢产物
            # qual_quant_blood = qual_quant_blood[qual_quant_blood["DB"] == "PRODUCT"]
            # rm prefix PRODUCT: TODO: only do this when compound id with prefix 'product:'
            # if using_product:
            #     qual_quant_blood["Compound ID"] = qual_quant_blood["Compound ID"].apply(
            #                                 lambda x: x.split("PRODUCT:", 1)[1] if x.__contains__('PRODUCT:') else x)
            # add annotation level
            qual_quant_blood["Level"] = qual_quant_blood.apply(
                lambda x: "level4" if x["DB"] == "LuMet-Animal-LCMS-Auto" else "level3" if x["DB"] == "HERB" else "level1" if x["Fragmentation Score"] >= 50 else "level2",
                axis=1)
            # TODO: retain interested compounds
            if interested_cpd:
                _tmp = qual_quant_blood.query("Description in @interested_cpd")
                if peak_id:
                    _nmap = dict(zip(interested_cpd, peak_id))
                    for i, row in _tmp.iterrows():
                        if _nmap.get(row["Description"]) != row["Compound"]:
                            _tmp.drop(index=i, inplace=True)
                if _tmp.shape[0]:
                    _tmp["Ion mode"] = mode
                    _tmp["TCM均值"] = _tmp[tcm_cols].mean(axis=1)
                    _tmp["TCM均值（峰面积的比值%）"] = 1
                interested_records.append(_tmp)
            # LC-MS-Process
            qual_quant_auto = qual_quant_blood.query("DB == 'LuMet-Animal-LCMS-Auto'")
            qual_quant_auto = qual_quant_auto.query("Score >= 36")

            # TCM
            qual_quant_tcm = qual_quant_blood.query("DB == 'TCM'")
            # qual_quant_tcm = qual_quant_tcm.query("Score >= 50 | (Score >= 40 & `Fragmentation Score` >= 50)")
            qual_quant_tcm = qual_quant_tcm.query("Score >= 40")
            # Animal
            qual_quant_animal = qual_quant_blood.query("DB == 'ANIMAL'")
            qual_quant_animal = qual_quant_animal.query(
                "(Score >= 50 & `Fragmentation Score` >= 50) | (Score >= 40 & `Fragmentation Score` >= 50)")
            # Product
            qual_quant_product = qual_quant_blood.query("DB == 'PRODUCT'")
            # qual_quant_product = qual_quant_product.query("Score >= 50 | (Score >= 40 & `Fragmentation Score` >= 50)")
            qual_quant_product = qual_quant_product.query("Score >= 40")
            # filter with rules and Sort
            if using_herb:
                # Herb
                qual_quant_herb = qual_quant_blood.query("DB == 'HERB'")
                qual_quant_herb = qual_quant_herb.query("Score >= 50 & `Fragmentation Score` >= 50")

                qual_quant_blood = pd.concat([qual_quant_tcm, qual_quant_animal, qual_quant_herb]).sort_values(
                    "Compound").reset_index(drop=True)
            else:
                qual_quant_blood = pd.concat([qual_quant_tcm, qual_quant_animal]).sort_values(
                    "Compound").reset_index(drop=True)
            if using_product:
                    qual_quant_blood = pd.concat([qual_quant_blood, qual_quant_product]).sort_values(
                        "Compound").reset_index(drop=True)

            best_index1 = self.selectIndex(qual_quant_blood, "Compound")
            qual_quant_blood = qual_quant_blood.loc[best_index1, :]

            qual_quant_blood["Ion mode"] = mode
            # add M mean
            qual_quant_blood["TCM均值"] = qual_quant_blood[tcm_cols].mean(axis=1)
            qual_quant_blood["TCM均值（峰面积的比值%）"] = 100 * qual_quant_blood["TCM均值"] / qual_quant_blood["TCM均值"].sum()

            data_processed_blood.append(qual_quant_blood)

            # lc-ms
            if len(qual_quant_auto) != 0:
                qual_quant_auto = qual_quant_auto.sort_values("Compound").reset_index(drop=True)
                best_index2 = self.selectIndex(qual_quant_auto, "Compound")
                qual_quant_auto = qual_quant_auto.loc[best_index2, :]
                qual_quant_auto["Ion mode"] = mode
                qual_quant_auto["TCM均值"] = qual_quant_auto[tcm_cols].mean(axis=1)
                qual_quant_auto["TCM均值（峰面积的比值%）"] = 100 * qual_quant_auto["TCM均值"] / qual_quant_auto["TCM均值"].sum()
                data_processed_auto.append(qual_quant_auto)

        res_blood = self.fill_qi_df(data_processed_blood)
        res_blood["theoretical m/z"] = (1e6 * res_blood["m/z"]) / (1e6 + res_blood["Mass Error (ppm)"])

        # lc-ms--------------------------------------------------------------------------------
        res_auto = None
        rel_cols = ["Metabolites"]
        if len(data_processed_auto) != 0:
            res_auto = self.fill_qi_df1(data_processed_auto)
            res_auto["theoretical m/z"] = (1e6 * res_auto["m/z"]) / (1e6 + res_auto["Mass Error (ppm)"])
            res_auto = res_auto[
                ["Compound", "Compound ID", "cid", "Adducts", "Formula", "Score", "Fragmentation Score", "Mass Error (ppm)",
                 "Description", "theoretical m/z", "m/z", "Retention time (min)", "Ion mode",
                 *control_cols, *blood_cols, #*tcm_cols, "TCM均值", "TCM均值（峰面积的比值%）",
                 "InChIKey", "SMILES", "HMDB",
                 "METLIN", "LipidMaps", "KEGG", "ChEBI", "PubChem", "CAS","Super Class", "Class", "Sub Class"]]
            # lc-ms level
            res_auto.loc[:, "level"] = 'Level 4'
            res_auto.loc[res_auto["Fragmentation Score"] > 45, "level"] = 'Level 3'
            res_auto["rtscore"] = res_auto["Score"] - res_auto["Fragmentation Score"]/5 -40
            res_auto.loc[res_auto["rtscore"] > 0, "level"] = "Level 2"
            res_auto.loc[(res_auto["rtscore"] > 0) & (res_auto["Fragmentation Score"] > 45), "level"] = "Level 1"
            res_auto.drop(columns=["rtscore", "Compound ID"], inplace=True)
            res_auto["Super Class"] = res_auto["Super Class"].fillna("Unclassified")
            res_auto["Class"] = res_auto["Class"].fillna("Unclassified")
            res_auto["Sub Class"] = res_auto["Sub Class"].fillna("Unclassified")

            # rename
            ren_dict = {}
            new_tcm_cols = []
            res_auto.columns = ["ID", "cid", "Adducts", "Formula", "Score", "Fragmentation Score", "Mass Error (ppm)",
                 "Metabolites", "theoretical m/z", "m/z", "Retention time (min)", "Ion mode",
                 *control_cols, *blood_cols, #*tcm_cols, "TCM均值", "TCM均值（峰面积的比值%）",
                                "InChIKey", "SMILES", "HMDB",
                 "METLIN", "LipidMaps", "KEGG", "ChEBI", "PubChem", "CAS","Super Class", "Class", "Sub Class", "level"]
            for i, ctl in enumerate(control_cols, start=1):
                ren_dict[ctl] = f"Control-{i}"
                rel_cols.append(f"Control-{i}")
            for i, exp in enumerate(blood_cols, start=1):
                ren_dict[exp] = f"Treatment-{i}"
                rel_cols.append(f"Treatment-{i}")
            if len(tcm_cols) == 1:
                ren_dict[tcm_cols[0]] = "TCM"
                new_tcm_cols.append("TCM")
            else:
                for i, tcm in enumerate(tcm_cols, start=1):
                    ren_dict[tcm] = f"TCM-{i}"
                    new_tcm_cols.append(f"TCM-{i}")
            res_auto = res_auto.sort_values(by="level", ascending=True)

            res_auto.rename(columns=ren_dict, inplace=True)
            res_auto = res_auto[~(res_auto.loc[:, rel_cols[1:]] == 0).all(axis=1)].drop_duplicates(subset="Metabolites").reset_index(drop=True)
        # lc-ms --------------------------------------------------------------------------------

        if interested_records:
            if not pd.concat(interested_records).empty:
                interested_records= self.fill_qi_df(interested_records, dedup=False)
                interested_records["theoretical m/z"] = (1e6 * interested_records["m/z"]) / (1e6 + interested_records["Mass Error (ppm)"])
                ctl_blood = interested_records[control_cols].mean(axis=1)
                exp_blood = interested_records[blood_cols].mean(axis=1)
                interested_records["是否为入血成分"] = "否"
                interested_records.loc[(interested_records["DB"] != "PRODUCT") &
                              (interested_records["TCM均值"] >= treshold_tcm) &
                              (ctl_blood <= treshold_ctl) &
                              (exp_blood >= min_treshold_exp) &
                              (exp_blood / ctl_blood >= FC) &
                              (interested_records["Score"] >= score_exp), "是否为入血成分"] = "是"
                interested_records["No."] = pd.Series(data=[f"cpdInterested{x:0>3d}" for x in range(1, interested_records.shape[0] + 1)])
                interested_records = interested_records[
                    ["No.", "Compound", "Compound ID", "Adducts", "Formula", "Score", "Fragmentation Score",
                     "Mass Error (ppm)",
                     "Description", "theoretical m/z", "m/z", "Retention time (min)", "Ion mode", "是否为入血成分",
                     *control_cols, *blood_cols, *tcm_cols, "TCM均值", "TCM均值（峰面积的比值%）", "InChIKey", "SMILES", "HMDB",
                     "METLIN", "LipidMaps", "KEGG", "PubChem", "CAS", "中文名", "中文大类", "中文子类", "英文大类", "英文子类", "厂家", "货号",
                     "纯度",
                     "DB", "Level"]]

                interested_records["lower_name"] = interested_records['Description'].apply(lambda x: x.lower())
                interested_records = interested_records.sort_values(by="lower_name").reset_index(drop=True)
                indexOfBestName = self.selectIndexOfBestName(interested_records, 'lower_name')
                interested_records = interested_records.loc[indexOfBestName, :].reset_index(drop=True).applymap(
                    lambda x: x.strip() if isinstance(x, str) else x)
                interested_records.sort_values(by="No.").to_excel(str(self.midware / "interested_cpd.xlsx"), index=False)
        # 排除 0值
        # res_blood = res_blood[res_blood["TCM均值"] > treshold_tcm]
        # 扣除空白 FC>=4
        ctl_blood = res_blood[control_cols].mean(axis=1)
        exp_blood = res_blood[blood_cols].mean(axis=1)
        res_blood["是否为入血成分"] = "否"

        # 原型入血成分
        res_blood.loc[(res_blood["DB"] != "PRODUCT") &
                      (res_blood["TCM均值"] >= treshold_tcm) &
                      (ctl_blood <= treshold_ctl) &
                      (exp_blood >= min_treshold_exp) &
                      (exp_blood / ctl_blood >= FC) &
                      (res_blood["Score"] >= score_exp), "是否为入血成分"] = "是"
        # 入血代谢产物
        self._product : bool = True if using_product else False
        if using_product:
            # TODO:按照入血成分标准来卡
            res_blood.loc[(res_blood["DB"] == "PRODUCT") &
                          (res_blood["TCM均值"] < treshold_tcm) &
                          (ctl_blood <= treshold_ctl) &
                          (exp_blood >= min_treshold_exp) &
                          (exp_blood / ctl_blood >= FC)
                          & (res_blood["Score"] >= score_exp)
                          ,"是否为入血成分"] = "是"
            # 过滤掉 不符合要求的代谢产物
            res_blood = res_blood.loc[(res_blood["DB"] != "PRODUCT") | ((res_blood["DB"] == "PRODUCT") & (res_blood["是否为入血成分"] == "是"))]
        res_blood_is = res_blood.query("`是否为入血成分` == '是'")

        # 中药原方
        res_blood_not = res_blood.query("`是否为入血成分` == '否'")
        res_blood_not = res_blood_not[(res_blood_not["TCM均值"] > treshold_tcm) & (res_blood_not["DB"] != "PRODUCT")]

        res_blood = pd.concat([res_blood_is, res_blood_not], ignore_index=True)
        # 根据名称去重
        res_blood["lower_name"] = res_blood['Description'].apply(lambda x: x.lower())
        res_blood = res_blood.sort_values(by="lower_name").reset_index(drop=True)
        indexOfBestName = self.selectIndexOfBestName(res_blood, 'lower_name')
        res_blood = res_blood.loc[indexOfBestName, :].reset_index(drop=True).applymap(lambda x:x.strip() if isinstance(x, str) else x)
        #  add cpd number
        res_blood["No."] = pd.Series(data=[f"compound{x:0>4d}" for x in range(1, res_blood.shape[0] + 1)])
        # 去列+排序
        res_blood = res_blood[
            ["No.", "Compound", "Compound ID", "Adducts", "Formula", "Score", "Fragmentation Score", "Mass Error (ppm)",
             "Description", "theoretical m/z", "m/z", "Retention time (min)", "Ion mode", "是否为入血成分",
             *control_cols, *blood_cols, *tcm_cols, "TCM均值", "TCM均值（峰面积的比值%）", "InChIKey", "SMILES", "HMDB",
             "METLIN", "LipidMaps", "KEGG", "PubChem", "CAS", "中文名", "中文大类", "中文子类", "英文大类", "英文子类", "厂家", "货号", "纯度",
             "DB", "Level"]]

        # merge & filter
        # export
        dirname, project_info = self.parse_project()
        # 删除动物+lipid, 如果中药是完全植物源
        is_del_animal = False if project_info.loc[0, "复方中是否包含动物"] == '是' else True
        if is_del_animal:
            with self.ENGINE.connect() as tx:
                df = pd.read_sql('herb_annotation', tx)
                dels = df[(df.is_animal_source.astype(int) + df.is_lipid.astype(int)) >= 1]["inchikey"].to_list()
            res_blood = res_blood.query('`Compound ID` not in @dels').reset_index(drop=True)
        res_blood = res_blood.sort_values(by=["是否为入血成分", "Level"], ascending=[False, True]).reset_index(drop=True)
        SAVE_PATH = Path(self.data) / dirname

        SAVE_PATH.mkdir(exist_ok=True)

        # (SAVE_PATH / "鹿明信息").mkdir(exist_ok=True)
        (SAVE_PATH / "实验内容").mkdir(exist_ok=True)
        (SAVE_PATH / "数据分析").mkdir(exist_ok=True)
        (SAVE_PATH / "数据分析" / "1.紫外吸收图").mkdir(exist_ok=True)
        (SAVE_PATH / "数据分析" / "2.色谱图").mkdir(exist_ok=True)
        (SAVE_PATH / "数据分析" / "3.质谱图").mkdir(exist_ok=True)
        (SAVE_PATH / "数据分析" / "3.质谱图" / "中药入血成分").mkdir(exist_ok=True)
        (SAVE_PATH / "数据分析" / "3.质谱图" / "中药入血成分" / "中药原型" / "png").mkdir(exist_ok=True, parents=True)
        (SAVE_PATH / "数据分析" / "3.质谱图" / "中药入血成分" / "中药原型" / "pdf").mkdir(exist_ok=True, parents=True)
        (SAVE_PATH / "数据分析" / "3.质谱图" / "中药原方成分").mkdir(exist_ok=True, parents=True)
        (SAVE_PATH / "数据分析" / "3.质谱图" / "中药原方成分" / "LuMet-TCM" / "png").mkdir(exist_ok=True, parents=True)
        (SAVE_PATH / "数据分析" / "3.质谱图" / "中药原方成分" / "LuMet-TCM" / "pdf").mkdir(exist_ok=True, parents=True)
        if using_herb:
            (SAVE_PATH / "数据分析" / "3.质谱图" / "中药原方成分" / "HerbDB" / "png").mkdir(exist_ok=True, parents=True)
            (SAVE_PATH / "数据分析" / "3.质谱图" / "中药原方成分" / "HerbDB" / "pdf").mkdir(exist_ok=True, parents=True)
        if using_product:
            (SAVE_PATH / "数据分析" / "3.质谱图" / "中药入血成分" / "中药代谢产物" / "png").mkdir(exist_ok=True, parents=True)
            (SAVE_PATH / "数据分析" / "3.质谱图" / "中药入血成分" / "中药代谢产物" / "pdf").mkdir(exist_ok=True, parents=True)
            (SAVE_PATH / "数据分析" / "6.代谢网络图" / "png").mkdir(exist_ok=True, parents=True)
            (SAVE_PATH / "数据分析" / "6.代谢网络图" / "pdf").mkdir(exist_ok=True, parents=True)

        (SAVE_PATH / "数据分析" / "4.定性定量结果").mkdir(exist_ok=True)
        (SAVE_PATH / "数据分析" / "5.中药成分分类").mkdir(exist_ok=True)

        res_blood.to_excel(str(MIDWARE_PATH / "数据矩阵.xlsx"), sheet_name="data", encoding="utf-8",
                           index=False)
        res_blood.to_excel(str(MIDWARE_PATH / "preprocessed.xlsx"), sheet_name="data", encoding="utf-8",
                           index=False)
        # lc-ms
        if res_auto is not None:
            # res_auto[rel_cols].to_excel(str(MIDWARE_PATH / "relations.xlsx"), index=False)
            self.prepare_table(res_auto, str(MIDWARE_PATH / "lcms-preprocessed.xlsx"), ncomp)
            self.diff_analysis()
        if export_unknown:
            unknown_bloods.to_excel(str(SAVE_PATH / "数据分析" / "4.定性定量结果" / "数据矩阵-未知入血成分.xlsx"), sheet_name="data",
                                    encoding="utf-8",
                                    index=False)

    def get_max_intensity(self, used_exps: List[MSExperiment], _mz: float, mass_tolerance: float):
        """deprecated methods"""
        res = []
        for exp in used_exps:
            intensities = []
            for spec in exp:
                if spec.getMSLevel() == 1:
                    filtered_int = []
                    for mz, i in zip(*spec.get_peaks()):
                        if _mz - mass_tolerance < mz < _mz + mass_tolerance:
                            filtered_int.append(i)
                    ints = sum(filtered_int)
                    intensities.append(ints)
            intensities = np.asarray(intensities)
            res.append(intensities.max())
        return max(res)

    def plot_eic(self, ax: Axes, exp: MSExperiment, _mz: float, mass_tolerance: float, label: str, color: str,
                 is_blood: bool=False):
        """the real method of plot EIC"""
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

        # 设置 轴刻度
        ax.xaxis.set_major_locator(MultipleLocator(4))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        if not is_blood:
            # ax.yaxis.set_major_locator(MultipleLocator(20))
            # ax.yaxis.set_minor_locator(MultipleLocator(5))
            # ax.axis([0, 16, 0, 100])
            # ax.plot(retention_times, intensities, lw=0.8, label=label, color=color)
            # ax.set_ylabel('Relative Abundance')
            ax.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
            ax.set_ylabel('Abundance')
            ax.set_xlim([0, 16])
            ax.plot(retention_times, intensities, lw=0.8, label=label, color=color)
            if intensities.max() <= 1000:
                ax.set_ylim([0, 1000])
            else:
                yticks = ax.get_yticks()
                ax.set_ylim([0, yticks[-1]])
        else:
            ax.set_xlim([0, 16])
            # TODO: align y labels vertically
            ax.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
            if intensities.max() <= 1000:
                ax.set_ylim([0, 1000])
            # ax.spines['bottom'].set_position(('data', 0))
            # ax.spines['left'].set_position(('data', 0))
            ax.plot(retention_times, intensities, lw=0.8, label=label, color=color)
            ax.set_ylabel('Abundance')
            yticks = ax.get_yticks()
            ax.set_ylim([0, yticks[-1]])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        if is_blood:
            ax.legend(fontsize=6)
        return retention_times, intensities * (100 / intensities.max())

    def query_plot_mirror(self, adduct, refMS2, ion_mode, _id, inchikey, formula, name, ax1, _mz, img_size, db):
        """deprecated at 2023/5/26
            为了兼容中药代谢产物，在 get_parent_smiles, readMS2FromDB, get_smiles中添加 is_product 参数
            与old_version的plot_msms配套
        """
        is_product = True if db == "PRODUCT" else False

        adducts = adduct.split(", ")
        expMS2 = refMS2[ion_mode].get(_id, "")
        # 前体一定能查到!
        precursor = self.get_parent_smiles(inchikey, is_product=is_product)
        fragment_ions_str = None
        # USING HERB
        if db == "HERB":
            # herb库无二级，plot_defect
            self.plot_defect(ax1, precursor, _mz, img_size=img_size)
        else:
            if len(expMS2) > 0:
                exp_MS2 = expMS2[0] if isinstance(expMS2, list) else expMS2
                exp_spec = MSSpectrum()
                exp_spec.set_peaks((exp_MS2.mz, exp_MS2.intensity))
                exp_spec = self.spec_normalize(exp_spec)

                ref_spec = self.readMS2FromDB(inchikey, ion_mode.lower(), adducts, exp_spec=exp_spec,
                                              is_product=is_product)

                if ref_spec:
                    ref_spec = self.spec_normalize(ref_spec)
                    _smiles: Union[None, list] = self.get_smiles(ref_spec, self.ENGINE, formula, inchikey, name,
                                                                 ion_mode, adducts, is_product=is_product)
                    if _smiles:
                        # ref在上面并绘制smiles
                        fragment_ions_str = self.plot_mirror_spectrum_with_smiles(ax1, exp_spec, ref_spec, _smiles,
                                                                                  precursor, _mz,
                                                                                  img_size=img_size)
                    else:
                        fragment_ions_str = self.plot_mirror_spectrum_with_smiles(ax1, exp_spec, ref_spec, [],
                                                                                  precursor,
                                                                                  _mz,
                                                                                  img_size=img_size)
                else:
                    if is_product:
                        self.plot_defect(ax1, precursor, _mz, img_size=img_size)
                    else:
                        fragment_ions_str = "delete"

            else:
                ref_spec = self.readMS2FromDB(inchikey, ion_mode.lower(), adducts, is_product=is_product)
                if ref_spec:
                    self.plot_defect(ax1, precursor, _mz, img_size=img_size)
                else:
                    if is_product:
                        self.plot_defect(ax1, precursor, _mz, img_size=img_size)
                    else:
                        fragment_ions_str = "delete"
        return fragment_ions_str

    def plot_msms(self, project,
                  _id, inchikey, adduct, formula, name, _mz, _rt, ion_mode, is_blood,
                  used_exps, refMS2, ppm,
                  labels, colormap, img_size, dpi, fmt, db, water_mask: bool = False):
        """deprecated at 2023/5/26"""
        fragment_ions_str = None
        if is_blood:
            fig, axes = plt.subplots(4, 1, figsize=(8, 11.4), dpi=dpi, gridspec_kw={'height_ratios': [1, 1, 1, 2.7]})
            # 设置白底
            fig.patch.set_alpha(1.)
            if water_mask:
                fig.text(0.5, 0.5, 'created by lumingbio',
                         fontsize=40, color='gray', alpha=0.5,
                         ha='center', va='center', rotation=30)
            ax0 = axes[0:3]
            ax1 = axes[3]
            exps = used_exps

            _rt, _mz = float(_rt), float(_mz)
            mass_tolerance = (_mz / 1e6) * ppm

            for i, ax in enumerate(ax0):
                self.plot_eic(ax, exps[i], _mz, mass_tolerance, labels[i], colormap[i], is_blood)

            title = self.filename_handler(name, escape=False)
            ax0[0].set_title(f"{title}: Mode_{ion_mode}  RT_{_rt:.2f} min  MZ_{_mz:.4f}")

            # 3.2 query data & plot mirror or MS2
            fragment_ions_str = self.query_plot_mirror(adduct, refMS2, ion_mode, _id, inchikey, formula, name, ax1, _mz,
                                                       img_size, db)
            if fragment_ions_str != 'delete':
                # save fig
                pic_name = self.filename_handler(name, escape=True)
                save_path_png = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "png")
                save_path_pdf = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "pdf")
                plt.savefig(f"{save_path_png}/{pic_name}.{fmt}")
                plt.savefig(f"{save_path_pdf}/{pic_name}.pdf")
            plt.close()
        else:
            fig, axes = plt.subplots(2, 1, figsize=(8, 7.4), dpi=dpi, gridspec_kw={'height_ratios': [1, 2.7]})
            # 设置白底
            fig.patch.set_alpha(1.)
            if water_mask:
                fig.text(0.5, 0.5, 'created by lumingbio',
                         fontsize=40, color='gray', alpha=0.5,
                         ha='center', va='center', rotation=30)

            ax0 = axes[0]
            ax1 = axes[1]
            exps = used_exps
            _rt, _mz = float(_rt), float(_mz)
            mass_tolerance = (_mz / 1e6) * ppm

            self.plot_eic(ax0, exps[2], _mz, mass_tolerance, labels[2], colormap[0], is_blood)

            title = self.filename_handler(name, escape=False)
            ax0.set_title(f"{title}: Mode_{ion_mode}  RT_{_rt:.2f} min  MZ_{_mz:.4f}")

            # 3.2 query data & plot mirror or MS2
            fragment_ions_str = self.query_plot_mirror(adduct, refMS2, ion_mode, _id, inchikey, formula, name, ax1, _mz,
                                                       img_size, db)
            if fragment_ions_str != 'delete':
                # save fig
                pic_name = self.filename_handler(name, escape=True)
                save_path_png = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "png")
                save_path_pdf = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "pdf")
                plt.savefig(f"{save_path_png}/{pic_name}.{fmt}")
                plt.savefig(f"{save_path_pdf}/{pic_name}.pdf")
            plt.close()
        return fragment_ions_str if fragment_ions_str else ''

    def query_plot_spec(self, adduct, refMS2, ion_mode, _id, inchikey, formula, name, ax1, _mz, img_size, db):
        """
            为了兼容中药代谢产物，在 get_parent_smiles, readMS2FromDB, get_smiles中添加 is_product 参数
            2023/5/26: 与 newversion的 plot_ms2配套
                TCM+Animal 处理方法相同 [会额外剔除标准品库中无二级，而实测有二级的物质]
                HERB+PRODUCT 处理相同
        """
        is_product = True if db == "PRODUCT" else False

        adducts = adduct.split(", ")
        expMS2 = refMS2[ion_mode].get(_id, "")
        precursor = self.get_parent_smiles(inchikey, is_product=is_product)
        fragment_ions_str = None
        # logic
        if len(expMS2) > 0:
            exp_MS2 = expMS2[0] if isinstance(expMS2, list) else expMS2
            exp_spec = MSSpectrum()
            exp_spec.set_peaks((exp_MS2.mz, exp_MS2.intensity))
            exp_spec = self.spec_normalize(exp_spec)

            ref_spec = self.readMS2FromDB(inchikey, ion_mode.lower(), adducts, exp_spec=exp_spec,
                                          is_product=is_product)
            if ref_spec:
                # TODO: 在HERB MS2库出来之前临时添加，后续需要删除
                if db == "HERB":
                    # herb找到std二级，说明RT不对，需要删除
                    fragment_ions_str = "delete"
                else:
                    ref_spec = self.spec_normalize(ref_spec)
                    _smiles: Union[None, list] = self.get_smiles(ref_spec, self.ENGINE, formula, inchikey, name,
                                                                 ion_mode, adducts, is_product=is_product)
                    if _smiles:
                        # ref在上面并绘制smiles
                        fragment_ions_str = self.plot_mirror_spectrum_with_smiles(ax1, exp_spec, ref_spec, _smiles,
                                                                                      precursor, _mz,
                                                                                      img_size=img_size)
                    else:
                        fragment_ions_str = self.plot_mirror_spectrum_with_smiles(ax1, exp_spec, ref_spec, [],
                                                                                  precursor,
                                                                                  _mz,
                                                                                  img_size=img_size)
            else:
                # std无MS2，实测有 MS2
                if db in ["TCM", "ANIMAL"]:
                    fragment_ions_str = "delete"
                else:
                    # 非标品库绘制
                    _smiles: Union[None, list] = self.get_smiles(exp_spec, self.ENGINE, formula, inchikey, name,
                                                                 ion_mode, adducts, is_product=is_product)
                    if _smiles:
                        fragment_ions_str = self.plot_spectrum_with_smiles(ax1, exp_spec, _smiles, precursor, _mz,
                                                                           img_size=img_size)
                    else:
                        fragment_ions_str = self.plot_spectrum(ax1, exp_spec, precursor, _mz, img_size=img_size)
        # 未检测到二级
        else:
            self.plot_defect(ax1, precursor, _mz, img_size=img_size)
            fragment_ions_str = ""
        return fragment_ions_str if fragment_ions_str else ""

    def query_plot_spec1(self, adduct, refMS2, ion_mode, _id, inchikey, formula, name, ax1, _mz, img_size, db):
        """
            only for add interested cpd ...
            为了兼容中药代谢产物，在 get_parent_smiles, readMS2FromDB, get_smiles中添加 is_product 参数
            2023/5/26: 与 newversion的 plot_ms2配套
                TCM+Animal 处理方法相同 [会额外剔除标准品库中无二级，而实测有二级的物质]
                HERB+PRODUCT 处理相同
        """
        is_product = True if db == "PRODUCT" else False

        adducts = adduct.split(", ")
        expMS2 = refMS2[ion_mode].get(_id, "")
        precursor = self.get_parent_smiles(inchikey, is_product=is_product)
        fragment_ions_str = None
        # logic
        if len(expMS2) > 0:
            exp_MS2 = expMS2[0] if isinstance(expMS2, list) else expMS2
            exp_spec = MSSpectrum()
            exp_spec.set_peaks((exp_MS2.mz, exp_MS2.intensity))
            exp_spec = self.spec_normalize(exp_spec)

            ref_spec = self.readMS2FromDB(inchikey, ion_mode.lower(), adducts, exp_spec=exp_spec,
                                          is_product=is_product)
            if ref_spec:
                # TODO: 在HERB MS2库出来之前临时添加，后续需要删除
                ref_spec = self.spec_normalize(ref_spec)
                _smiles: Union[None, list] = self.get_smiles(ref_spec, self.ENGINE, formula, inchikey, name,
                                                             ion_mode, adducts, is_product=is_product)
                if _smiles:
                    # ref在上面并绘制smiles
                    fragment_ions_str = self.plot_mirror_spectrum_with_smiles(ax1, exp_spec, ref_spec, _smiles,
                                                                                  precursor, _mz,
                                                                                  img_size=img_size)
                else:
                    fragment_ions_str = self.plot_mirror_spectrum_with_smiles(ax1, exp_spec, ref_spec, [],
                                                                              precursor,
                                                                              _mz,
                                                                              img_size=img_size)
            else:
                # std无MS2，实测有 MS2
                #if db in ["TCM", "ANIMAL"]:
                # 非标品库绘制
                _smiles: Union[None, list] = self.get_smiles(exp_spec, self.ENGINE, formula, inchikey, name,
                                                             ion_mode, adducts, is_product=is_product)
                if _smiles:
                    fragment_ions_str = self.plot_spectrum_with_smiles(ax1, exp_spec, _smiles, precursor, _mz,
                                                                       img_size=img_size)
                else:
                    fragment_ions_str = self.plot_spectrum(ax1, exp_spec, precursor, _mz, img_size=img_size)
        # 未检测到二级
        else:
            self.plot_defect(ax1, precursor, _mz, img_size=img_size)
            fragment_ions_str = ""
        return fragment_ions_str if fragment_ions_str else ""

    def plot_ms2(self, cpdn, project,
                 _id, inchikey, adduct, formula, name, _mz, _rt, ion_mode, is_blood,
                 used_exps, refMS2, ppm,
                 labels, colormap, img_size, dpi, fmt, db, water_mask: bool = False):

        """2023/5/26： 根据新的需求对 plot_msms 方法进行修改， 后续在plot_xic中使用此方法"""

        fragment_ions_str = None
        is_blood_plot_save = False
        if is_blood:
            fig, axes = plt.subplots(4, 1, figsize=(8, 11.4), dpi=dpi, gridspec_kw={'height_ratios': [1, 1, 1, 2.7]})
            # 设置白底
            fig.patch.set_alpha(1.)
            if water_mask:
                fig.text(0.5, 0.5, 'created by lumingbio',
                         fontsize=40, color='gray', alpha=0.5,
                         ha='center', va='center', rotation=30)
            ax0 = axes[0:3]
            ax1 = axes[3]
            exps = used_exps

            _rt, _mz = float(_rt), float(_mz)
            mass_tolerance = (_mz / 1e6) * ppm

            for i, ax in enumerate(ax0):
                self.plot_eic(ax, exps[i], _mz, mass_tolerance, labels[i], colormap[i], is_blood)

            title = self.filename_handler(name, escape=False)
            ax0[0].set_title(f"{title}: Mode_{ion_mode}  RT_{_rt:.2f} min  MZ_{_mz:.4f}")

            # 3.2 query data & plot mirror or MS2
            fragment_ions_str = self.query_plot_spec(adduct, refMS2, ion_mode, _id, inchikey, formula, name, ax1, _mz,
                                                     img_size, db)
            if fragment_ions_str != 'delete':
                is_blood_plot_save = True
                # save fig
                pic_name = cpdn
                save_path_png = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药原型" / "png") if db != "PRODUCT" else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药代谢产物" / "png")
                save_path_pdf = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药原型" / "pdf") if db != "PRODUCT" else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药代谢产物" / "pdf")
                plt.savefig(f"{save_path_png}/{pic_name}.{fmt}")
                plt.savefig(f"{save_path_pdf}/{pic_name}.pdf")
            plt.close()

        # 非入血成分全部绘制
        # BUG: 无法提前计算出图形的长宽比 8,5.2   1,1.6
        if db == "PRODUCT":
            # 代谢产物不需画TCM图
            return fragment_ions_str
        fig, axes = plt.subplots(2, 1, figsize=(8, 7.4), dpi=dpi, gridspec_kw={'height_ratios': [1, 2.7]})
        # 设置白底
        fig.patch.set_alpha(1.)
        if water_mask:
            fig.text(0.5, 0.5, 'created by lumingbio',
                     fontsize=40, color='gray', alpha=0.5,
                     ha='center', va='center', rotation=30)

        ax0 = axes[0]
        ax1 = axes[1]

        exps = used_exps
        _rt, _mz = float(_rt), float(_mz)
        mass_tolerance = (_mz / 1e6) * ppm
        # rt in min, relative its in 0-100
        rts, its = self.plot_eic(ax0, exps[2], _mz, mass_tolerance, labels[2], colormap[0])
        # new: check ms1 20230623
        is_ms1_valid = self.check_ms1(rts, its, _rt)

        # calc
        title = self.filename_handler(name, escape=False)
        ax0.set_title(f"{title}: Mode_{ion_mode}  RT_{_rt:.2f} min  MZ_{_mz:.4f}")
        # todo: tmp save bad peaks
        if not is_ms1_valid:
            # pic_name = self.filename_handler(name, escape=True, full=True)
            # save_path_png = Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "delete"
            # if not save_path_png.exists():
            #     save_path_png.mkdir(parents=True, exist_ok=True)
            # # ax0.axhline(noise, color='red')
            # # ax0.text(_rt + 0.1, noise + 1, f'Noise={noise}')
            # ax0.axvline(_rt, color='green')
            # plt.savefig(f"{str(save_path_png)}/{pic_name}.{fmt}")
            fragment_ions_str = 'delete'
        else:
            # 3.2 query data & plot mirror or MS2
            fragment_ions_str = self.query_plot_spec(adduct, refMS2, ion_mode, _id, inchikey, formula, name, ax1, _mz,
                                                     img_size, db)
        # save fig
        pic_name = cpdn
        if fragment_ions_str != 'delete':
            save_path_png = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "LuMet-TCM" / "png") if db in ["TCM", "ANIMAL"] \
                            else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "HerbDB" / "png")
            save_path_pdf = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "LuMet-TCM" / "pdf") if db in ["TCM", "ANIMAL"] \
                            else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "HerbDB" / "pdf")

            plt.savefig(f"{save_path_png}/{pic_name}.{fmt}")
            plt.savefig(f"{save_path_pdf}/{pic_name}.pdf")
        else:
            if is_blood_plot_save:
                save_path_png = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药原型" / "png") if db != "PRODUCT" else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药代谢产物" / "png")
                save_path_pdf = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药原型" / "pdf") if db != "PRODUCT" else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药代谢产物" / "pdf")
                os.unlink(f"{save_path_png}/{pic_name}.{fmt}")
                os.unlink(f"{save_path_pdf}/{pic_name}.pdf")
        plt.close()
        return fragment_ions_str if fragment_ions_str else ''

    def plot_ms2_nocheck(self, cpdn, project,
                 _id, inchikey, adduct, formula, name, _mz, _rt, ion_mode, is_blood,
                 used_exps, refMS2, ppm,
                 labels, colormap, img_size, dpi, fmt, db, water_mask: bool = False):

        """2023/5/26： 根据新的需求对 plot_msms 方法进行修改， 后续在plot_xic中使用此方法"""

        fragment_ions_str = None
        is_blood_plot_save = False
        if is_blood:
            fig, axes = plt.subplots(4, 1, figsize=(8, 11.4), dpi=dpi, gridspec_kw={'height_ratios': [1, 1, 1, 2.7]})
            # 设置白底
            fig.patch.set_alpha(1.)
            if water_mask:
                fig.text(0.5, 0.5, 'created by lumingbio',
                         fontsize=40, color='gray', alpha=0.5,
                         ha='center', va='center', rotation=30)
            ax0 = axes[0:3]
            ax1 = axes[3]
            exps = used_exps

            _rt, _mz = float(_rt), float(_mz)
            mass_tolerance = (_mz / 1e6) * ppm

            for i, ax in enumerate(ax0):
                self.plot_eic(ax, exps[i], _mz, mass_tolerance, labels[i], colormap[i], is_blood)

            title = self.filename_handler(name, escape=False)
            ax0[0].set_title(f"{title}: Mode_{ion_mode}  RT_{_rt:.2f} min  MZ_{_mz:.4f}")

            # 3.2 query data & plot mirror or MS2
            fragment_ions_str = self.query_plot_spec(adduct, refMS2, ion_mode, _id, inchikey, formula, name, ax1, _mz,
                                                     img_size, db)
            if fragment_ions_str != 'delete':
                is_blood_plot_save = True
                # save fig
                pic_name = cpdn
                save_path_png = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药原型" / "png") if db != "PRODUCT" else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药代谢产物" / "png")
                save_path_pdf = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药原型" / "pdf") if db != "PRODUCT" else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药代谢产物" / "pdf")
                plt.savefig(f"{save_path_png}/{pic_name}.{fmt}")
                plt.savefig(f"{save_path_pdf}/{pic_name}.pdf")
            plt.close()

        # 非入血成分全部绘制
        # BUG: 无法提前计算出图形的长宽比 8,5.2   1,1.6
        if db == "PRODUCT":
            # 代谢产物不需画TCM图
            return fragment_ions_str
        fig, axes = plt.subplots(2, 1, figsize=(8, 7.4), dpi=dpi, gridspec_kw={'height_ratios': [1, 2.7]})
        # 设置白底
        fig.patch.set_alpha(1.)
        if water_mask:
            fig.text(0.5, 0.5, 'created by lumingbio',
                     fontsize=40, color='gray', alpha=0.5,
                     ha='center', va='center', rotation=30)

        ax0 = axes[0]
        ax1 = axes[1]

        exps = used_exps
        _rt, _mz = float(_rt), float(_mz)
        mass_tolerance = (_mz / 1e6) * ppm
        # rt in min, relative its in 0-100
        rts, its = self.plot_eic(ax0, exps[2], _mz, mass_tolerance, labels[2], colormap[0])
        # new: check ms1 20230623
        is_ms1_valid = True

        # calc
        title = self.filename_handler(name, escape=False)
        ax0.set_title(f"{title}: Mode_{ion_mode}  RT_{_rt:.2f} min  MZ_{_mz:.4f}")
        # todo: tmp save bad peaks
        if not is_ms1_valid:
            # pic_name = self.filename_handler(name, escape=True, full=True)
            # save_path_png = Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "delete"
            # if not save_path_png.exists():
            #     save_path_png.mkdir(parents=True, exist_ok=True)
            # # ax0.axhline(noise, color='red')
            # # ax0.text(_rt + 0.1, noise + 1, f'Noise={noise}')
            # ax0.axvline(_rt, color='green')
            # plt.savefig(f"{str(save_path_png)}/{pic_name}.{fmt}")
            fragment_ions_str = 'delete'
        else:
            # 3.2 query data & plot mirror or MS2
            fragment_ions_str = self.query_plot_spec1(adduct, refMS2, ion_mode, _id, inchikey, formula, name, ax1, _mz,
                                                     img_size, db)
        # save fig
        pic_name = cpdn
        if fragment_ions_str != 'delete':
            save_path_png = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "LuMet-TCM" / "png") if db in ["TCM", "ANIMAL"] \
                            else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "HerbDB" / "png")
            save_path_pdf = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "LuMet-TCM" / "pdf") if db in ["TCM", "ANIMAL"] \
                            else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药原方成分" / "HerbDB" / "pdf")

            plt.savefig(f"{save_path_png}/{pic_name}.{fmt}")
            plt.savefig(f"{save_path_pdf}/{pic_name}.pdf")
        else:
            if is_blood_plot_save:
                save_path_png = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药原型" / "png") if db != "PRODUCT" else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药代谢产物" / "png")
                save_path_pdf = str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药原型" / "pdf") if db != "PRODUCT" else str(Path(self.data) / project / "数据分析" / "3.质谱图" / "中药入血成分" / "中药代谢产物" / "pdf")
                os.unlink(f"{save_path_png}/{pic_name}.{fmt}")
                os.unlink(f"{save_path_pdf}/{pic_name}.pdf")
        plt.close()
        return fragment_ions_str if fragment_ions_str else ''

    def plot_xic(self, ppm: int = 5, fmt: str = "png", water_mask: bool = False, interested_cpd=None, **kwargs):
        dpi = 600  # 200
        img_size = 250  # 300
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
        labels = self.group
        colormap = plt.get_cmap('Set1')(range(3))
        # new
        blood_cols, control_cols, tcm_cols = self.parse_exp_col()

        blood_fs = list(map(lambda x: x + ".mzML", blood_cols))
        control_fs = list(map(lambda x: x + ".mzML", control_cols))
        tcm_fs = list(map(lambda x: x + ".mzML", tcm_cols))

        cb_cols = blood_cols + control_cols
        mzml_files = control_fs + blood_fs + tcm_fs
        for mode in mode_used:
            # TODO: 确保 mzML PATH 中没有中文, 否则会 error.
            MSExps[mode] = []
            files = []
            for file in mzml_files:
                abs_f = Path(self.data) / "mzmldata" / mode / file
                files.append(abs_f)
            # load mzml
            for file in files:
                exp = MSExperiment()
                MzMLFile().load(str(file), exp)
                MSExps[mode].append(exp)
        # 2.MS2 container for ref MS/MS
        msp_path = Path(self.data) / "QIdata" / "MSP"
        refMS2: dict = self.loadMS2FromMsp(msp_path, mode_used)

        # 3.iter DataFrame of metabolites
        data_path = str(self.midware / "数据矩阵.xlsx") if not interested_cpd else str(self.midware / "interested_cpd.xlsx")
        # load peak table
        peaks = pd.read_excel(data_path, sheet_name=0)
        # plot each line
        fragment_ions = []
        # progress bar
        with tqdm(peaks.iterrows(), total=peaks.shape[0], desc="INFO - Plot MS2") as pb:
            for i, row in pb:
                cpdn, _id, inchikey, adduct, formula, name, _mz, _rt, ion_mode, is_blood, db = \
                    row["No."], row["Compound"], row["Compound ID"], row["Adducts"], row["Formula"], row["Description"], row["m/z"], \
                    row["Retention time (min)"], row["Ion mode"], row["是否为入血成分"], row["DB"]
                # logger
                pb.set_postfix(compound=name)
                controls = np.asarray(row[control_cols].tolist())
                bloods = np.asarray(row[blood_cols].tolist())
                tcms = np.asarray(row[tcm_cols].tolist())

                index_control = np.argmin(controls)
                index_blood = np.argmax(bloods) + len(controls)
                index_tcm = np.argmax(tcms) + len(cb_cols)

                used_exp = [MSExps[ion_mode][index_control], MSExps[ion_mode][index_blood], MSExps[ion_mode][index_tcm]]
                # 3.1 extract data from mzML & plot XIC
                # pyplot top plot
                is_blood = True if is_blood == "是" else False
                is_product = True if db == "PRODUCT" else False
                # print(self.is_cpd_valid(inchikey, is_product=is_product))
                if self.is_cpd_valid(inchikey, is_product=is_product):
                    func = self.plot_ms2 if not interested_cpd else self.plot_ms2_nocheck
                    fragment_ions_str = func(cpdn, project,
                                                       _id, inchikey, adduct, formula, name, _mz, _rt, ion_mode, is_blood,
                                                       used_exp, refMS2, ppm,
                                                       labels, colormap, img_size, dpi, fmt, db, water_mask=water_mask)
                else:
                    fragment_ions_str = 'delete'
                # print(fragment_ions_str)
                fragment_ions.append(fragment_ions_str)
        peaks["Fragment Ions"] = pd.Series(fragment_ions)
        peaks = peaks.query("`Fragment Ions` != 'delete'").reset_index(drop=True)
        # float handle
        peaks["m/z"] = peaks["m/z"].apply(lambda x: f"{x:.4f}")
        peaks["Retention time (min)"] = peaks["Retention time (min)"].apply(lambda x: f"{x:.2f}")
        peaks["Mass Error (ppm)"] = peaks["Mass Error (ppm)"].apply(lambda x: f"{x:.2f}")
        peaks.to_excel(data_path, index=False)

    def rename_excel(self, interested_cpd=None):
        logger.info("正在处理入血成分组的零值，修改数据矩阵的列名...")

        project_name, project_info = self.parse_project()
        blood_cols, control_cols, tcm_cols = self.parse_exp_col()
        data_path = str(self.midware / "数据矩阵.xlsx") if not interested_cpd else str(self.midware / "interested_cpd.xlsx")

        # load peak table
        peaks = pd.read_excel(data_path, sheet_name=0).reset_index(drop=True)

        for i in range(peaks.shape[0]):
            is_blood = peaks.loc[i, "是否为入血成分"]
            # 0值填充
            if is_blood == "是":
                ctrls = peaks.loc[i, control_cols].tolist()
                controls = np.asarray(ctrls, np.float64)
                if controls.mean() <= 1000:
                    for col in control_cols:
                        peaks.loc[i, col] = 0
        # excel 列名修改
        ren_dict = {"Compound": "ID", "Description": "Metabolites", }
        peaks.rename(columns=ren_dict, inplace=True)
        # force to 100
        peaks["TCM均值（峰面积的比值%）"] = 100 * peaks["TCM均值"] / peaks["TCM均值"].sum()

        # sort
        peaks.sort_values("是否为入血成分", ascending=False, inplace=True)
        # 重排 -1 Fragment Ions
        idx = peaks.columns.get_loc("Metabolites")
        peaks = peaks[[*peaks.columns[:idx], peaks.columns[-1], *peaks.columns[idx:-1]]]
        # del
        peaks.drop(columns=["厂家"], inplace=True)
        peaks.to_excel(data_path, sheet_name="data", encoding="utf-8", index=False)

    def network_analysis(self):
        # do this step after process qi
        logger.info("正在绘制代谢网络图...")

        project_name, project_info = self.parse_project()
        data_path = str(self.midware / "数据矩阵.xlsx")

        # load peak table
        df = pd.read_excel(data_path, sheet_name=0)
        # filter product
        peaks = df[df.DB == "PRODUCT"].reset_index(drop=True)
        # 来自同一母体的代谢物分组 按照绘制网络图
        peaks["group"] = peaks.apply(
            lambda x: x["Compound ID"].split("_", 1)[0] if x.DB == "PRODUCT" else x["Description"], axis=1)
        peaks = peaks.set_index(["group", "Description"]).sort_index()
        peaks["Mid"] = ""
        peaks["SMILES"] = ""
        peaks["Parent Compound"] = ""
        peaks["Transformations"] = ""
        peaks["Metabolism Type"] = ""

        # load parent info
        engine = create_engine(_engine)
        conn = engine.connect()
        ha = pd.read_sql("herb_annotation", conn).set_index(['compound_name']).rename(
            columns={"smiles": "SMILES", "inchikey": "InChIKey", "inchi": "InChI", "rt": "RT"})

        # read herb_rules
        rules = pd.read_sql("herb_rules", conn).set_index(['reaction'])

        # iter cpd group
        with tqdm(peaks.index.levels[0], total=peaks.index.levels[0].size, desc="INFO - Plot Network") as pb:
            for g in pb:
                pb.set_postfix(cpd=g)
                sub = peaks.loc[g]
                if (sub["DB"] == "PRODUCT").any():
                    # prepare network input
                    parent_name = sub.index[0].split("_", 1)[0]
                    p = ha.loc[parent_name, ["SMILES", "InChIKey", "InChI", "RT"]]
                    p = p.iloc[0, :] if isinstance(p, pd.DataFrame) else p
                    parent = pd.DataFrame(data=p.to_dict(), index=[-1])
                    db_idx = sub.columns.get_loc('DB')
                    _ = tuple(sub.iloc[1:, db_idx].index) if sub.iloc[0, db_idx] != "PRODUCT" else tuple(
                        sub.iloc[:, db_idx].index)
                    if len(_) == 1:
                        _ = f"""("{_[0]}")"""
                    sql = text(f"""select * from herb_network WHERE network_id in {_}""")
                    products = pd.read_sql(sql, conn).drop(
                        columns=["pid", "ion_mode", "parent_smiles", "parent_inchi", "parent_inchikey", "parent_name",
                                 "batch"]) \
                        .rename(
                        columns={"select": "筛选", "ReactionClass": "Reaction Class", "Intermediate1": "Intermediate 1",
                                 "Intermediate2": "Intermediate 2"})
                    data = pd.concat([products, parent]).sort_index()
                    # plot func
                    cpd_name = self.filename_handler(parent_name, full=True)
                    save_path = str(Path(self.data) / project_name / "数据分析" / "6.代谢网络图")
                    new_id: pd.DataFrame = plot_networkx(data, save_path, cpd_name)
                    #
                    ids = new_id.Mid.apply(lambda x: x.split('-')[0])
                    clear_cpd: pd.Series = ids[~ids.duplicated(keep=False)]
                    # met = peaks.index.levels[1]
                    for i, row in new_id.iterrows():
                        if i != 0:
                            # if parent_name in met:
                            #     peaks.loc[(g, parent_name), "Mid"] = "M0"
                            # else:
                            # ["network_id"], row["Mid"], row["Reaction Class"], row["SMILES"]
                            net_id, mid, rc, smiles = row
                            mid_fixed = mid.split("-")[0]
                            mid = mid_fixed if mid_fixed in clear_cpd.values else mid
                            peaks.loc[(g, net_id), "Mid"] = mid
                            peaks.loc[(g, net_id), "Parent Compound"] = parent_name
                            df.loc[df["Compound ID"] == net_id, "Description"] = parent_name + "_" + mid
                            peaks.loc[(g, net_id), "SMILES"] = smiles
                            # 处理反应类型 并添加 Transformations & Metabolism Type
                            transformations = []
                            _type = []
                            for r in rc.split("+", -1):
                                _ = r.strip()
                                en, phase = rules.loc[_, ["en", "phase"]]

                                transformations.append(en)
                                _type.append(phase)
                            peaks.loc[(g, net_id), "Transformations"] = ",".join(transformations)
                            peaks.loc[(g, net_id), "Metabolism Type"] = ",".join(sorted(set(_type)))

        conn.close()
        peaks.reset_index(drop=False).sort_values(by=["group", "Mid"]).drop(columns=["group", "Description"]) \
            .to_excel(str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "代谢产物数据矩阵.xlsx"), index=False)
        df.to_excel(data_path, index=False)

    def diff_analysis(self):
        logger.info("start diff_analysis!")
        project_path, project_info = self.parse_project()
        cwd = str(Path(self.data).cwd())
        if "物种" not in project_info.columns:
            raise IndexError("column `物种` was not found in Sampleinfo.xlsx!")
        species = project_info.loc[0, "物种"]
        dest = f"{project_path}/数据分析/7.差异分析结果" if self._product else f"{project_path}/数据分析/6.差异分析结果"
        args = ["flow_diff_Analyze",
                    "-pr", dest,
                    "-om", "M",
                    "-or", species,
                    "-rf", "midware/lcms-preprocessed.xlsx",
                    "-ff", "0",
                    "-vf", "1",
                    "-pf", "0.05"
                ]
        r = sp.run(args, stderr=sp.PIPE, stdout=sp.PIPE, cwd=cwd, encoding='utf8', text=True)
        if r.returncode != 0:
            logger.info("diff_analysis is failed!")
            logger.info(r.stderr)
        logger.info("diff_analysis is completed!")
        if cwd.split(os.sep)[1] == "data":
            if os.path.exists(os.path.join(cwd, ".snakemake")):
                shutil.rmtree(os.path.join(cwd, ".snakemake"))
            if os.path.exists(os.path.join(cwd, "oecloud")):
                shutil.rmtree(os.path.join(cwd, "oecloud"))
        if len(list(glob(f"{cwd}/{project_path}/数据分析/*差异分析结果/report.html"))):
            os.unlink(list(glob(f"{cwd}/{project_path}/数据分析/*差异分析结果/report.html"))[0])

    def relation_analysis(self):
        def table_handler(dataframe):
            dataframe = dataframe.drop_duplicates(subset="Metabolites").reset_index(drop=True)
            remove_zeros = dataframe[~(dataframe.iloc[:, 1:] == 0).all(axis=1)]
            return remove_zeros

        project_path, project_info = self.parse_project()
        if len(list(glob(f"{str(Path(project_path).absolute())}/数据分析/*差异分析结果"))):
            logger.info("relation_analysis")
            cols = ["Metabolites"]
            cols_old = ["Description"]
            blood_cols, control_cols, tcm_cols = self.parse_exp_col()
            for i, ctl in enumerate(control_cols, start=1):
                cols_old.append(ctl)
                cols.append(f"Control-{i}")
            for i, exp in enumerate(blood_cols, start=1):
                cols_old.append(exp)
                cols.append(f"Treatment-{i}")
            data_path = str(self.midware / "数据矩阵.xlsx")

            df = pd.read_excel(data_path, sheet_name=0)
            df = df.query("`是否为入血成分` == '是'")[cols_old]
            if len(df) == 0:
                logger.info("No blood cpds: relation_analysis FAILED!")
                return
            df.columns = cols
            df = table_handler(df)
            used = df["Metabolites"].to_list()
            df.to_excel(str(self.midware / "bloodIngredient.xlsx"), index=False)


            diff_path = list(glob(f"{str(Path(project_path).absolute())}/数据分析/*差异分析结果/*差异代谢物/*差异表达矩阵.xlsx"))[0]
            diff = pd.read_excel(diff_path, sheet_name=0)[cols]
            diff = table_handler(diff)
            diff = diff.query("Metabolites not in @used")
            diff.to_excel(str(self.midware / "diff.xlsx"), index=False)

            dest = list(glob(f"{str(Path(project_path).absolute())}/数据分析/*差异分析结果/"))[0]
            dest = os.path.join(dest, "5.入血代谢物与差异代谢物相关性分析")

            args = ["map_common_corrnetwork",
                    "-f",  str(self.midware / "diff.xlsx"),
                    "-fy", str(self.midware / "bloodIngredient.xlsx"),
                    "-xn", "Endogenous",
                    "-yn", "Ingrediet",
                    "-pf", "0.01",
                    "-cf", "0.99",
                    "-s", dest,
                    "-mn", "入血代谢物与差异代谢物相关性网络"]

            r = sp.run(args, stdout=sp.PIPE, stderr=sp.PIPE, text=True, encoding="utf8")
            if r.returncode != 0:
                logger.info(r.stderr)
                logger.info("relation_analysis FAILED!")
            else:
                logger.info("relation_analysis SUCCESS!")
                if os.path.exists(os.path.join(dest, "说明.txt")):
                    os.unlink(os.path.join(dest, "说明.txt"))

    def extend_network_df(self):
        """补齐所有候选smiles, 导出所有候选"""
        project_name, project_info = self.parse_project()
        data_path = str(self.midware / "raw代谢产物数据矩阵.xlsx")
        df = pd.read_excel(data_path, sheet_name=0)

        df["rns"] = df["Compound ID"].apply(lambda x: x.split('_')[1])
        df["select"] = df["Compound ID"].apply(lambda x: '候选' if x.split('_')[2].__contains__('候选') else '是')
        df["rt"] = df["Compound ID"].apply(lambda x: x.split('_')[3])
        df = df.drop_duplicates(subset=["Parent Compound", "rns", "select", "rt"])

        sql = """SELECT SMILES from herb_network WHERE parent_name="%s" and ReactionClass="%s" and RT="%s" and `select` LIKE '候选_'"""
        conn = self.ENGINE.connect()
        for i, row in df.iterrows():
            if row['select'].strip() != "是":
                _ = text(sql % (row["Parent Compound"], row["rns"], row["rt"]))
                r = pd.read_sql(_, conn)
                if r.shape[0] !=0:
                    smi = ", ".join(r["SMILES"].to_list())
                    df.loc[i, "SMILES"] = smi
        conn.close()
        df.rename(columns={"Compound": "ID"}, inplace=True)
        col_order = ['ID', 'Mid', 'Parent Compound', 'Transformations', 'Adducts', 'Formula', 'Score',
                     'Fragmentation Score',
                     'Mass Error (ppm)', 'm/z',
                     'Retention time (min)', "Metabolites", 'Ion mode', 'Level', 'SMILES', 'Fragment Ions',
                     'Metabolism Type']
        new_products = df[col_order]
        # float handle
        new_products["m/z"] = new_products["m/z"].apply(lambda x: f"{x:.4f}")
        new_products["Retention time (min)"] = new_products["Retention time (min)"].apply(lambda x: f"{x:.2f}")
        new_products["Mass Error (ppm)"] = new_products["Mass Error (ppm)"].apply(lambda x: f"{x:.2f}")
        new_products.to_excel(str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "代谢产物数据矩阵-所有候选.xlsx"), index=False)

    def oebio_report_deprecated(self, version="/data/hstore1/database/report/2023-03-04~2023/src", company="LM",
                     using_product=False, using_herb=False, cloud=True):

        from lmbio.basic.reportinfo import getreportinfo
        os.environ['OEBIO'] = version
        reportinfo = getreportinfo(sys_logo=company)
        if cloud:
            os.environ['CLOUD'] = "True"
            os.environ['CLOUDUSER'] = "lmbioinformatics@lumingbio.com"
            os.environ['PASSWORD'] = "Lmbio@123"
        logger.info("正在生成中药入血成分报告...")

        src_path = Path("/data/hstore1/database/database/tcm/")

        project_name, project_info = self.parse_project()
        # -------------------do before report-------------------
        # 批量重命名
        data_path = str(self.midware / "数据矩阵.xlsx")
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

        # -------------------------------------------------------
        blood_cols, control_cols, tcm_cols = self.parse_exp_col()
        ren_dict = {}
        new_tcm_cols = []
        for i, exp in enumerate(blood_cols, start=1):
            ren_dict[exp] = f"Treatment-{i}"
        for i, ctl in enumerate(control_cols, start=1):
            ren_dict[ctl] = f"Control-{i}"
        if len(tcm_cols) == 1:
            ren_dict[tcm_cols[0]] = "TCM"
            new_tcm_cols.append("TCM")
        else:
            for i, tcm in enumerate(tcm_cols, start=1):
                ren_dict[tcm] = f"TCM-{i}"
                new_tcm_cols.append(f"TCM-{i}")
        # 切分表格
        data_path = str(self.midware / "数据矩阵.xlsx")
        ingredient_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "中药成分鉴定数据矩阵.xlsx")
        blood_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "入血成分分析数据矩阵.xlsx")

        peaks = pd.read_excel(data_path, sheet_name=0).reset_index(drop=True)

        # rename
        peaks.rename(columns=ren_dict, inplace=True)
        #
        # -----------------------export products df----------------------------
        if using_product:
            product_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "代谢产物数据矩阵.xlsx")
            raw_product_path = str(self.midware / "raw代谢产物数据矩阵.xlsx")
            products = pd.read_excel(product_path, sheet_name=0).reset_index(drop=True)
            if 'Fragment Ions' not in products.columns:
                # add frag to products
                new_products = pd.merge(products, peaks[["Compound ID", "Metabolites", "Fragment Ions"]], how="inner",
                                        on="Compound ID")
                # drop_cols = ['Compound ID', '是否为入血成分', 'HMDB', 'METLIN', 'LipidMaps', 'KEGG', 'PubChem', 'CAS',
                # '中文名', '中文大类',
                # '中文子类', '英文大类', '英文子类', '厂家', '货号', '纯度', 'DB', 'TCM均值（峰面积的比值%）'] + blood_cols + control_cols + tcm_cols
                # TODO:临时导出， 用于从MYSQL查找候选
                new_products.to_excel(raw_product_path, index=False)

                new_products.rename(columns={"Compound": "ID"}, inplace=True)
                col_order = ["No.", 'ID', 'Mid', 'Parent Compound', 'Transformations', 'Adducts', 'Formula', 'Score', 'Fragmentation Score',
                             'Mass Error (ppm)', 'm/z',
                             'Retention time (min)', "Metabolites", 'Ion mode', 'Level', 'SMILES', 'Fragment Ions',
                             'Metabolism Type']
                new_products = new_products[col_order]
                # float handle
                new_products["m/z"] = new_products["m/z"].apply(lambda x: f"{x:.4f}")
                new_products["Retention time (min)"] = new_products["Retention time (min)"].apply(lambda x: f"{x:.2f}")
                new_products["Mass Error (ppm)"] = new_products["Mass Error (ppm)"].apply(lambda x: f"{x:.2f}")
                new_products.to_excel(product_path, index=False)
        # -----------------------export ingredients df----------------------------
        used_cols = ["No.", "ID", "Adducts", "Formula", "Score", "Fragmentation Score", "Mass Error (ppm)", "Metabolites",
                     "Fragment Ions", "theoretical m/z", "m/z", "Retention time (min)", "Ion mode",
                     *new_tcm_cols, 'TCM均值（峰面积的比值%）', "InChIKey", "SMILES", "HMDB", "METLIN", "LipidMaps", "KEGG", "PubChem", "CAS", "中文名",
                     "中文大类", "中文子类", "英文大类", "英文子类", "货号", "纯度", 'Level']
        # used_cols = ["ID", "Adducts", "Formula", "Score", "Fragmentation Score", "Mass Error (ppm)", "Metabolites",
        #              "Fragment Ions", "Retention time (min)", "Ion mode",
        #              *new_tcm_cols, 'TCM均值（峰面积的比值%）', "HMDB", "METLIN", "LipidMaps", "KEGG", "PubChem", "CAS", "中文名",
        #              "中文大类", "中文子类", "英文大类", "英文子类", "货号", "纯度", 'Level'] if ['theoretical m/z', 'InChIKey', 'SMILES'] not in used_cols else used_cols
        tcm_ingredients = peaks.query("DB != 'PRODUCT'")[used_cols]
        tcm_ingredients = tcm_ingredients.loc[tcm_ingredients[new_tcm_cols].mean(axis=1) >= 1000, :]
        tcm_ingredients["TCM均值（峰面积的比值%）"] = 100 * tcm_ingredients["TCM均值（峰面积的比值%）"] / tcm_ingredients[
            "TCM均值（峰面积的比值%）"].sum()
        tcm_ingredients.sort_values(by=["Level"]).to_excel(ingredient_path, index=False)

        # -----------------------export blood df----------------------------
        peaks.drop(columns=["Compound ID", "DB"]).sort_values(by=["是否为入血成分", "Level"], ascending=[False, True]).to_excel(blood_path, index=False)

        # split
        # os.unlink(data_path)
        # -------------------------------------------------------
        # -------------------do before report-------------------

        # 加载配置文件
        config = yaml.load(open(str(src_path / "blood" / "config.yaml"), encoding="utf8"), Loader=yaml.FullLoader)
        # 拷贝文件到实验内容
        report_xls = [path for path in (src_path / "tables").rglob("*.xlsx")] + \
                     [src_path / "blood" / "数据矩阵字段说明.xlsx"]
        if using_product:
            report_xls = report_xls + [src_path / "blood" / "中药代谢反应中英文对照及缩写详细说明表.xlsx"] + \
                         [src_path / "blood" / "中药代谢产物数据矩阵字段说明.xlsx"]
        for xls in report_xls:
            shutil.copy(str(xls), f"{str(Path(self.data) / project_name / '实验内容')}/{xls.name}")

        # 拷贝UV data
        uv_path = [path for path in (Path(self.data) / "UVdata").rglob("*紫外吸收图*")]
        for uv in uv_path:
            shutil.copy(str(uv), f"{str(Path(self.data) / project_name / '数据分析' / '1.紫外吸收图')}/{uv.name}")

        # 设置项目信息
        config["header_info"]["项目名称"] = '中药入血成分分析'
        # config["header_info"]["客户单位"] = project_info.loc[0, "客户单位"]
        config["header_info"]["任务单号"] = project_info.loc[0, "项目编号"] if "-b" in project_info.loc[0, "项目编号"] else project_info.loc[0, "项目编号"] + "-b1"
        config["header_info"]["客户名称"] = project_info.loc[0, "客户名称"]
        config["header_info"]["联系人名称"] = project_info.loc[0, "联系人名称"]
        config["header_info"]["项目编号"] = project_info.loc[0, "项目编号"]
        config["header_info"]["样本"] = project_info.loc[0, "样本"]
        # config["header_info"]["完成时间"] = datetime.datetime.now().strftime("%Y-%m-%d")
        # 默认 jiang tao
        # config["header_info"]["执行编码"] = 'LM0460'

        header_info = config.get("header_info")
        # OUT PATH
        report_dir = str(Path(self.data) / project_name)

        # 实例化报告对象
        report = oebioReport('中药入血成分分析', title='中药入血成分分析', header_info=header_info, oe_welcome = reportinfo.oe_weclome)

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
        # 紫外吸收图
        a = analysis.add_section('紫外吸收图', description=config["数据分析"]["紫外吸收图"]["description"])
        a.add_fig(f"{report_dir}/{config['数据分析']['紫外吸收图']['imgpath1']}", caption=config['数据分析']['紫外吸收图']['imgname1'])
        a.add_fig(f"{report_dir}/{config['数据分析']['紫外吸收图']['imgpath2']}", caption=config['数据分析']['紫外吸收图']['imgname2'])
        # BPC
        a1 = analysis.add_section('基峰图', description=config["数据分析"]['基峰图']['description'])
        a1.add_fig(f"{report_dir}/{config['数据分析']['基峰图']['imgpath1']}", caption=config['数据分析']['基峰图']['imgname1'])
        a1.add_fig(f"{report_dir}/{config['数据分析']['基峰图']['imgpath2']}", caption=config['数据分析']['基峰图']['imgname2'])
        # 数据预处理
        a2 = analysis.add_section('数据预处理')
        a2.add_section("Progenesis QI v3.0的定性分析",
                       description=config["数据分析"]['数据预处理']['Progenesis QI v3.0的定性分析']['description'])
        # 定性定量结果
        son3 = a2.add_section("定性定量结果", description=config["数据分析"]['数据预处理']['定性定量结果']['description'])

        son3.add_table(f"{report_dir}/{config['数据分析']['数据预处理']['定性定量结果']['tablepath1']}",
                       caption=config['数据分析']['数据预处理']['定性定量结果']['tablename1'])
        # 中药成分鉴定
        son3.add_table(f"{report_dir}/{config['数据分析']['数据预处理']['定性定量结果']['tablepath2']}",
                       caption=config['数据分析']['数据预处理']['定性定量结果']['tablename2'])

        son3.add_fig(f"{report_dir}/{config['数据分析']['数据预处理']['中药成分分类']['imgpath0']}",
                     caption=config['数据分析']['数据预处理']['中药成分分类']['imgname0'])

        # 中药入血
        son3.add_table(f"{report_dir}/{config['数据分析']['数据预处理']['定性定量结果']['tablepath3']}",
                       caption=config['数据分析']['数据预处理']['定性定量结果']['tablename3'], show_search=True)

        # 代谢产物
        if using_product:
            son3.add_table(f"{report_dir}/{config['数据分析']['数据预处理']['中药代谢产物网络图']['tablepath1']}",
                           caption=config['数据分析']['数据预处理']['中药代谢产物网络图']['tablename1'])
            son3.add_table(f"{report_dir}/{config['数据分析']['数据预处理']['定性定量结果']['tablepath4']}",
                           caption=config['数据分析']['数据预处理']['定性定量结果']['tablename4'], show_search=True)

        # MS/MS
        # 中药原方成分
        son4 = a2.add_section("中药原方成分的EIC图及其与标准品比对情况")
        # LuMet-TCM
        lumet_tcm = son4.add_section("中药原方成分的EIC图及其与标准品库LuMet-TCM的比对情况")
        lumet_tcm.add_fig(f"{config['数据分析']['数据预处理']['中药原方成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与标准品库LuMet-TCM的比对情况']['imgpath']}",
                     caption=config['数据分析']['数据预处理']['中药原方成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与标准品库LuMet-TCM的比对情况']['imgname'],
                     path=report_dir)
        if using_herb:
            # HerbDB
            herb_db = son4.add_section("中药原方成分的EIC图及其与公共库HerbDB的比对情况")
            herb_db.add_fig(f"{config['数据分析']['数据预处理']['中药原方成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与公共库HerbDB的比对情况']['imgpath']}",
                         caption=config['数据分析']['数据预处理']['中药原方成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与公共库HerbDB的比对情况']['imgname'],
                         path=report_dir)
        # 中药入血成分
        son41 = a2.add_section("中药入血成分的EIC图及其与标准品比对情况",
                               description=config['数据分析']['数据预处理']['中药入血成分的EIC图及其与标准品比对情况']["description"])
        son411 = son41.add_section("中药入血原型成分")
        # 中药原型
        son411.add_fig(f"{config['数据分析']['数据预处理']['中药入血成分的EIC图及其与标准品比对情况']['中药入血原型成分的EIC图及其与参考化合物比对情况']['imgpath']}",
                      caption=config['数据分析']['数据预处理']['中药入血成分的EIC图及其与标准品比对情况']['中药入血原型成分的EIC图及其与参考化合物比对情况']['imgname'],
                      path=report_dir)
        # 中药代谢产物
        if using_product:
            son412 = son41.add_section("中药代谢产物")
            # product
            son412.add_fig(f"{config['数据分析']['数据预处理']['中药入血成分的EIC图及其与标准品比对情况']['中药代谢产物的EIC图及其与参考化合物比对情况']['imgpath']}",
                         caption=config['数据分析']['数据预处理']['中药入血成分的EIC图及其与标准品比对情况']['中药代谢产物的EIC图及其与参考化合物比对情况']['imgname'],
                         path=report_dir)
        # pie
        son6 = a2.add_section("中药成分分类", description=config['数据分析']['数据预处理']['中药成分分类']['description'])
        son6.add_fig(f"{report_dir}/{config['数据分析']['数据预处理']['中药成分分类']['imgpath1']}",
                     caption=config['数据分析']['数据预处理']['中药成分分类']['imgname1'], description="为了显示美观，自动隐藏了占比<1%的数据标签")
        son6.add_fig(f"{report_dir}/{config['数据分析']['数据预处理']['中药成分分类']['imgpath2']}",
                     caption=config['数据分析']['数据预处理']['中药成分分类']['imgname2'], description="为了显示美观，自动隐藏了占比<1%的数据标签")

        # 代谢物网络图
        if using_product:
            son5 = a2.add_section("中药代谢产物网络图", description=config['数据分析']['数据预处理']['中药代谢产物网络图']["description"])
            son5.add_table(f"{report_dir}/{config['数据分析']['数据预处理']['中药代谢产物网络图']['tablepath']}",
                           caption=config['数据分析']['数据预处理']['中药代谢产物网络图']['tablename'])
            son5.add_fig(f"{config['数据分析']['数据预处理']['中药代谢产物网络图']['imgpath']}",
                         caption=config['数据分析']['数据预处理']['中药代谢产物网络图']['imgname'],
                         path=report_dir)

        # 公司简介
        reportinfo.companyinfo(report=report, yamlpath=os.environ['OEBIO'])

        # export
        report.write_to(f"{report_dir}/项目报告.html", zip_report_name=report_dir+"-OECloud.zip")

        # 添加与报告到 /public/项目报告检查空间/预报告文件夹/
        pnum = Path(self.data).name # project_info.loc[0, "项目编号"]
        pname = project_info.loc[0, "项目名称"]

        # shutil.copy(f"{report_dir}/项目报告.html", f"/public/项目报告检查空间/预报告文件夹/{pnum}-{pname}-LCMS-1.0.0-分析预报告.html")
        shutil.copy(f"{report_dir}/项目报告.html", f"/data/nas/177/代谢/项目报告检查空间/预报告文件夹/{pnum}-{pname}-LCMS-1.0.0-分析预报告.html")

    def oebio_report(self, version="/data/hstore1/database/report/2023-03-04~2023/src", company="LM",
                     using_product=False, using_herb=False, cloud=True, skip=False):
        import subprocess as sp
        script_path = str(Path(__file__).parent / "oecloud.py")
        local = "0" if self.local else "1"
        _skip = "1" if skip else "0"
        args = ["python", script_path, "-d", self.data, "-t", "blood", "-l", local, "-k", _skip]
        if using_herb:
            args.append("--herb")
        if using_product:
            args.append("--product")
        if cloud:
            args.append("--cloud")
        r = sp.run(args, stderr=sp.PIPE, stdout=sp.PIPE, text=True, encoding='utf8')
        if r.returncode != 0:
            logger.info(r.stderr)
            logger.info("oebio_report Blood Failed!")
        else:
            logger.info("oebio_report Blood SUCCESS!")

    def zip_report(self):
        project_name, info = self.parse_project()
        zip_in = self.data + os.sep + project_name
        zip_out = zip_in + "_local.zip"
        logger.info("正在生成压缩文件%s...", zip_out)

        zip_file = zipfile.ZipFile(zip_out, 'w', zipfile.ZIP_DEFLATED)
        for dir_path, dir_names, file_names in os.walk(zip_in):
            save_path = dir_path.replace(zip_in, '')
            # nas virus dir
            if Path(dir_path).name != "@eaDir" and Path(dir_path).name != "delete":
                for filename in file_names:
                    # nas virus file
                    # add "数据矩阵.xlsx", "raw代谢产物数据矩阵.xlsx" for this proj
                    if filename not in ["Thumbs.db", "数据矩阵.xlsx", "raw代谢产物数据矩阵.xlsx", "delCpds.xlsx", "delCpds1.xlsx", "preprocessed.xlsx"]\
                            and not filename.startswith('~$') and not filename.endswith(".tmp"):
                        zip_file.write(os.path.join(dir_path, filename), os.path.join(save_path, filename))
        zip_file.close()
        logger.info("压缩文件完成%s...", zip_out)

    def run_pipeline(self, using_herb=False, using_product=False, FC=10):
        logger.info("正在执行中药入血成分分析的完整pipeline...")
        self.raw2mzml()
        self.raw2obs()
        self.preprocess_qi(using_herb=using_herb, using_product=using_product, FC=FC)
        self.plot_chromatograms()
        if using_product:
            self.network_analysis()
        self.plot_xic()
        self.rename_excel()
        self.add_src()

    def manual_inspect(self, skip=False, using_herb=False, using_product=False, del_blood=True, cloud=True):
        if skip:
            project_name, project_info = self.parse_project()
            del_path = str(self.midware / "delCpds.xlsx")
            del_path_yf = str(self.midware / "delCpds1.xlsx")
            # tmp file
            tmp_path = str(self.midware / "数据矩阵.xlsx")

            data = pd.read_excel(tmp_path, sheet_name=0).set_index(keys=["No.", ], drop=False)
            # namesmap = dict(zip(data["Metabolites"], data.index))
            # ------------
            pdf, products, conn, ha, rules = None, None, None, None, None
            if using_product:
                # 代谢产物数据矩阵
                pdf = pd.read_excel(str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "代谢产物数据矩阵.xlsx"), sheet_name=0)
                # products
                products = data[data.DB == "PRODUCT"].reset_index(drop=True)
                products["group"] = products.apply(
                    lambda x: x["Compound ID"].split("_", 1)[0] if x.DB == "PRODUCT" else x["Metabolites"], axis=1)
                products = products.set_index(["group", "Metabolites"]).sort_index()
                products["Mid"] = ""
                products["SMILES"] = ""
                products["Parent Compound"] = ""
                products["Transformations"] = ""
                products["Metabolism Type"] = ""
                # load parent info
                engine = create_engine(_engine)
                conn = engine.connect()
                ha = pd.read_sql("herb_annotation", conn).set_index(['compound_name']).rename(
                    columns={"select": "筛选", "smiles": "SMILES", "inchikey": "InChIKey", "inchi": "InChI", "rt": "RT"})

                # read herb_rules
                rules = pd.read_sql("herb_rules", conn).set_index(['reaction'])
            if using_product:
                # -------------------重新绘制MS2 加载mzml
                logger.info("正在读取原方mzML...")
                dpi = 600  # 200
                img_size = 250  # 300
                # 获取项目信息
                project, project_info = self.parse_project()
                # 1.MSExp container for ex xic
                MSExps = dict()
                # mzML path
                mode_used, sample_info = self.parse_meta()
                labels = self.group
                colormap = plt.get_cmap('Set1')(range(3))
                # new
                blood_cols, control_cols, tcm_cols = self.parse_exp_col()

                blood_fs = list(map(lambda x: x + ".mzML", blood_cols))
                control_fs = list(map(lambda x: x + ".mzML", control_cols))
                tcm_fs = list(map(lambda x: x + ".mzML", tcm_cols))

                cb_cols = blood_cols + control_cols
                mzml_files = control_fs + blood_fs + tcm_fs
                for mode in mode_used:
                    MSExps[mode] = []
                    files = []
                    for file in mzml_files:
                        abs_f = Path(self.data) / "mzmldata" / mode / file
                        files.append(abs_f)
                    # load mzml
                    for file in files:
                        exp = MSExperiment()
                        MzMLFile().load(str(file), exp)
                        MSExps[mode].append(exp)
                # 2.MS2 container for ref MS/MS
                msp_path = Path(self.data) / "QIdata" / "MSP"
                refMS2: dict = self.loadMS2FromMsp(msp_path, mode_used)
            #  ----------------------------------------
            if del_blood:
                delCpds = pd.read_excel(del_path, sheet_name=0)["Compound ID"].apply(lambda x: x.strip()).drop_duplicates()
                data.query("`No.` in @delCpds").to_excel(str(self.midware / "delCpds_info.xlsx"), index=False)
                logger.info("正在删除人工检查去除的入血代谢物...")
                """人工检查去除"""
                pic_path = str(Path(self.data) / project_name / "数据分析" / "3.质谱图" / "中药入血成分" / "中药原型" / "png")
                pdf_path = str(Path(self.data) / project_name / "数据分析" / "3.质谱图" / "中药入血成分" / "中药原型" / "pdf")
                pic_path1 = str(Path(self.data) / project_name / "数据分析" / "3.质谱图" / "中药入血成分" / "中药代谢产物" / "png")
                pdf_path1 = str(Path(self.data) / project_name / "数据分析" / "3.质谱图" / "中药入血成分" / "中药代谢产物" / "pdf")
                pic_net = str(Path(self.data) / project_name / "数据分析" / "6.代谢网络图" / "png")
                pdf_net = str(Path(self.data) / project_name / "数据分析" / "6.代谢网络图" / "pdf")
                # 解压缩
                # zip_path = self.data + os.sep + project_name + ".zip"
                # zip_file = zipfile.ZipFile(zip_path)
                # zip_file.extractall(str(Path(self.data) / project_name))
                # products = pd.read_excel(product_path, sheet_name=0)
                # network_infected = []
                for cpdid in delCpds:
                    # TODO: tmp 避免报错
                    if cpdid not in data.index.to_list():
                        continue
                    cpd = data.loc[cpdid, "Metabolites"]
                    db_id = data.loc[cpdid, "Compound ID"]
                    # tmp
                    db = data.loc[cpdid, "DB"]
                    logger.info(f"正在删除入血代谢物 - {cpd}...")
                    # 修改 是否为入血成分 为 否
                    data.loc[cpdid, "是否为入血成分"] = "否"

                    filename = cpdid
                    if db == "PRODUCT":
                        g = cpd.split("_M", 1)[0]
                        net_name = self.filename_handler(g, escape=True, full=True)
                        sub = products.loc[g]
                        if sub.shape[0] == 1:
                            # 删除代谢产物数据矩阵中的行
                            idx = pdf.loc[pdf["Compound ID"] == db_id, :].index.values[0]
                            pdf.drop(index=idx, inplace=True)
                            # 删除数据矩阵中的行
                            data.drop(index=cpdid, inplace=True)
                            try:
                                # 删除MS2
                                os.unlink(pic_path1 + '/' + filename + '.png')
                                os.unlink(pdf_path1 + '/' + filename + '.pdf')
                                # 删除网络图
                                os.unlink(pic_net + '/' + net_name + '.png')
                                os.unlink(pdf_net + '/' + net_name + '.pdf')
                            except:
                                pass
                        elif sub.shape[0] > 1:
                            # 删除代谢产物数据矩阵中的行
                            idx = pdf.loc[pdf["Compound ID"] == db_id, :].index.values[0]
                            pdf.drop(index=idx, inplace=True)
                            # 删除数据矩阵中的行
                            data.drop(index=cpdid, inplace=True)
                            # 删除网络图
                            try:
                                os.unlink(pic_net + '/' + net_name + '.png')
                                os.unlink(pdf_net + '/' + net_name + '.pdf')
                            except:
                                pass
                            # 删除所有MS
                            for i, row in sub.iterrows():
                                fn = row["No."]
                                try:
                                    os.unlink(pic_path1 + '/' + fn + '.png')
                                    os.unlink(pdf_path1 + '/' + fn + '.pdf')
                                except:
                                    pass
                            # 删除 in sub
                            index = sub[sub["Compound ID"] == db_id].index[0]
                            sub.drop(index=index, inplace=True)
                            # update products
                            ind = products[products["Compound ID"] == db_id].index
                            products.drop(index=ind, inplace=True)

                            logger.info("重新绘制代谢网络图...")
                            # prepare network input
                            parent_name = sub.index[0].split("_", 1)[0]
                            p = ha.loc[parent_name, ["SMILES", "InChIKey", "InChI", "RT"]]
                            p = p.iloc[0, :] if isinstance(p, pd.DataFrame) else p
                            parent = pd.DataFrame(data=p.to_dict(), index=[-1])
                            db_idx = sub.columns.get_loc('DB')
                            _idx = sub.columns.get_loc('Compound ID')
                            _ = tuple(sub.iloc[1:, _idx]) if sub.iloc[0, db_idx] != "PRODUCT" else tuple(sub.iloc[:, _idx])
                            if len(_) == 1:
                                _ = f"""("{_[0]}")"""
                            sql = text(f"""select * from herb_network WHERE network_id in {_}""")
                            prod = pd.read_sql(sql, conn).drop(
                                columns=["pid", "ion_mode", "parent_smiles", "parent_inchi", "parent_inchikey",
                                         "parent_name",
                                         "batch"]) \
                                .rename(
                                columns={"select": "筛选", "ReactionClass": "Reaction Class",
                                         "Intermediate1": "Intermediate 1",
                                         "Intermediate2": "Intermediate 2"})
                            data_net = pd.concat([prod, parent]).sort_index()
                            # import pdb
                            # pdb.set_trace()
                            # plot func
                            cpd_name = self.filename_handler(parent_name, full=True)
                            save_path = str(Path(self.data) / project_name / "数据分析" / "6.代谢网络图")
                            new_id: pd.DataFrame = plot_networkx(data_net, save_path, cpd_name)
                            ids = new_id.Mid.apply(lambda x: x.split('-')[0])
                            clear_cpd: pd.Series = ids[~ids.duplicated(keep=False)]
                            sub = sub.reset_index()
                            for i, row in new_id.iterrows():
                                if i != 0:
                                    net_id, mid, rc, smiles = row
                                    mid_fixed = mid.split("-")[0]
                                    mid = mid_fixed if mid_fixed in clear_cpd.values else mid
                                    pdf.loc[pdf["Compound ID"] == net_id, "Mid"] = mid
                                    pdf.loc[pdf["Compound ID"] == net_id, "Parent Compound"] = parent_name
                                    data.loc[data["Compound ID"] == net_id, "Metabolites"] = parent_name + "_" + mid
                                    pdf.loc[pdf["Compound ID"] == net_id, "SMILES"] = smiles
                                    sub.loc[sub["Compound ID"] == net_id, "Metabolites"] = parent_name + "_" + mid
                                    # 处理反应类型 并添加 Transformations & Metabolism Type
                                    transformations = []
                                    _type = []
                                    for r in rc.split("+", -1):
                                        _ = r.strip()
                                        en, phase = rules.loc[_, ["en", "phase"]]

                                        transformations.append(en)
                                        _type.append(phase)
                                    pdf.loc[pdf["Compound ID"] == net_id, "Transformations"] = ",".join(transformations)
                                    pdf.loc[pdf["Compound ID"] == net_id, "Metabolism Type"] = ",".join(sorted(set(_type)))
                            # 重新绘制MS2
                            logger.info("重新绘制MS2...")
                            for i, row in sub.iterrows():
                                cpdn, _id, inchikey, adduct, formula, name, _mz, _rt, ion_mode, is_blood, db = \
                                    row["No."], row["ID"], row["Compound ID"], row["Adducts"], row["Formula"], row["Metabolites"], row["m/z"], row["Retention time (min)"], row["Ion mode"], row["是否为入血成分"], row["DB"]
                                if cpdn is None:
                                    raise ValueError(f"不存在该化合物: {name}")
                                controls = np.asarray(row[control_cols].tolist())
                                bloods = np.asarray(row[blood_cols].tolist())
                                tcms = np.asarray(row[tcm_cols].tolist())

                                index_control = np.argmin(controls)
                                index_blood = np.argmax(bloods) + len(controls)
                                index_tcm = np.argmax(tcms) + len(cb_cols)

                                used_exp = [MSExps[ion_mode][index_control], MSExps[ion_mode][index_blood],
                                            MSExps[ion_mode][index_tcm]]
                                # 3.1 extract data from mzML & plot XIC
                                # pyplot top plot
                                is_blood = True if is_blood == "是" else False
                                self.plot_ms2(cpdn, project, _id, inchikey, adduct, formula, name,
                                              _mz, _rt, ion_mode, is_blood, used_exp, refMS2, 5, labels, colormap, img_size, dpi, 'png', db)
                        else:
                            raise ValueError(f"不存在该id: {cpdid}")
                    else:
                        try:
                            os.unlink(pic_path + '/' + filename + '.png')
                            os.unlink(pdf_path + '/' + filename + '.pdf')
                        except:
                            pass
                conn.close() if conn is not None else ''

            else:
                # 删除原方代谢物
                delCpds1 = pd.read_excel(del_path_yf, sheet_name=0)["Compound ID"].apply(lambda x: x.strip()).drop_duplicates()
                logger.info("正在删除人工检查去除的原方代谢物...")
                for cpdid in delCpds1:
                    cpd = data.loc[cpdid, "Metabolites"]
                    db = data.loc[cpdid, "DB"]
                    logger.info(f"正在删除原方代谢物 - {cpd}...")
                    folder = "LuMet-TCM" if db in ["TCM", "ANIMAL"] else "HerbDB" if db == "HERB" else ValueError('删除的化合物不支持该数据库来源')
                    pic_path = str(Path(self.data) / project_name / "数据分析" / "3.质谱图" / "中药原方成分" / folder / "png")
                    pdf_path = str(Path(self.data) / project_name / "数据分析" / "3.质谱图" / "中药原方成分" / folder / "pdf")
                    data.drop(index=[cpdid, ], inplace=True)
                    filename = cpdid
                    os.unlink(pic_path + '/' + filename + '.png')
                    os.unlink(pdf_path + '/' + filename + '.pdf')

            # force to 100
            data["TCM均值（峰面积的比值%）"] = 100 * data["TCM均值（峰面积的比值%）"] / data["TCM均值（峰面积的比值%）"].sum()
            # sort
            data.sort_values(by=["是否为入血成分", "Level"], ascending=[False, True], inplace=True)
            # export
            data.reset_index(drop=True).to_excel(tmp_path, index=False)
            if using_product:
                pdf.to_excel(str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "代谢产物数据矩阵.xlsx"), index=False)
            # --------------
            # os.unlink(del_path)
        else:
            self.relation_analysis()
            self.plot_pies(water_mask=False)
            self.oebio_report(using_herb=using_herb, using_product=using_product, cloud=cloud)
            self.zip_report()

    def split_df_only(self, using_product=False):
        project_name, project_info = self.parse_project()
        # -------------------------------------------------------
        blood_cols, control_cols, tcm_cols = self.parse_exp_col()
        ren_dict = {}
        new_tcm_cols = []
        for i, exp in enumerate(blood_cols, start=1):
            ren_dict[exp] = f"Treatment-{i}"
        for i, ctl in enumerate(control_cols, start=1):
            ren_dict[ctl] = f"Control-{i}"
        if len(tcm_cols) == 1:
            ren_dict[tcm_cols[0]] = "TCM"
            new_tcm_cols.append("TCM")
        else:
            for i, tcm in enumerate(tcm_cols, start=1):
                ren_dict[tcm] = f"TCM-{i}"
                new_tcm_cols.append(f"TCM-{i}")
        # 切分表格
        data_path = str(self.midware / "数据矩阵.xlsx")
        ingredient_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "中药成分鉴定数据矩阵.xlsx")
        blood_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "入血成分分析数据矩阵.xlsx")

        peaks = pd.read_excel(data_path, sheet_name=0).reset_index(drop=True)

        # rename
        peaks.rename(columns=ren_dict, inplace=True)

        # -----------------------export products df----------------------------
        if using_product:
            product_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "代谢产物数据矩阵.xlsx")
            raw_product_path = str(self.midware / "raw代谢产物数据矩阵.xlsx")
            products = pd.read_excel(product_path, sheet_name=0).reset_index(drop=True)
            # add frag to products
            new_products = pd.merge(products, peaks[["Compound ID", "Metabolites", "Fragment Ions"]], how="inner",
                                    on="Compound ID")
            # drop_cols = ['Compound ID', '是否为入血成分', 'HMDB', 'METLIN', 'LipidMaps', 'KEGG', 'PubChem', 'CAS',
            # '中文名', '中文大类',
            # '中文子类', '英文大类', '英文子类', '厂家', '货号', '纯度', 'DB', 'TCM均值（峰面积的比值%）'] + blood_cols + control_cols + tcm_cols
            # TODO:临时导出， 用于从MYSQL查找候选
            new_products.to_excel(raw_product_path, index=False)

            new_products.rename(columns={"Compound": "ID"}, inplace=True)
            col_order = ['ID', 'Mid', 'Parent Compound', 'Transformations', 'Adducts', 'Formula', 'Score', 'Fragmentation Score',
                         'Mass Error (ppm)', 'm/z',
                         'Retention time (min)', "Metabolites", 'Ion mode', 'Level', 'SMILES', 'Fragment Ions',
                         'Metabolism Type']
            new_products = new_products[col_order]
            # float handle
            new_products["m/z"] = new_products["m/z"].apply(lambda x: f"{x:.4f}")
            new_products["Retention time (min)"] = new_products["Retention time (min)"].apply(lambda x: f"{x:.2f}")
            new_products["Mass Error (ppm)"] = new_products["Mass Error (ppm)"].apply(lambda x: f"{x:.2f}")
            new_products.to_excel(product_path, index=False)

        # -----------------------export ingredients df----------------------------
        used_cols = ["ID", "Adducts", "Formula", "Score", "Fragmentation Score", "Mass Error (ppm)", "Metabolites",
                     "Fragment Ions", "m/z", "Retention time (min)", "Ion mode",
                     *new_tcm_cols, 'TCM均值（峰面积的比值%）', "HMDB", "METLIN", "LipidMaps", "KEGG", "PubChem", "CAS", "中文名",
                     "中文大类", "中文子类", "英文大类", "英文子类", "货号", "纯度", 'Level']
        used_cols = used_cols if '来源' in peaks.columns else used_cols+['来源']
        tcm_ingredients = peaks.query("DB != 'PRODUCT'")[used_cols]
        tcm_ingredients = tcm_ingredients.loc[tcm_ingredients[new_tcm_cols].mean(axis=1) >= 1000, :]
        tcm_ingredients["TCM均值（峰面积的比值%）"] = 100 * tcm_ingredients["TCM均值（峰面积的比值%）"] / tcm_ingredients[
            "TCM均值（峰面积的比值%）"].sum()
        tcm_ingredients.sort_values(by=["Level"]).to_excel(ingredient_path, index=False)

        # -----------------------export blood df----------------------------
        peaks.drop(columns=["Compound ID", "DB"]).sort_values(by=["是否为入血成分", "Level"],
                                                              ascending=[False, True]).to_excel(blood_path, index=False)





if __name__ == '__main__':
    ...
