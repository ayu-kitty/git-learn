# File: ProductResolve
# Author: jiang tao
# Time: 2023/7/11 9:42
# Desc: 解析用户所提交的代谢产物的SOP
import logging
import os
import shutil
from pathlib import Path
from typing import List, Union

import pandas as pd
import subprocess as sp

import yaml
from pyopenms import MSSpectrum
from sqlalchemy import text, Engine
from tqdm import tqdm

from .Blood2report import BloodReport
from .exception import DatabaseNotImplementError, SqlExecuteError
from .network import plot_networkx
from oebio.report import Report as oebioReport

logger = logging.getLogger("TCMreport")


class ProductResolve(BloodReport):

    sql1 = """SELECT
                    InChIKey,
                    SMILES
                FROM
                    herb_network_local
                WHERE
                    network_id = "%s" and batch = "%s"
            """

    sql2 = """SELECT
                    InChIKey,
                    SMILES
                FROM
                    herb_network_local
                WHERE
                    InChIKey = "%s" and batch = "%s"
            """

    def is_cpd_valid(self, inchikey, is_product=False):
        with self.ENGINE.connect() as conn:
            # 从标准品库查
            sql = f"""select * from herb_network_local where network_id={repr(inchikey)}""" if is_product else f"""select * from herb_annotation where inchikey={repr(inchikey)}"""
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
            with open(str(src_path / "metfrag" / "localdb.csv"), "w", encoding="utf-8") as d:
                d.write(
                    f""""Identifier","MonoisotopicMass","MolecularFormula","SMILES","InChI","InChIKey1","InChIKey2","InChIKey3","Name","InChIKey"\n""")
                # 从标准品库查
                # sql = f"""select * from herb_product where `NAME`={repr(inchikey)}""" if is_product else f"""select * from herb_msms where inchikey={repr(inchikey)}"""
                sql = f"""SELECT MonoisotopicMass AS exact_mass, SMILES AS smiles, InChIKey AS inchikey, InChI AS inchi FROM herb_network_local where network_id={repr(inchikey)}""" if is_product else f"""select * from herb_annotation where inchikey={repr(inchikey)}"""
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
                fetch_data = res if res.shape[0] else _res # if _res.shape[0] else _res1
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
            with open(str(src_path / "metfrag" / "MS2.txt"), "w", encoding="utf-8") as f:
                for mz, ints in zip(_mz, _intensity):
                    f.write(f"{mz}\t{ints}\n")
            with open(str(src_path / "metfrag" / "parameter_template.txt"), "r", encoding="utf-8") as p:
                content = p.read().replace("NeutralPrecursorMolecularFormula = C10H10O3",
                                           f"NeutralPrecursorMolecularFormula = {formula}"). \
                    replace("PeakListPath = MS2.txt", f"PeakListPath = {str(src_path / 'metfrag' / 'MS2.txt')}"). \
                    replace("LocalDatabasePath = localdb.csv",
                            f"LocalDatabasePath = {str(src_path / 'metfrag' / 'localdb.csv')}") \
                    .replace("ResultsPath = .", f"ResultsPath = {str(src_path / 'metfrag')}")
                if ion_mode == "NEG":
                    content = content.replace("IsPositiveIonMode = True", f"IsPositiveIonMode = False").replace(
                        "PrecursorIonMode = 1", f"PrecursorIonMode = {index}")
                else:
                    content = content.replace("PrecursorIonMode = 1", f"PrecursorIonMode = {index}")
                with open(str(src_path / "metfrag" / "parameter_file.txt"), "w", encoding="utf-8") as ff:
                    ff.write(content)

            # run MetFrag in cmd using java >= 8
            # status = os.system(
            #     f"java -jar /root/MetFrag2.4.5-CL.jar {str(src_path / 'metfrag' / 'parameter_file.txt')}")
            # use subprocess instead
            r = sp.run(['java', '-jar', '/root/MetFrag2.4.5-CL.jar', f"{str(src_path / 'metfrag' / 'parameter_file.txt')}"],
                            stdout=sp.PIPE, stderr=sp.STDOUT)
            if r.returncode == 0:
                with open(str(src_path / "metfrag" / "FragmentSmilesExpl.psv"), "r", encoding="utf-8") as file:
                    content = file.read().strip()
                if content:
                    df = pd.read_csv(str(src_path / "metfrag" / "FragmentSmilesExpl.psv"), sep="|", header=0,
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
        conn.close()
        return None

    def get_parent_smiles(cls, inchikey: str, is_product=False) -> Union[None, str]:
        engine = cls.ENGINE
        conn = engine.connect()
        # 从标准品库查
        sql = f"""select * from herb_network_local where network_id={repr(inchikey)}""" if is_product else f"""select * from herb_annotation where inchikey={repr(inchikey)}"""
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

    def fill_qi_df(self, data_processed: List[pd.DataFrame]) -> pd.DataFrame:
        """post_process steps in preprocess_qi"""
        # concat res in two ion modes
        _ , proj_info = self.parse_project()
        batch = proj_info.loc[0, "项目编号"]
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
                sql = self.sql1 % (cid, batch)
                res_proxy = conn.execute(text(sql)).fetchall()
                if len(res_proxy) != 0:
                    inchikey, smiles = res_proxy[0][0], res_proxy[0][1]
                    res.loc[i, "InChIKey"] = inchikey
                    res.loc[i, "SMILES"] = smiles
                    continue
            elif db == "HERB":
                sql = self.sql2 % (cid, batch)
                res_proxy = conn.execute(text(sql)).fetchall()
                if len(res_proxy) != 0:
                    inchikey, smiles = res_proxy[0][0], res_proxy[0][1]
                    res.loc[i, "InChIKey"] = inchikey
                    res.loc[i, "SMILES"] = smiles
                    continue
            else:
                raise DatabaseNotImplementError("{} : 代谢产物解析不支持此数据库类型".format(db))
        conn.close()

        # del duplicated cpd and select best
        res = res.sort_values("Compound ID").reset_index(drop=True)
        best_index = self.selectIndex(res, "Compound ID")
        res = res.loc[best_index, :].sort_values("Compound").reset_index(drop=True)
        return res

    def preprocess_qi(self, FC=10, treshold_ctl=5000, min_treshold_exp=1000, treshold_tcm=2500, score_exp=48,
                      peak_width: tuple = (0.075, 0.40),
                      export_unknown=False, using_herb=False, using_product=False, interested_cpd=None):

        MIDWARE_PATH = Path(self.data) / "midware"
        MIDWARE_PATH.mkdir(exist_ok=True)
        logger.info("正在预处理QI结果，添加代谢物注释...")
        logger.info(f"差异倍数： {FC}")
        mode_used, sample_info = self.parse_meta()

        blood_cols, control_cols, tcm_cols = self.parse_exp_col()
        # calc n
        n = len(blood_cols + control_cols + tcm_cols) + 1

        data_processed_blood = []

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

            # 只分析代谢产物[PRODUCT] + 母体[HERB]
            qual_quant_blood = qual_quant_blood[(qual_quant_blood["DB"] == "PRODUCT") | (
                        "HERB" == qual_quant_blood["DB"])]

            # add annotation level
            qual_quant_blood["Level"] = qual_quant_blood.apply(
                lambda x: "level3" if x["DB"] == "HERB" else "level1" if x["Fragmentation Score"] >= 50 else "level2", axis=1)

            if interested_cpd:
                pass

            # 不卡得分

            best_index1 = self.selectIndex(qual_quant_blood, "Compound")
            qual_quant_blood = qual_quant_blood.loc[best_index1, :]
            qual_quant_blood["Ion mode"] = mode
            # add M mean
            qual_quant_blood["TCM均值"] = qual_quant_blood[tcm_cols].mean(axis=1)
            qual_quant_blood["TCM均值（峰面积的比值%）"] = 100 * qual_quant_blood["TCM均值"] / qual_quant_blood["TCM均值"].sum()

            data_processed_blood.append(qual_quant_blood)

        res_blood = self.fill_qi_df(data_processed_blood)
        res_blood["theoretical m/z"] = (1e6 * res_blood["m/z"]) / (1e6 + res_blood["Mass Error (ppm)"])

        # 排除 0值
        # res_blood = res_blood[res_blood["TCM均值"] > treshold_tcm]
        # 扣除空白 FC>=4
        ctl_blood = res_blood[control_cols].mean(axis=1)
        exp_blood = res_blood[blood_cols].mean(axis=1)
        res_blood["是否为入血成分"] = "否"

        # 原型入血成分
        res_blood.loc[(res_blood["DB"] != "PRODUCT") &
                      # (res_blood["TCM均值"] >= treshold_tcm) &
                      (ctl_blood <= treshold_ctl) &
                      (exp_blood >= min_treshold_exp) &
                      (exp_blood / ctl_blood >= FC)
                      # & (res_blood["Score"] >= score_exp)
                        , "是否为入血成分"] = "是"
        # 入血代谢产物
        if using_product:
            res_blood.loc[(res_blood["DB"] == "PRODUCT") &
                          # (res_blood["TCM均值"] < treshold_tcm) &
                          (ctl_blood <= treshold_ctl) &
                          (exp_blood >= min_treshold_exp) &
                          (exp_blood / ctl_blood >= FC)
                          # & (res_blood["Score"] >= score_exp)
                          ,"是否为入血成分"] = "是"
            # 过滤掉 不符合要求的代谢产物
            res_blood = res_blood.loc[(res_blood["DB"] != "PRODUCT") | ((res_blood["DB"] == "PRODUCT") & (res_blood["是否为入血成分"] == "是"))]
        res_blood_is = res_blood.query("`是否为入血成分` == '是'")

        # 中药原方
        res_blood_not = res_blood.query("`是否为入血成分` == '否'")
        res_blood_not = res_blood_not[#(res_blood_not["TCM均值"] > treshold_tcm) &
                                      (res_blood_not["DB"] != "PRODUCT")]

        res_blood = pd.concat([res_blood_is, res_blood_not], ignore_index=True)
        # 根据名称去重
        res_blood["lower_name"] = res_blood['Description'].apply(lambda x: x.lower())
        res_blood = res_blood.sort_values(by="lower_name").reset_index(drop=True)
        indexOfBestName = self.selectIndexOfBestName(res_blood, 'lower_name')
        res_blood = res_blood.loc[indexOfBestName, :].reset_index(drop=True).applymap(lambda x:x.strip() if isinstance(x, str) else x)
        # add compound No.
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
        # # 删除动物+lipid, 如果中药是完全植物源
        # is_del_animal = False if project_info.loc[0, "复方中是否包含动物"] == '是' else True
        # if is_del_animal:
        #     with self.ENGINE.connect() as tx:
        #         df = pd.read_sql('herb_annotation', tx)
        #         dels = df[(df.is_animal_source.astype(int) + df.is_lipid.astype(int)) >= 1]["inchikey"].to_list()
        #     res_blood = res_blood.query('`Compound ID` not in @dels').reset_index(drop=True)

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
        # (SAVE_PATH / "数据分析" / "3.质谱图" / "中药原方成分" / "LuMet-TCM" / "png").mkdir(exist_ok=True, parents=True)
        # (SAVE_PATH / "数据分析" / "3.质谱图" / "中药原方成分" / "LuMet-TCM" / "pdf").mkdir(exist_ok=True, parents=True)
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

        res_blood.to_excel(str(self.midware / "数据矩阵.xlsx"), sheet_name="data", encoding="utf-8",
                           index=False)
        res_blood.to_excel(str(self.midware / "preprocessed.xlsx"), sheet_name="data", encoding="utf-8",
                           index=False)

    def adjustRTByDetected(self):
        project_path, project_info = self.parse_project()
        batch = project_info.loc[0, "项目编号"]

        pp = str(self.midware / "数据矩阵.xlsx")
        peak = pd.read_excel(pp, sheet_name=0)

        sql = """UPDATE herb_network_local SET RT="%s", network_id="%s" WHERE SMILES="%s" and batch="%s" """

        sql1 = """SELECT InChIKey FROM herb_network_local WHERE SMILES="%s" and batch="%s" """
        conn = self.ENGINE.connect()
        for i, row in peak.iterrows():
            smi, rt, net_id = row["SMILES"], row["Retention time (min)"], row["Compound ID"]
            _id = net_id.rsplit('_', 1)[0] + '_' + str(rt)
            _sql = text(sql % (rt, _id, smi, batch))

            _sql1 = text(sql1 % (smi, batch))
            inchikey = pd.read_sql(_sql1, conn).iloc[0, 0]

            peak.loc[i, "Compound ID"] = _id
            peak.loc[i, "Description"] = _id

            peak.loc[i, "InChIKey"] = inchikey
            try:
                conn.execute(_sql)
                conn.commit()
            except Exception as e:
                print(e)
                conn.rollback()
                raise SqlExecuteError("{}: 保留时间更新失败".format(net_id))
        conn.close()
        peak.to_excel(pp, index=False)

    def network_analysis(self):
        # do this step after process qi
        logger.info("正在绘制代谢网络图...")

        project_name, project_info = self.parse_project()
        batch = project_info.loc[0, "项目编号"]
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
        conn = self.ENGINE.connect()
        sql = text("""SELECT * FROM herb_network_local WHERE batch="%s" """ % batch)
        ha = pd.read_sql(sql, conn).set_index(['StructureID'])

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
                    sql = text(f"""select * from herb_network_local WHERE network_id in {_} and batch={repr(batch)}""")
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
        #
        # # 批量重命名
        # data_path = str(self.midware / "数据矩阵.xlsx")
        # df = pd.read_excel(data_path, sheet_name=0).applymap(lambda x: x.strip() if isinstance(x, str) else x)
        # new_names = pd.Series(data=[f"compound{x:0>5d}" for x in range(1, df.shape[0] + 1)])
        # name_map = dict(zip(df["No."], new_names))
        # fig_path = Path(self.data) / project_name / "数据分析" / "3.质谱图"
        # formats = ("pdf", "png")
        # for fm in formats:
        #     for p in fig_path.rglob(f"*.{fm}"):
        #         _ = name_map.get(p.stem, None)
        #         if _:
        #             dst = Path(p.parent) / (_ + f".{fm}")
        #             os.rename(p, dst)
        #         else:
        #             logger.info(f"文件名{p.stem} 映射出错!")
        # df["No."] = new_names
        # df.to_excel(data_path, index=False)
        #
        # # -------------------------------------------------------
        # blood_cols, control_cols, tcm_cols = self.parse_exp_col()
        # ren_dict = {}
        # new_tcm_cols = []
        # for i, exp in enumerate(blood_cols, start=1):
        #     ren_dict[exp] = f"Treatment-{i}"
        # for i, ctl in enumerate(control_cols, start=1):
        #     ren_dict[ctl] = f"Control-{i}"
        # if len(tcm_cols) == 1:
        #     ren_dict[tcm_cols[0]] = "TCM"
        #     new_tcm_cols.append("TCM")
        # else:
        #     for i, tcm in enumerate(tcm_cols, start=1):
        #         ren_dict[tcm] = f"TCM-{i}"
        #         new_tcm_cols.append(f"TCM-{i}")
        # # 切分表格
        # data_path = str(self.midware / "数据矩阵.xlsx")
        # ingredient_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "中药成分鉴定数据矩阵.xlsx")
        # blood_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "入血成分分析数据矩阵.xlsx")
        #
        # peaks = pd.read_excel(data_path, sheet_name=0).reset_index(drop=True)
        #
        # # rename
        # peaks.rename(columns=ren_dict, inplace=True)
        # #
        # # -----------------------export products df----------------------------
        # if using_product:
        #     product_path = str(Path(self.data) / project_name / "数据分析" / "4.定性定量结果" / "代谢产物数据矩阵.xlsx")
        #     raw_product_path = str(self.midware / "raw代谢产物数据矩阵.xlsx")
        #     products = pd.read_excel(product_path, sheet_name=0).reset_index(drop=True)
        #     if 'Fragment Ions' not in products.columns:
        #         # add frag to products
        #         new_products = pd.merge(products, peaks[["Compound ID", "Metabolites", "Fragment Ions"]], how="inner",
        #                                 on="Compound ID")
        #         # drop_cols = ['Compound ID', '是否为入血成分', 'HMDB', 'METLIN', 'LipidMaps', 'KEGG', 'PubChem', 'CAS',
        #         # '中文名', '中文大类',
        #         # '中文子类', '英文大类', '英文子类', '厂家', '货号', '纯度', 'DB', 'TCM均值（峰面积的比值%）'] + blood_cols + control_cols + tcm_cols
        #         # TODO:临时导出， 用于从MYSQL查找候选
        #         new_products.to_excel(raw_product_path, index=False)
        #
        #         new_products.rename(columns={"Compound": "ID"}, inplace=True)
        #         col_order = ["No.", 'ID', 'Mid', 'Parent Compound', 'Transformations', 'Adducts', 'Formula', 'Score', 'Fragmentation Score',
        #                      'Mass Error (ppm)', 'm/z',
        #                      'Retention time (min)', "Metabolites", 'Ion mode', 'Level', 'SMILES', 'Fragment Ions',
        #                      'Metabolism Type']
        #         new_products = new_products[col_order]
        #         # float handle
        #         new_products["m/z"] = new_products["m/z"].apply(lambda x: f"{x:.4f}")
        #         new_products["Retention time (min)"] = new_products["Retention time (min)"].apply(lambda x: f"{x:.2f}")
        #         new_products["Mass Error (ppm)"] = new_products["Mass Error (ppm)"].apply(lambda x: f"{x:.2f}")
        #         new_products.to_excel(product_path, index=False)
        # # -----------------------export ingredients df----------------------------
        # used_cols = ["No.", "ID", "Adducts", "Formula", "Score", "Fragmentation Score", "Mass Error (ppm)", "Metabolites",
        #              "Fragment Ions", "theoretical m/z", "m/z", "Retention time (min)", "Ion mode",
        #              *new_tcm_cols, 'TCM均值（峰面积的比值%）', "InChIKey", "SMILES", "HMDB", "METLIN", "LipidMaps", "KEGG", "PubChem", "CAS", "中文名",
        #              "中文大类", "中文子类", "英文大类", "英文子类", "货号", "纯度", 'Level']
        # # used_cols = ["ID", "Adducts", "Formula", "Score", "Fragmentation Score", "Mass Error (ppm)", "Metabolites",
        # #              "Fragment Ions", "Retention time (min)", "Ion mode",
        # #              *new_tcm_cols, 'TCM均值（峰面积的比值%）', "HMDB", "METLIN", "LipidMaps", "KEGG", "PubChem", "CAS", "中文名",
        # #              "中文大类", "中文子类", "英文大类", "英文子类", "货号", "纯度", 'Level'] if ['theoretical m/z', 'InChIKey', 'SMILES'] not in used_cols else used_cols
        # tcm_ingredients = peaks.query("DB != 'PRODUCT'")[used_cols]
        # tcm_ingredients = tcm_ingredients.loc[tcm_ingredients[new_tcm_cols].mean(axis=1) >= 1000, :]
        # tcm_ingredients["TCM均值（峰面积的比值%）"] = 100 * tcm_ingredients["TCM均值（峰面积的比值%）"] / tcm_ingredients[
        #     "TCM均值（峰面积的比值%）"].sum()
        # tcm_ingredients.sort_values(by=["Level"]).to_excel(ingredient_path, index=False)
        #
        # # -----------------------export blood df----------------------------
        # peaks.drop(columns=["Compound ID", "DB"]).sort_values(by=["是否为入血成分", "Level"], ascending=[False, True]).to_excel(blood_path, index=False)

        # split
        # os.unlink(data_path)
        # -------------------------------------------------------
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

        # son3.add_fig(f"{report_dir}/{config['数据分析']['数据预处理']['中药成分分类']['imgpath0']}",
        #              caption=config['数据分析']['数据预处理']['中药成分分类']['imgname0'])

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
        # son6 = a2.add_section("中药成分分类", description=config['数据分析']['数据预处理']['中药成分分类']['description'])
        # son6.add_fig(f"{report_dir}/{config['数据分析']['数据预处理']['中药成分分类']['imgpath1']}",
        #              caption=config['数据分析']['数据预处理']['中药成分分类']['imgname1'], description="为了显示美观，自动隐藏了占比<1%的数据标签")
        # son6.add_fig(f"{report_dir}/{config['数据分析']['数据预处理']['中药成分分类']['imgpath2']}",
        #              caption=config['数据分析']['数据预处理']['中药成分分类']['imgname2'], description="为了显示美观，自动隐藏了占比<1%的数据标签")

        # 代谢物网络图
        if using_product:
            son5 = a2.add_section("中药代谢产物网络图", description=config['数据分析']['数据预处理']['中药代谢产物网络图']["description"])
            son5.add_table(f"{report_dir}/{config['数据分析']['数据预处理']['中药代谢产物网络图']['tablepath']}",
                           caption=config['数据分析']['数据预处理']['中药代谢产物网络图']['tablename'])
            son5.add_fig(
                "数据分析/6.代谢网络图/*",
                #f"{config['数据分析']['数据预处理']['中药代谢产物网络图']['imgpath']}",
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
        script_path = str(Path(__file__).parent / "oecloud.py")
        local = "0" if self.local else "1"
        args = ["python", script_path, "-d", self.data, "-t", "product", "-l", local]
        if using_herb:
            args.append("--herb")
        if using_product:
            args.append("--product")
        if cloud:
            args.append("--cloud")
        r = sp.run(args, stderr=sp.PIPE, stdout=sp.PIPE, text=True, encoding='utf8')
        if r.returncode != 0:
            logger.info(r.stderr)
            logger.info("oebio_report Product Failed!")
        else:
            logger.info("oebio_report Product SUCCESS!")

