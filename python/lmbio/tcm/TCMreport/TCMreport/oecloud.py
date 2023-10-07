# File: oecloud
# Author: jiang tao
# Time: 2023/7/25 10:59
import functools
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Tuple
import subprocess as sp
import getpass
import datetime

import requests
from oebio.report import Report as oebioReport
from glob import glob
import pandas as pd
import yaml

from Constants import BASE_URL, TOKEN
from exception import SaveOnlineRegisterTableError


class FilteredPrinter(object):
    def __init__(self, filtered_print, stdout, kw):
        self._filtered_print = filtered_print
        self._stderr = stdout
        self._kw = kw

    def _write(self, string):
        self._filtered_print(string, self._stderr, self._kw)

    def __getattr__(self, attr):
        if attr == 'write':
            return self._write
        return getattr(self._stderr, attr)

def filtered_print(string, stderr, kw):
    if kw == 0:
        if not string.__contains__("font"):
            stderr.write(string)
    else:
        if not string.__contains__("DEBUG") and\
                not string.__contains__("oebio.utils.sdk:upload_file") and\
                not string.__contains__("oebio.report:save_fig"):
            stderr.write(string)

sys.stderr = FilteredPrinter(filtered_print, sys.stderr, kw=0)
sys.stdout = FilteredPrinter(filtered_print, sys.stdout, kw=1)

NAME_MAP = {'shuting.ge': '葛淑婷',
 'shuang.fu': '付双',
 'yicheng.xu': '许以成',
 'longhai.ding': '丁龙海',
 'qianqian.kang': '康倩倩',
 'jing.pan': '潘静',
 'ayu.li': '李阿雨'}


def write_log(proj_path: Path, file="/data/hstore4/project/项目执行记录.txt"):
    try:
        FILE = Path(file)
        if not FILE.exists():
            with open(file=file, mode='w', encoding='utf8') as f:
                ...
        table = str(proj_path / "分析确认单" / "分析确认单.xlsx")
        base: dict = pd.read_excel(table, index_col=0, sheet_name='分析基本信息').fillna('NaN').to_dict()['value']

        proj_id = base.get("项目编号", "NaN")
        species = base.get("样本物种", "NaN")
        proj_type = base.get("项目类别", "NaN")
        user = getpass.getuser()
        CN_USER = NAME_MAP.get(user, user)
        # TODO: oecloudProduct
        sample_num = "1" if proj_type == '中药成分鉴定' else "6"
        comp_grp = "NaN" if proj_type == '中药成分鉴定' else "1"
        time = datetime.datetime.now().strftime("%Y/%m/%d")

        info = [proj_id, species, proj_type, CN_USER, sample_num, comp_grp, time, "\n"]
        log_line = "\t".join(info)
        with open(file=file, mode='a', encoding='utf8') as f:
            f.write(log_line)
    except Exception as e:
        print(e)


def add_loginfo(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        data = args[0]
        local = kwargs.get("local", 0)
        report_dir, proj_id = func(*args, **kwargs)
        # copy file to dest
        time = datetime.datetime.now().strftime("%Y%m")
        dest = Path(f"/data/nas/177/代谢/项目报告检查空间/预报告文件夹/{time}")
        if not dest.exists():
            dest.mkdir(parents=False,exist_ok=True)
        if local == 1 or local == "1":
            shutil.copy(f"{report_dir}/{proj_id}_report.html", f"/data/nas/177/代谢/项目报告检查空间/预报告文件夹/{time}/{proj_id}_report.html")
        # write log
        project_name, project_info = parse_project(data, local)
        proj_path = Path(data) / project_name
        write_log(proj_path, file="/data/hstore4/project/项目执行记录.txt")
    return wrapper


def wrap(p: list):
    return dict(zip(p[:10], [None]*10))


def filtered_print(string, stderr, kw):
    if kw == 0:
        if not string.__contains__("font"):
            stderr.write(string)
    else:
        if not string.__contains__("DEBUG"):
            stderr.write(string)

sys.stderr = FilteredPrinter(filtered_print, sys.stderr, kw=0)
sys.stdout = FilteredPrinter(filtered_print, sys.stdout, kw=1)


@functools.lru_cache(maxsize=1000)
def parse_meta(data) -> Tuple[list, pd.DataFrame]:
    """parse sample info"""
    mode_used = []

    SAMPLE_INFO = Path(data) / "sampleinfo.xlsx"
    sample_info = pd.read_excel(SAMPLE_INFO, sheet_name="样品信息", index_col=0)

    is_null = sample_info.isnull().any(axis=1)
    if not is_null["POS"]:
        mode_used.append("POS")
    if not is_null["NEG"]:
        mode_used.append("NEG")

    return mode_used, sample_info


@functools.lru_cache(maxsize=1000)
def parse_exp_col(data):
    """
    返回样本的实际列名
    """
    mode_used, sample_info = parse_meta(data)
    sample_info = sample_info.T
    mode = mode_used[0]
    tcm_cols = sample_info[mode].drop(index=["QI定性结果", "QI定量结果", "QI二级MSP"])\
        .apply(lambda x: x.rsplit(".", 1)[0]).tolist()
    return tcm_cols


@functools.lru_cache(maxsize=1000)
def parse_meta_blood(data) -> Tuple[list, pd.DataFrame]:
    """override parse sample info"""
    mode_used = []

    xls_path = Path(data) / "sampleinfo.xlsx"
    sample_info = pd.read_excel(str(xls_path), sheet_name="样品信息", index_col=[0, 1]).T

    is_null = sample_info.isnull().any(axis=1)
    if not is_null["POS"]:
        mode_used.append("POS")
    if not is_null["NEG"]:
        mode_used.append("NEG")

    return mode_used, sample_info


@functools.lru_cache(maxsize=1000)
def parse_exp_col_blood(data):
    """
    返回样本的实际列名
    """
    mode_used, sample_info = parse_meta_blood(data)
    mode = mode_used[0]
    blood_cols = sample_info["Blood"].drop(columns=["QI定性结果", "QI定量结果", "QI二级MSP"]).loc[mode, :].apply(
        lambda x: x.rsplit(".", 1)[0]).tolist()
    control_cols = sample_info["Control"].loc[mode, :].apply(lambda x: x.rsplit(".", 1)[0]).tolist()
    tcm_cols = sample_info["TCM"].loc[mode, :].apply(lambda x: x.rsplit(".", 1)[0]).tolist()

    return blood_cols, control_cols, tcm_cols


@functools.lru_cache(maxsize=1000)
def get_online_register_form(analysis_id: str) -> dict:
    resp = requests.post(BASE_URL, data=dict(analysis_id=analysis_id, token=TOKEN), timeout=10)
    if resp.status_code != 200:
        raise ValueError(f"不存在的分析编号：{analysis_id}")
    r: dict = resp.json()
    resp.close()

    key = r['分析基本信息']['key'].values()
    values = r['分析基本信息']['value'].values()
    parse_data = dict(zip(key, values))
    return parse_data


@functools.lru_cache(maxsize=1000)
def parse_project(data, local) -> Tuple[str, pd.DataFrame]:
    """parse project info"""
    PROJECT_INFO = Path(data) / "sampleinfo.xlsx"
    project_info = pd.read_excel(PROJECT_INFO, sheet_name="项目信息").applymap(
        lambda x: x.strip() if isinstance(x, str) else x)
    # 使用线上登记单
    if local == 1 or local == "1":
        aid: str
        obj = re.search('-a[a-z]', project_info.loc[0, "项目编号"])
        if isinstance(obj, re.Match):
            aid = project_info.loc[0, "项目编号"]
        else:
            aid = project_info.loc[0, "项目编号"] + "-aa"
        form = get_online_register_form(aid)
        project_path = Path(data).cwd().name + "-" + form.get("客户名称", "***") + "-" + project_info.loc[
            0, "项目名称"] + "结题报告"
        r = sp.run(["tool_GetAnalystInfo", "-ad", aid, "-sp", str(Path(data) / project_path / "分析确认单"), "-fn", "分析确认单.xlsx",
                    "--overwrite"])
        if r.returncode != 0:
            raise SaveOnlineRegisterTableError("保存线上登记单出错!")
    else:
        form = dict()

    project_info.loc[0, "客户名称"] = form.get("客户名称", project_info.loc[0, "客户名称"])
    project_info.loc[0, "联系人名称"] = form.get("联系人", project_info.loc[0, "联系人名称"])
    project_info.loc[0, "样本"] = form.get("样本类型", project_info.loc[0, "样本"])
    project_path = Path(data).cwd().name + "-" + project_info.loc[0, "客户名称"] + "-" + project_info.loc[
        0, "项目名称"] + "结题报告"
    return project_path, project_info


@add_loginfo
def oecloudTcm(data, version="/data/hstore1/database/report/2023-03-04~2023/src", company="LM", using_herb=False, cloud=False, local=0):

    from lmbio.basic.reportinfo import getreportinfo
    os.environ['OEBIO'] = version
    reportinfo = getreportinfo(sys_logo=company)
    if cloud:
        os.environ['CLOUD'] = "True"
        os.environ['CLOUDUSER'] = "lmbioinformatics@lumingbio.com"
        os.environ['PASSWORD'] = "Lmbio@123"


    src_path = Path("/data/hstore1/database/database/tcm/")

    project_name, project_info = parse_project(data, local)
    # 批量重命名
    data_path = str(Path(data) / project_name / "数据分析" / "4.定性定量结果" / "数据矩阵.xlsx")
    df = pd.read_excel(data_path, sheet_name=0).applymap(lambda x: x.strip() if isinstance(x, str) else x)
    new_names = pd.Series(data=[f"compound{x:0>5d}" for x in range(1, df.shape[0] + 1)])
    name_map = dict(zip(df["No."], new_names))
    fig_path = Path(data) / project_name / "数据分析" / "3.质谱图"
    formats = ("pdf", "png")
    for fm in formats:
        for p in fig_path.rglob(f"*.{fm}"):
            _ = name_map.get(p.stem, None)
            if _:
                dst = Path(p.parent) / (_ + f".{fm}")
                os.rename(p, dst)
            else:
                print(f"文件名{p.stem} 映射出错!")
    df["No."] = new_names
    df.to_excel(data_path, index=False)

    # 加载配置文件
    config = yaml.load(open(str(src_path / "tcm" / "config.yaml"), encoding="utf8"),
                       Loader=yaml.FullLoader)
    # 拷贝文件到实验内容
    report_xls = [path for path in (src_path / "tables").rglob("*.xlsx")] + [src_path / "tcm" / "数据矩阵字段说明.xlsx"]
    for xls in report_xls:
        shutil.copy(str(xls), f"{str(Path(data) / project_name / '实验内容')}/{xls.name}")

    # 拷贝UV data
    uv_path = [path for path in (Path(data) / "UVdata").rglob("*紫外吸收图*")]
    for uv in uv_path:
        shutil.copy(str(uv), f"{str(Path(data) / project_name / '数据分析' / '1.紫外吸收图')}/{uv.name}")

    # 设置项目信息
    config["header_info"]["项目名称"] = '中药成分鉴定分析'

    task_id = None
    try:
        task_id = project_info.loc[0, "任务单号"]
    except:
        task_id = project_info.loc[0, "项目编号"]

    config["header_info"]["任务单号"] = task_id if "-b" in task_id else task_id + "-b1"
    config["header_info"]["客户名称"] = project_info.loc[0, "客户名称"] if project_info.loc[0, "客户名称"] != '' and pd.notna(project_info.loc[0, "客户名称"]) else project_info.loc[0, "联系人名称"]
    config["header_info"]["联系人名称"] = project_info.loc[0, "联系人名称"] if project_info.loc[0, "联系人名称"] != '' and pd.notna(project_info.loc[0, "联系人名称"]) else project_info.loc[0, "客户名称"]
    config["header_info"]["项目编号"] = Path(data).cwd().name
    config["header_info"]["样本"] = project_info.loc[0, "样本"]
    if config["header_info"]["客户名称"] == config["header_info"]["联系人名称"]:
        config["header_info"].pop("客户名称")
    header_info = config.get("header_info")

    # OUT PATH
    report_dir = str(Path(data) / project_name)

    # 实例化报告对象
    report = oebioReport('中药成分鉴定分析', title='中药成分鉴定分析', header_info=header_info, oe_welcome = reportinfo.oe_weclome)

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
    # 中药成分鉴定
    son3.add_table(f"{report_dir}/{config['数据分析']['数据预处理']['定性定量结果']['tablepath2']}",
                   caption=config['数据分析']['数据预处理']['定性定量结果']['tablename2'], show_rows=10)
    son3.add_comment("""仅展示部分结果，详细结果请见支持文件：[中药成分鉴定数据矩阵.xlsx](数据分析/4.定性定量结果)""")
    son3.add_fig(f"{report_dir}/{config['数据分析']['数据预处理']['中药成分分类']['imgpath0']}",
                 caption=config['数据分析']['数据预处理']['中药成分分类']['imgname0'])

    # MS/MS
    son4 = a2.add_section("中药主要化学成分的EIC图及其与标准品比对情况",
                          description=config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['description'])
    # son4.add_fig(wrap(glob(f"{report_dir}/{config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['imgpath']}")),
    #              caption=config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['imgname'])
    # son4.add_comment("""仅展示部分结果，详细结果请见支持文件：[质谱图](数据分析/3.质谱图)""")

    # LuMet-TCM
    lumet_tcm = son4.add_section("中药原方成分的EIC图及其与标准品库LuMet-TCM的比对情况")
    lumet_tcm.add_fig(wrap(glob(f"{report_dir}/{config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与标准品库LuMet-TCM的比对情况']['imgpath']}")),
                 caption=config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与标准品库LuMet-TCM的比对情况']['imgname'])
    lumet_tcm.add_comment("""仅展示部分结果，详细结果请见支持文件：[LuMet-TCM](数据分析/3.质谱图/中药原方成分/LuMet-TCM)""")

    if using_herb:
        # HerbDB
        herb_db = son4.add_section("中药原方成分的EIC图及其与公共库HerbDB的比对情况")
        herb_db.add_fig(wrap(glob(f"{report_dir}/{config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与公共库HerbDB的比对情况']['imgpath']}")),
                     caption=config['数据分析']['数据预处理']['中药主要化学成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与公共库HerbDB的比对情况']['imgname'])
        herb_db.add_comment("""仅展示部分结果，详细结果请见支持文件：[HerbDB](数据分析/3.质谱图/中药原方成分/HerbDB)""")
    # pie
    son6 = a2.add_section("中药成分分类", description=config['数据分析']['数据预处理']['中药成分分类']['description'])
    son6.add_fig(f"{report_dir}/{config['数据分析']['数据预处理']['中药成分分类']['imgpath1']}",
                 caption=config['数据分析']['数据预处理']['中药成分分类']['imgname1'], description="为了显示美观，自动隐藏了占比<1%的数据标签")
    son6.add_fig(f"{report_dir}/{config['数据分析']['数据预处理']['中药成分分类']['imgpath2']}",
                 caption=config['数据分析']['数据预处理']['中药成分分类']['imgname2'], description="为了显示美观，自动隐藏了占比<1%的数据标签")

    # 公司简介
    reportinfo.companyinfo(report=report,
                           yamlpath=os.environ['OEBIO'])

    proj_id = config["header_info"]["项目编号"]
    # export
    report.write_to(f"{report_dir}/{proj_id}_report.html", zip_report_name=report_dir+f".zip")

    return report_dir, proj_id
    # # 添加与报告到 /public/项目报告检查空间/预报告文件夹/
    # pnum = Path(data).name # project_info.loc[0, "项目编号"]
    # pname = project_info.loc[0, "项目名称"]
    # # shutil.copy(f"{report_dir}/项目报告.html", f"/public/项目报告检查空间/预报告文件夹/{pnum}-{pname}-LCMS-1.0.0-分析预报告.html")
    # shutil.copy(f"{report_dir}/项目报告.html", f"/data/nas/177/代谢/项目报告检查空间/预报告文件夹/{pnum}-{pname}-LCMS-1.0.0-分析预报告.html")


@add_loginfo
def oecloudBlood(data,
                 version="/data/hstore1/database/report/2023-03-04~2023/src", company="LM",
                 using_product=False, using_herb=False, cloud=True, local=0, skip=False):

    project_name, project_info = parse_project(data, local)
    blood_cols, control_cols, tcm_cols = parse_exp_col_blood(data)
    midware = Path(data) / "midware"

    from lmbio.basic.reportinfo import getreportinfo
    os.environ['OEBIO'] = version
    reportinfo = getreportinfo(sys_logo=company)
    if cloud:
        os.environ['CLOUD'] = "True"
        os.environ['CLOUDUSER'] = "lmbioinformatics@lumingbio.com"
        os.environ['PASSWORD'] = "Lmbio@123"

    src_path = Path("/data/hstore1/database/database/tcm/")
    if not skip:
    # -------------------do before report-------------------
        # 批量重命名
        data_path = str(midware / "数据矩阵.xlsx")
        df = pd.read_excel(data_path, sheet_name=0).applymap(lambda x: x.strip() if isinstance(x, str) else x)
        new_names = pd.Series(data=[f"compound{x:0>5d}" for x in range(1, df.shape[0] + 1)])
        name_map = dict(zip(df["No."], new_names))
        fig_path = Path(data) / project_name / "数据分析" / "3.质谱图"
        formats = ("pdf", "png")
        for fm in formats:
            for p in fig_path.rglob(f"*.{fm}"):
                _ = name_map.get(p.stem, None)
                if _:
                    dst = Path(p.parent) / (_ + f".{fm}")
                    os.rename(p, dst)
                else:
                    print(f"文件名{p.stem} 映射出错!")
        df["No."] = new_names
        df.to_excel(data_path, index=False)

        # -------------------------------------------------------
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
        data_path = str(midware / "数据矩阵.xlsx")
        ingredient_path = str(Path(data) / project_name / "数据分析" / "4.定性定量结果" / "中药成分鉴定数据矩阵.xlsx")
        blood_path = str(Path(data) / project_name / "数据分析" / "4.定性定量结果" / "入血成分分析数据矩阵.xlsx")

        peaks = pd.read_excel(data_path, sheet_name=0).reset_index(drop=True)
        peaks["TCM均值（峰面积的比值%）"] = 100 * peaks["TCM均值"] / peaks["TCM均值"].sum()

        # rename
        peaks.rename(columns=ren_dict, inplace=True)
        #
        # -----------------------export products df----------------------------
        if using_product:
            product_path = str(Path(data) / project_name / "数据分析" / "4.定性定量结果" / "代谢产物数据矩阵.xlsx")
            raw_product_path = str(midware / "raw代谢产物数据矩阵.xlsx")
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
                col_order = ["No.", 'ID', 'Mid', 'Parent Compound', 'Transformations', 'Adducts', 'Formula', 'Score',
                             'Fragmentation Score',
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
                     *new_tcm_cols, 'TCM均值', 'TCM均值（峰面积的比值%）', "InChIKey", "SMILES", "HMDB", "METLIN", "LipidMaps", "KEGG",
                     "PubChem", "CAS", "中文名",
                     "中文大类", "中文子类", "英文大类", "英文子类", "货号", "纯度", 'Level']

        tcm_ingredients = peaks.query("DB != 'PRODUCT'")[used_cols]
        tcm_ingredients = tcm_ingredients.loc[tcm_ingredients[new_tcm_cols].mean(axis=1) >= 1000, :]
        tcm_ingredients["TCM均值（峰面积的比值%）"] = 100 * tcm_ingredients["TCM均值"] / tcm_ingredients["TCM均值"].sum()
        tcm_ingredients.sort_values(by=["Level"]).to_excel(ingredient_path, index=False)

        # -----------------------export blood df----------------------------
        peaks.drop(columns=["Compound ID", "DB"]).sort_values(by=["是否为入血成分", "Level"], ascending=[False, True]).to_excel(
            blood_path, index=False)
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
        shutil.copy(str(xls), f"{str(Path(data) / project_name / '实验内容')}/{xls.name}")

    # 拷贝UV data
    uv_path = [path for path in (Path(data) / "UVdata").rglob("*紫外吸收图*")]
    for uv in uv_path:
        shutil.copy(str(uv), f"{str(Path(data) / project_name / '数据分析' / '1.紫外吸收图')}/{uv.name}")

    # 设置项目信息
    config["header_info"]["项目名称"] = '中药入血成分分析'
    task_id = None
    try:
        task_id = project_info.loc[0, "任务单号"]
    except:
        task_id = project_info.loc[0, "项目编号"]
    config["header_info"]["任务单号"] = task_id if "-b" in task_id else task_id + "-b1"
    config["header_info"]["客户名称"] = project_info.loc[0, "客户名称"] if project_info.loc[0, "客户名称"] != '' and pd.notna(project_info.loc[0, "客户名称"]) else project_info.loc[0, "联系人名称"]
    config["header_info"]["联系人名称"] = project_info.loc[0, "联系人名称"] if project_info.loc[0, "联系人名称"] != '' and pd.notna(project_info.loc[0, "联系人名称"]) else project_info.loc[0, "客户名称"]
    config["header_info"]["项目编号"] = Path(data).cwd().name
    config["header_info"]["样本"] = project_info.loc[0, "样本"]
    if config["header_info"]["客户名称"] == config["header_info"]["联系人名称"]:
        config["header_info"].pop("客户名称")

    header_info = config.get("header_info")
    report_dir = str(Path(data) / project_name)

    # 实例化报告对象
    report = oebioReport('中药入血成分分析', title='中药入血成分分析', header_info=header_info, oe_welcome=reportinfo.oe_weclome)

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
                   caption=config['数据分析']['数据预处理']['定性定量结果']['tablename2'], show_rows=10)
    son3.add_comment("""仅展示部分结果，详细结果请见支持文件：[中药成分鉴定数据矩阵.xlsx](数据分析/4.定性定量结果)""")
    son3.add_fig(f"{report_dir}/{config['数据分析']['数据预处理']['中药成分分类']['imgpath0']}",
                 caption=config['数据分析']['数据预处理']['中药成分分类']['imgname0'])


    # 中药入血
    son3.add_table(f"{report_dir}/{config['数据分析']['数据预处理']['定性定量结果']['tablepath3']}",
                   caption=config['数据分析']['数据预处理']['定性定量结果']['tablename3'], show_rows=10)
    son3.add_comment("""仅展示部分结果，详细结果请见支持文件：[入血成分分析数据矩阵.xlsx](数据分析/4.定性定量结果)""")
    # 代谢产物
    if using_product:
        son3.add_table(f"{report_dir}/{config['数据分析']['数据预处理']['中药代谢产物网络图']['tablepath1']}",
                       caption=config['数据分析']['数据预处理']['中药代谢产物网络图']['tablename1'])
        son3.add_table(f"{report_dir}/{config['数据分析']['数据预处理']['定性定量结果']['tablepath4']}",
                       caption=config['数据分析']['数据预处理']['定性定量结果']['tablename4'], show_rows=10)
        son3.add_comment("""仅展示部分结果，详细结果请见支持文件：[代谢产物数据矩阵.xlsx](数据分析/4.定性定量结果)""")
    # MS/MS
    # 中药原方成分
    son4 = a2.add_section("中药原方成分的EIC图及其与标准品比对情况")
    # LuMet-TCM
    lumet_tcm = son4.add_section("中药原方成分的EIC图及其与标准品库LuMet-TCM的比对情况")

    lumet_tcm.add_fig(wrap(glob(f"{report_dir}/{config['数据分析']['数据预处理']['中药原方成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与标准品库LuMet-TCM的比对情况']['imgpath']}")),
                 caption=config['数据分析']['数据预处理']['中药原方成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与标准品库LuMet-TCM的比对情况']['imgname'])
    lumet_tcm.add_comment("""仅展示部分结果，详细结果请见支持文件：[LuMet-TCM](数据分析/3.质谱图/中药原方成分/LuMet-TCM)""")
    if using_herb:
        # HerbDB
        herb_db = son4.add_section("中药原方成分的EIC图及其与公共库HerbDB的比对情况")
        herb_db.add_fig(wrap(glob(f"{report_dir}/{config['数据分析']['数据预处理']['中药原方成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与公共库HerbDB的比对情况']['imgpath']}")),
                     caption=config['数据分析']['数据预处理']['中药原方成分的EIC图及其与标准品比对情况']['中药原方成分的EIC图及其与公共库HerbDB的比对情况']['imgname'])
        herb_db.add_comment("""仅展示部分结果，详细结果请见支持文件：[HerbDB](数据分析/3.质谱图/中药原方成分/HerbDB)""")
    # 中药入血成分
    son41 = a2.add_section("中药入血成分的EIC图及其与标准品比对情况",
                           description=config['数据分析']['数据预处理']['中药入血成分的EIC图及其与标准品比对情况']["description"])
    son411 = son41.add_section("中药入血原型成分")
    # 中药原型
    son411.add_fig(wrap(glob(f"{report_dir}/{config['数据分析']['数据预处理']['中药入血成分的EIC图及其与标准品比对情况']['中药入血原型成分的EIC图及其与参考化合物比对情况']['imgpath']}")),
                  caption=config['数据分析']['数据预处理']['中药入血成分的EIC图及其与标准品比对情况']['中药入血原型成分的EIC图及其与参考化合物比对情况']['imgname'])
    son411.add_comment("""仅展示部分结果，详细结果请见支持文件：[中药入血原型成分](数据分析/3.质谱图/中药入血成分/中药原型)""")
    # 中药代谢产物
    if using_product:
        son412 = son41.add_section("中药代谢产物")
        # product
        son412.add_fig(wrap(glob(f"{report_dir}/{config['数据分析']['数据预处理']['中药入血成分的EIC图及其与标准品比对情况']['中药代谢产物的EIC图及其与参考化合物比对情况']['imgpath']}")),
                     caption=config['数据分析']['数据预处理']['中药入血成分的EIC图及其与标准品比对情况']['中药代谢产物的EIC图及其与参考化合物比对情况']['imgname'])
        son412.add_comment("""仅展示部分结果，详细结果请见支持文件：[中药代谢产物](数据分析/3.质谱图/中药入血成分/中药代谢产物)""")
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
        son5.add_fig(wrap(glob(f"{report_dir}/{config['数据分析']['数据预处理']['中药代谢产物网络图']['imgpath']}")),
                     caption=config['数据分析']['数据预处理']['中药代谢产物网络图']['imgname'])
        son5.add_comment("""仅展示部分结果，详细结果请见支持文件：[中药代谢产物网络图](数据分析/6.代谢网络图)""")

    # ---------------------内源代谢物分析 start---------------
    def mv_top(p):
        p = Path(p)
        return str(p.relative_to(p.parts[0]))
    endo_dir = f'{report_dir}/数据分析/*差异分析结果'
    if list(glob(endo_dir)):
        report.add_yaml_config(os.path.join("/data/nas/174/研发-项目部-why/生信研发-生信研发部-ljw/database/report/2023-03-04~2023/全谱代谢-云平台", 'Report.yaml'))
        Analysis_module = analysis.add_section('内源代谢物分析')
        # 多元统计分析
        if list(glob(f'{endo_dir}/*多元统计分析/')):
            diff_ana1 = Analysis_module.add_section('多元统计分析', link=mv_top(list(glob(f'{endo_dir}/*数据矩阵/数据矩阵.xlsx'))[0]))
        if list(glob(f'{endo_dir}/*多元统计分析/PCA')):
            link = list(glob(f'{endo_dir}/*多元统计分析/PCA'))
            diff_ana1_1 = diff_ana1.add_section("主成分分析(PCA)", link=mv_top(link[0]))
            if len(list(glob(f'{endo_dir}/*多元统计分析/PCA/PCA*.png'))) == 0:
                diff_ana1_plot = diff_ana1_1.add_plot('PCA-example.png',
                                                      caption='PCA图',
                                                      description='图片说明：不同颜色形状组合表示不同样本分组，椭圆区域代表95%的置信区间。')
            else:
                diff_ana1_plot = diff_ana1_1.add_plot([(f, '') for f in glob(f'{endo_dir}/*多元统计分析/PCA/PCA-score-*.png')][:10],
                                                      caption='PCA图',
                                                      description='图片说明：横坐标PC1为第一主成分解释率，纵坐标PC2为第二主成分解释率，图中每个点代表一个样品,不同颜色形状组合表示不同样本分组，椭圆区域代表95%的置信区间。')
                diff_ana1_plot = diff_ana1_1.add_plot([(f, '') for f in glob(f'{endo_dir}/*多元统计分析/PCA/PCA-3Dscore-*.png')][:10],
                                                      caption='PCA-3D图',
                                                      description='图片说明：横坐标PC1为第一主成分解释率，纵坐标PC2为第二主成分解释率，图中每个点代表一个样品,不同颜色形状组合表示不同样本分组。')
        if list(glob(f'{endo_dir}/*多元统计分析/PLS-DA')):
            link = list(glob(f'{endo_dir}/*多元统计分析/PLS-DA'))
            diff_ana1_2 = diff_ana1.add_section("偏最小二乘-判别分析(PLS-DA)", link=mv_top(link[0]))
            if len(list(glob(f'{endo_dir}/*多元统计分析/PLS-DA/PLS-DA*.png'))) == 0:
                diff_ana1_plot = diff_ana1_2.add_plot('PLS-example.png',
                                                      caption='PLS-DA图',
                                                      description='图片说明：解释率R2Y (cum)和预测率Q2 (cum)，两者越接近1，说明PLS-DA模型能更好地解释和预测两组样本之间的差异，代表模型预测能力越好。')
            else:
                diff_ana1_plot = diff_ana1_2.add_plot([(f, '') for f in glob(f'{endo_dir}/*多元统计分析/PLS-DA/PLS-DA-score-*.png')][:10],
                                                      caption='PLS-DA图',
                                                      description='图片说明：解释率R2Y (cum)和预测率Q2 (cum)，两者越接近1，说明PLS-DA模型能更好地解释和预测两组样本之间的差异，代表模型预测能力越好。')
                diff_ana1_plot = diff_ana1_2.add_plot([(f, '') for f in glob(f'{endo_dir}/*多元统计分析/PLS-DA/PLS-DA-3Dscore-*.png')][:10],
                                                      caption='PLS-DA-3D图',
                                                      description='图片说明：解释率R2Y (cum)和预测率Q2 (cum)，两者越接近1，说明PLS-DA模型能更好地解释和预测两组样本之间的差异，代表模型预测能力越好。')
                diff_ana1_plot = diff_ana1_2.add_plot([(f, '') for f in glob(f'{endo_dir}/*多元统计分析/PLS-DA/PLS-DA-loading-*.png')][:10],
                                                      caption='PLS-DA-loading图',
                                                      description='图片说明：使用载荷图可标识代谢物对比较组的影响强度。载荷范围可以为 -1 到 1。接近于 -1 或 1 的载荷表明变量对分量影响非常强。接近于 0 的载荷表明变量对分量的影响很弱。通过载荷有助于根据变量表示每个分量的特征。')
                diff_ana1_plot = diff_ana1_2.add_plot([(f, '') for f in glob(f'{endo_dir}/*多元统计分析/PLS-DA/PLS-DA-splot-*.png')][:10],
                                                      caption='PLS-DA-splot图',
                                                      description='图片说明：Splot图横坐标是代谢物对比较组影响的特征值，纵坐标是样品得分与代谢物之间的相关性。由于特征值与相关性同正同负，因此可视化的图中所有的点全部分布在第一和第三象限，类似于S形，被称为Splot，越靠近右上角和左下角的代谢物表示其差异越显著。')

        if list(glob(f'{endo_dir}/*多元统计分析/OPLS-DA')):
            link = list(glob(f'{endo_dir}/*多元统计分析/OPLS-DA'))
            diff_ana1_3 = diff_ana1.add_section("正交偏最小二乘方-判别分析(OPLS-DA)", link=mv_top(link[0]))
            if len(list(glob(f'{endo_dir}/*多元统计分析/OPLS-DA/OPLS-DA*.png'))) == 0:
                diff_ana1_plot = diff_ana1_3.add_plot('OPLS-example.png',
                                                      caption='OPLS-DA图',
                                                      description='图片说明：OPLS-DA将组间差异最大化反映在t1上，所以从t1上能直接区分组间变异，而正交主成分上则反映了组内变异。')
            else:
                diff_ana1_plot = diff_ana1_3.add_plot([(f, '') for f in glob(f'{endo_dir}/*多元统计分析/OPLS-DA/OPLS-DA-score-*.png')][:10],
                                                      caption='OPLS-DA图',
                                                      description='图片说明：OPLS-DA将组间差异最大化反映在t1上，所以从t1上能直接区分组间变异，而正交主成分上则反映了组内变异。')
                diff_ana1_plot = diff_ana1_3.add_plot(
                    [(f, '') for f in glob(f'{endo_dir}/*多元统计分析/OPLS-DA/OPLS-DA-3Dscore-*.png')][:10],
                    caption='OPLS-DA-3D图',
                    description='图片说明：OPLS-DA将组间差异最大化反映在t1上，所以从t1上能直接区分组间变异，而正交主成分上则反映了组内变异。')
                diff_ana1_plot = diff_ana1_3.add_plot(
                    [(f, '') for f in glob(f'{endo_dir}/*多元统计分析/OPLS-DA/OPLS-DA-loading-*.png')][:10],
                    caption='OPLS-DA-loading图',
                    description='图片说明：使用载荷图可标识代谢物对比较组的影响强度。载荷范围可以为 -1 到 1。接近于 -1 或 1 的载荷表明变量对分量影响非常强。接近于 0 的载荷表明变量对分量的影响很弱。通过载荷有助于根据变量表示每个分量的特征。')
                diff_ana1_plot = diff_ana1_3.add_plot([(f, '') for f in glob(f'{endo_dir}/*多元统计分析/OPLS-DA/OPLS-DA-splot-*.png')][:10],
                                                      caption='OPLS-DA-splot图',
                                                      description='图片说明：Splot图横坐标是代谢物对比较组影响的特征值，纵坐标是样品得分与代谢物之间的相关性。由于特征值与相关性同正同负，因此可视化的图中所有的点全部分布在第一和第三象限，类似于S形，被称为Splot，越靠近右上角和左下角的代谢物表示其差异越显著。')

        if list(glob(f'{endo_dir}/*多元统计分析/permutation')):
            link = list(glob(f'{endo_dir}/*多元统计分析/permutation'))
            diff_ana1_4 = diff_ana1.add_section("permutation", link=mv_top(link[0]))
            if len(list(glob(f'{endo_dir}/*多元统计分析/permutation/permutation*.png'))) == 0:
                diff_ana1_plot = diff_ana1_4.add_plot('Permutation-example.png',
                                                      caption='permutation图',
                                                      description='图片说明：模型有效性评估，主要参考两条标准：1.左边的所有绿色Q2值都低于右边的原始点；或者2.Q2点的绿色回归线与纵轴(左侧)相交于或低于零。')
            else:
                diff_ana1_plot = diff_ana1_4.add_plot(
                    [(f, '') for f in glob(f'{endo_dir}/*多元统计分析/permutation/permutation*.png')][:10],
                    caption='permutation图',
                    description='图片说明：模型有效性评估，主要参考两条标准：1.左边的所有绿色Q2值都低于右边的原始点；或者2.Q2点的绿色回归线与纵轴(左侧)相交于或低于零。')

        if list(glob(f'{endo_dir}/*多元统计分析/summarydata.xls')):
            diff_ana1_table = diff_ana1.add_table(f'{endo_dir}/*多元统计分析/summarydata.xls',
                                                  caption='模型参数',
                                                  description='''表头解释如下：
    
[1]. Modetype：建立的多元统计分析模型；

[2]. PRE：代表建模时主成分的个数；

[3]. ORT：代表建模时正交成分的个数；

[4]. N：代表建模时样本的个数；

[5]. R2X（cum）：代表多元统计分析建模时，在X轴方向模型的累积解释率（或可以理解为X轴方向保留原始数据信息百分比的平方），cum表示几个主成分累积的结果；

[6]. R2Y（cum）：代表在Y轴方向模型的累积解释率（或可以理解为Y轴方向保留原始数据信息百分比的平方）；

[7]. Q2（cum）：代表模型的累积预测率；

[8]. R2、Q2：响应排序检验的参数，用来衡量模型是否过拟合；''',show_rows=4)

        diff_ana = Analysis_module.add_section('差异比较分析', link=mv_top(list(glob(f'{endo_dir}/*数据矩阵/数据矩阵.xlsx'))[0]))
        # 差异代谢物
        if list(glob(f'{endo_dir}/*差异代谢物')):
            link = list(glob(f'{endo_dir}/*差异代谢物'))
            diff_ana3 = diff_ana.add_section("差异代谢物筛选", link=mv_top(link[0]))
            if list(glob(f'{endo_dir}/*差异代谢物/差异表达矩阵.xlsx')):
                diff_ana3_table = diff_ana3.add_table(f'{endo_dir}/*差异代谢物/差异表达矩阵.xlsx',
                                                      caption='差异代谢物',
                                                      description='''表头解释如下：

[1]. VIP：变量权重值，来自OPLS-DA模型的VIP值，VIP越大，说明该变量对分组的贡献越大;

[2]. p-value：T检验的结果，用来评价变量在两组样本之间差异是否显著，p<0.05表示显著;

[3]. q-value：错误发现率（FDR），用于控制差异分析中假阳性结果的比例,采用Benjamini-Hochberg方法计算;

[4]. FoldChange：代谢物在两组样本中平均表达量的比值;

[5]. log2FoldChange：对FC值取log2对数，正值表示上调，负值表示下调;

[6]. Average （#）：#组代谢产物平均表达量''', show_rows=4)

            if list(glob(f'{endo_dir}/*差异代谢物/Heatmap/heatmap-*.xls')):
                link1 = list(glob(f'{endo_dir}/*差异代谢物/Heatmap'))
                link2 = list(glob(f'{endo_dir}/*差异代谢物/Heatmap_top50/'))
                diff_ana3_1 = diff_ana3.add_section("差异代谢物聚类热图", link1=mv_top(link1[0]), link2=mv_top(link2[0]),
                                                    link0=mv_top(list(glob(f'{endo_dir}/*差异代谢物/*Heatmap/heatmap*.xls'))[0]))
                if len(list(glob(f'{endo_dir}/*差异代谢物/*Heatmap_top50/*cluster_none/*.png'))) == 0:
                    diff_ana3_plot = diff_ana3_1.add_plot('heatmap-example.png',
                                                          caption='聚类热图(heatmap)',
                                                          description='图片说明：横坐标表示样本名称，纵坐标表示差异代谢物。颜色从蓝到红表示代谢物的表达丰度从低到高，即越红表示差异代谢物的表达丰度越高。')
                else:
                    diff_ana3_plot = diff_ana3_1.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*差异代谢物/*Heatmap_top50/*cluster_none/*.png')][:10],
                        caption='聚类热图(heatmap)',
                        description='图片说明：横坐标表示样本名称，纵坐标表示差异代谢物。颜色从蓝到红表示代谢物的表达丰度从低到高，即越红表示差异代谢物的表达丰度越高。')

            if list(glob(f'{endo_dir}/*差异代谢物/Volcano/volcano-*.xls')):
                link = list(glob(f'{endo_dir}/*差异代谢物/Volcano'))
                diff_ana3_2 = diff_ana3.add_section("差异代谢物火山图", link=mv_top(link[0]),
                                                    link0=mv_top(list(glob(f'{endo_dir}/*差异代谢物/Volcano/volcano*.xls'))[0]))
                if len(list(glob(f'{endo_dir}/*差异代谢物/Volcano/*.png'))) == 0:
                    diff_ana3_plot = diff_ana3_2.add_plot('volcano-example-2.png',
                                                          caption='火山图(volcano)',
                                                          description='图片说明：其中红色圆点代表在实验组中显著上调的差异代谢物，蓝色圆点代表显著下调的差异代谢物，灰色点代表不显著差异的代谢产物。图中每个点代表一个代谢物，横坐标是两组比对的log2（FC）值，纵坐标为-log10（p-value）值，红色点为显著上调的差异代谢物（p<0.05, VIP>1 且FC>1），蓝色点为显著下调的差异代谢物（p<0.05, VIP>1 且FC<1）')
                else:
                    diff_ana3_plot = diff_ana3_2.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*差异代谢物/Volcano/*.png')][:10],
                        caption='火山图(volcano)',
                        description='图片说明：其中红色圆点代表在实验组中显著上调的差异代谢物，蓝色圆点代表显著下调的差异代谢物，灰色点代表不显著差异的代谢产物。图中每个点代表一个代谢物，横坐标是两组比对的log2（FC）值，纵坐标为-log10（p-value）值，红色点为显著上调的差异代谢物（p<0.05, VIP>1 且FC>1），蓝色点为显著下调的差异代谢物（p<0.05, VIP>1 且FC<1）')

            if list(glob(f'{endo_dir}/*差异代谢物/Violin_top50/*/*.png')):
                link = list(glob(f'{endo_dir}/*差异代谢物/Violin_top50'))
                diff_ana3_3 = diff_ana3.add_section("差异代谢物小提琴图", link=mv_top(link[0]), link0=
                mv_top(list(glob(f'{endo_dir}/*差异代谢物/Violin_top50/violin_top50-*.xls'))[0]))
                if len(list(glob(f'{endo_dir}/*差异代谢物/Violin_top50/*/Violin-*.png'))) == 0:
                    diff_ana3_plot = diff_ana3_3.add_plot('Violin_top50-example.png',
                                                          caption='小提琴图(Violin_top50)')
                else:
                    diff_ana3_plot = diff_ana3_3.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*差异代谢物/Violin_top50/*/Violin-*.png')][:10],
                        caption='小提琴图(Violin_top50)')

            if list(glob(f'{endo_dir}/*差异代谢物/Boxplot_top50/*/*.png')):
                link = list(glob(f'{endo_dir}/*差异代谢物/Boxplot_top50'))
                diff_ana3_3 = diff_ana3.add_section("差异代谢物箱型图", link=mv_top(link[0]),
                                                    link0=mv_top(list(glob(f'{endo_dir}/*差异代谢物/Boxplot_top50/*.xls'))[0]))
                if len(list(glob(f'{endo_dir}/*差异代谢物/Boxplot_top50/*/*.png'))) == 0:
                    diff_ana3_plot = diff_ana3.add_plot('Boxplot_top50-example.png',
                                                        caption='箱型图(boxchart-top50)')
                else:
                    diff_ana3_plot = diff_ana3.add_plot(f'{endo_dir}/*差异代谢物/Boxplot_top50/*/*.png',
                                                        caption='箱型图(Boxplot_top50)')

            if list(glob(f'{endo_dir}/*差异代谢物/Zscore_top20/*.png')):
                link = list(glob(f'{endo_dir}/*差异代谢物/Zscore_top20'))
                diff_ana3_4 = diff_ana3.add_section("Z-Score标准化", link=mv_top(link[0]))
                if len(list(glob(f'{endo_dir}/*差异代谢物/Zscore_top20/*.png'))) == 0:
                    diff_ana3_plot = diff_ana3_4.add_plot('z-score-top20-example.png',
                                                          caption='Zscore_top20',
                                                          description='图片说明：在 z-score 图中，横轴表示 z-score 值范围，也就是该样本相对于均值的偏离程度，纵轴表示代谢物名称，每个数据点代表每个样本中的该代谢物。')
                else:
                    diff_ana3_plot = diff_ana3_4.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*差异代谢物/Zscore_top20/*.png')][:10],
                        caption='Zscore_top20',
                        description='图片说明：在 z-score 图中，横轴表示 z-score 值范围，也就是该样本相对于均值的偏离程度，纵轴表示代谢物名称，每个数据点代表每个样本中的该代谢物。')
            if list(glob(f'{endo_dir}/*差异代谢物/Venn/*.png')):
                link = list(glob(f'{endo_dir}/*差异代谢物/Venn'))
                diff_ana3_5 = diff_ana3.add_section("Venn", link=mv_top(link[0]))
                # diff_ana3_5_1 = diff_ana3_5.add_section("Venn",link = link[0])
                if len(list(glob(f'{endo_dir}/*差异代谢物/Venn/*.png'))) == 0:
                    diff_ana3_plot = diff_ana3_5.add_plot('Venn.png', caption='venn图')
                else:
                    diff_ana3_plot1 = diff_ana3_5.add_plot(f'{endo_dir}/*差异代谢物/Venn/Venn*.png', caption='venn图')

                if list(glob(f'{endo_dir}/*差异代谢物/Venn/Upset*.png')):
                    diff_ana3_5_2 = diff_ana3_5.add_section("Upset", link=mv_top(link[0]))
                    diff_ana3_5_2_plot2 = diff_ana3_5_2.add_plot(f'{endo_dir}/*差异代谢物/Venn/Upset*.png', caption='Upset图')

                if list(glob(f'{endo_dir}/*差异代谢物/Venn/Flower*.png')):
                    diff_ana3_5_3 = diff_ana3_5.add_section("Flower", link=mv_top(link[0]))
                    diff_ana3_5_3_plot3 = diff_ana3_5_3.add_plot(f'{endo_dir}/*差异代谢物/Venn/Flower*.png',
                                                                 caption='Flower图')

            if list(glob(f'{endo_dir}/*差异代谢物/Lolipopmap/*.png')):
                link = list(glob(f'{endo_dir}/*差异代谢物/Lolipopmap'))
                diff_ana3_6 = diff_ana3.add_section("Lolipopmap", link=mv_top(link[0]))
                if len(list(glob(f'{endo_dir}/*差异代谢物/Lolipopmap/*.png'))) == 0:
                    diff_ana3_plot = diff_ana3_6.add_plot('Lolipopmap.png',
                                                          caption='Lolipopma图')
                else:
                    diff_ana3_plot = diff_ana3_6.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*差异代谢物/Lolipopmap/*.png')][:10],
                        caption='Lolipopmap图',
                        description='图片说明：图中横坐标为log2（FoldChange），纵坐标为差异代谢物，分别取上调和下调差异代谢物中VIP值最大的10个差异代谢物，共计20个差异代谢物绘制图形，红色表示上调，蓝色表示下调，星号表示差异代谢的显著性（`*`表示显著性小于0.05、大于0.01；`**`表示显著性小于0.01、大于0.001；`***`表示显著性小于0.001、大于0.0001；`****`表示显著性小于0.0001），圆点大小由VIP值决定')

            if list(glob(f'{endo_dir}/*差异代谢物/相关性分析')):
                link = list(glob(f'{endo_dir}/*差异代谢物/相关性分析'))
                diff_ana4 = diff_ana.add_section("相关性分析", link=mv_top(link[0]), link0=
                mv_top(list(glob(f'{endo_dir}/*差异代谢物/相关性分析/Correlation_expression_top20-*.xls'))[0]))

            if list(glob(f'{endo_dir}/*差异代谢物/相关性分析/Correlation_top20')):
                link = list(glob(f'{endo_dir}/*差异代谢物/相关性分析/Correlation_top20'))
                diff_ana4_1 = diff_ana4.add_section("相关性图", link=mv_top(link[0]))
                if len(list(glob(f'{endo_dir}/*差异代谢物/相关性分析/Correlation_top20/*.png'))) == 0:
                    diff_ana4_plot = diff_ana4.add_plot('correlation-top20-example1.png',
                                                        caption='相关性图(correlation)',
                                                        description='图片说明：相关性分析使用Pearson相关系数衡量两个代谢物之间的线性相关程度。红色表示正相关，蓝色表示负相关。圆点越大，表示两个变量之间的相关性系数越大。')
                else:
                    diff_ana4_plot = diff_ana4.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*差异代谢物/相关性分析/Correlation_top20/*.png')][:10],
                        caption='相关性图(correlation)',
                        description='图片说明：相关性分析使用Pearson相关系数衡量两个代谢物之间的线性相关程度。红色表示正相关，蓝色表示负相关。圆点越大，表示两个变量之间的相关性系数越大。')

            if list(glob(f'{endo_dir}/*差异代谢物/相关性分析/Corrnetwork_top20')):
                link = list(glob(f'{endo_dir}/*差异代谢物/相关性分析/Corrnetwork_top20'))
                diff_ana4_2 = diff_ana4.add_section("相关性网络图", link=mv_top(link[0]))
                if len(list(glob(f'{endo_dir}/*差异代谢物/相关性分析/Corrnetwork_top20/*.png'))) == 0:
                    diff_ana4_plot = diff_ana4.add_plot('cornetwork-top20-example.png',
                                                        caption='相关性网络图(cornetwork)',
                                                        description='图片说明：图中圆形为差异代谢物，形状大小为连接数量， 形状之间的连线代表差异代谢物与差异代谢物之间的关联，连线的粗细代表关联性的高低程度，连线为红色为正相关，连线为蓝色为负相关。')
                else:
                    diff_ana4_plot = diff_ana4.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*差异代谢物/相关性分析/Corrnetwork_top20/*.png')][:10],
                        caption='相关性网络图(cornetwork)',
                        description='图片说明：图中圆形为差异代谢物，形状大小为连接数量， 形状之间的连线代表差异代谢物与差异代谢物之间的关联，连线的粗细代表关联性的高低程度，连线为红色为正相关，连线为蓝色为负相关。')

            # 代谢通路富集
            diff_ana_Term = Analysis_module.add_section("通路富集分析", link=mv_top(link[0]))
            if list(glob(f'{endo_dir}/*通路富集')):
                link = list(glob(f'{endo_dir}/*通路富集'))
                diff_ana5 = diff_ana_Term.add_section("KEGG通路富集分析", link=mv_top(link[0]),
                                                      link0=mv_top(list(glob(f'{endo_dir}/*通路富集/KEGG/*/*-Total.xls'))[0]))

            if list(glob(f'{endo_dir}/*通路富集/*KEGG/*/*Total.xls')):
                listdiff = list(glob(f'{endo_dir}/*通路富集/*KEGG/*/*Total.xls'))
                diff_ana5_1 = diff_ana5.add_table(listdiff[0], show_rows=10, description='''表头解释如下：

[1]. id：KEGG pathway 编号

[2]. Classification_level1：KEGG第一层级分类

[3]. Classification_level2：KEGG第二层级分类

[4]. Term：该通路的描述

[5]. ListHits：该通路中总差异代谢物数

[6]. ListTotal：注释到KEGG的总差异代谢物数

[7]. PopHits：注释到该通路中的所有代谢物个数

[8]. PopTotal：注释到KEGG的总代谢物数

[9]. p-value：富集显著性pvalue值，P≤0.05表示显著富集

[10]. padjust：校正后的pvalue值

[11]. Enrichment_score：富集打分

[12]. Substances： 属于该条通路的代谢物
        ''')

                # level3柱状图
                diff_ana5_1 = diff_ana5.add_section('KEGG Level3')
                if len(list(glob(f'{endo_dir}/*通路富集/KEGG/*/KEGG.Bar*.png'))) == 0:
                    diff_ana5_plot = diff_ana5_1.add_plot('KEGG Level3-example.png',
                                                          caption='KEGG Level3',
                                                          description='图片说明：横坐标为每条通路的-log10 p值，纵坐标为不同通路名称，柱子上的数字为注释到该通路的差异代谢物数量，柱子不同颜色代表不同的KEGG Level1信息。')
                else:
                    diff_ana5_plot = diff_ana5_1.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*通路富集/KEGG/*/KEGG.Bar*.png')][:10],
                        caption='通路富集分析图',
                        description='图片说明：横坐标为每条通路的-log10 p值，纵坐标为不同通路名称，柱子上的数字为注释到该通路的差异代谢物数量，柱子不同颜色代表不同的KEGG Level1信息。')
                # kegg富集分析气泡
                diff_ana6 = diff_ana5.add_section("KEGG气泡图1")
                if len(list(glob(f'{endo_dir}/*通路富集/KEGG/*/KEGG.Bubble*.png'))) == 0:
                    diff_ana6_plot = diff_ana6.add_plot('KEGG.Bubble-example.png',
                                                        caption='KEGG气泡图',
                                                        description='图片说明：图中横坐标Enrichment Score为富集分值，纵坐标为top20的通路信息。气泡越大的通路包含的差异代谢物越多，气泡颜色由蓝-红变化，其富集pvalue值越小，显著程度越大。')
                else:
                    diff_ana6_plot = diff_ana6.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*通路富集/KEGG/*/KEGG.Bubble*.png')][:10],
                        caption='KEGG气泡图',
                        description='图片说明：图中横坐标Enrichment Score为富集分值，纵坐标为top20的通路信息。气泡越大的通路包含的差异代谢物越多，气泡颜色由蓝-红变化，其富集pvalue值越小，显著程度越大')

                # kegg富集分析和弦图
                diff_ana7 = diff_ana5.add_section("KEGG富集分析和弦图")
                if len(list(glob(f'{endo_dir}/*通路富集/KEGG/*/KEGG.chord.png'))) == 0:
                    diff_ana7_plot = diff_ana7.add_plot('KEGG.chord-example.png',
                                                        caption='KEGG富集分析和弦图',
                                                        description='图片说明：左面为代谢KEGG：C号，红色表示上调，蓝色表示下调，右面为所选KEGG通路。')
                else:
                    diff_ana7_plot = diff_ana7.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*通路富集/KEGG/*/KEGG.chord*.png')][:10],
                        caption='KEGG富集分析和弦图',
                        description='图片说明：左面为代谢KEGG：C号，红色表示上调，蓝色表示下调，右面为所选KEGG通路。')

                # kegg富集分析圈图
                diff_ana8 = diff_ana5.add_section("KEGG富集分析圈图")
                if len(list(glob(f'{endo_dir}/*通路富集/KEGG/*/KEGG.circos.png'))) == 0:
                    diff_ana8_plot = diff_ana8.add_plot('KEGG.circos-example.png',
                                                        caption='KEGG富集分析圈图',
                                                        description='''图片说明:从外到内共4圈：

第一圈：富集的分类，圈外为代谢物数目的坐标尺，不同的颜色代表不同分类；

第二圈：背景代谢中该分类的数目以及p-value。代谢物越多条形越长，值越小颜色越红，越大越蓝；      

第三圈：上下调代谢比例条形图，浅红色代表上调代谢物比例，浅蓝色代表下调代谢物比例；下方显示具体的数值；      

第四圈：各分类的RichFactor值，背景辅助线每个小格表示0.2。
        ''')
                else:
                    diff_ana8_plot = diff_ana8.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*通路富集/KEGG/*/KEGG.circos.png')][:10],
                        caption='KEGG富集分析圈图',
                        description='图片说明:从外到内共4圈：第一圈：富集的分类，圈外为代谢物数目的坐标尺，不同的颜色代表不同分类；第二圈：背景代谢中该分类的数目以及p-value。代谢物越多条形越长，值越小颜色越红，越大越蓝；第三圈：上下调代谢比例条形图，浅红色代表上调代谢物比例，浅蓝色代表下调代谢物比例；下方显示具体的数值；第四圈：各分类的RichFactor值，背景辅助线每个小格表示0.2。')

                # kegg上下调对比图
                diff_ana9 = diff_ana5.add_section("KEGG上下调对比图")
                if len(list(glob(f'{endo_dir}/*通路富集/KEGG/*/Up_Down.Comparison.png'))) == 0:
                    diff_ana9_plot = diff_ana9.add_plot('Up_Down.Comparison-example.png',
                                                        caption='KEGG上下调对比图',
                                                        description='图片说明：横坐标为该通路前景差异代谢物数量与背景代谢物数量之比（ListHits/TotalHits），纵坐标为KEGG通路名称，有的通路会在上下调均显著富集。')
                else:
                    diff_ana9_plot = diff_ana9.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*通路富集/KEGG/*/Up_Down.Comparison.png')][:10],
                        caption='KEGG上下调对比图',
                        description='图片说明：横坐标为该通路前景差异代谢物数量与背景代谢物数量之比（ListHits/TotalHits），纵坐标为KEGG通路名称，有的通路会在上下调均显著富集。')

                # kegg level2分布图
                diff_ana10 = diff_ana5.add_section("KEGG Level2水平分布图")
                if len(list(glob(f'{endo_dir}/*通路富集/KEGG/*/Up_vs_Down.KEGG_Classification.png'))) == 0:
                    diff_ana10_plot = diff_ana10.add_plot('Up_vs_Down.KEGG_Classification-example.png',
                                                          caption='KEGG Level2水平分布图',
                                                          description='图片说明：横坐标是注释到各Level2代谢通路的上调（下调）差异代谢物和所有注释到KEGG通路的上调（下调）差异表达代谢物总数的比值（%），纵轴表示Level2 pathway的名称，柱子右边数字代表注释到该Level2 pathway下的上调（下调）差异代谢物数量。')
                else:
                    diff_ana10_plot = diff_ana10.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*通路富集/KEGG/*/Up_vs_Down.KEGG_Classification.png')][:10],
                        caption='KEGG Level2水平分布图',
                        description='图片说明：横坐标是注释到各Level2代谢通路的上调（下调）差异代谢物和所有注释到KEGG通路的上调（下调）差异表达代谢物总数的比值（%），纵轴表示Level2 pathway的名称，柱子右边数字代表注释到该Level2 pathway下的上调（下调）差异代谢物数量。')

                # kegg网络通路图
                link = list(glob(f'{endo_dir}/*通路富集/*KEGG_map.zip'))
                diff_ana11 = diff_ana5.add_section("KEGG网络通路图", link=mv_top(link[0]))
                # if len(list(glob(f'{endo_dir}/*通路富集/KEGG_map/*/*.png'))) == 0:
                diff_ana11_plot = diff_ana5.add_plot('pathway.png',
                                                     caption='代谢通路图')

            # Reactome
            if list(glob(f'{endo_dir}/*通路富集/Reactome')):
                link = list(glob(f'{endo_dir}/*通路富集/Reactome'))
                Reactome_ana5 = diff_ana_Term.add_section("Reactome富集分析", link=mv_top(link[0]))

            if list(glob(f'{endo_dir}/*通路富集/*Reactome/*/*Total.xls')):
                listdiff = list(glob(f'{endo_dir}/*通路富集/*Reactome/*/*Total.xls'))
                Reactome_ana5_1 = Reactome_ana5.add_table(listdiff[0], show_rows=10, description='''表头解释如下：

[1]. id：Reactome 通路编号

[2]. Term：该通路的描述

[3]. ListHits：该通路中总差异代谢物数

[4]. ListTotal：注释到Reactome的总差异代谢物数

[5]. PopHits：注释到该通路中的所有代谢物个数

[6]. PopTotal：注释到Reactome的总代谢物数

[7]. p-value：富集显著性pvalue值，P≤0.05表示显著富集

[8]. padjust：校正后的pvalue值

[9]. Enrichment_score：富集打分

[10]. Substances： 属于该条通路的代谢物
        ''')

                ##Reactome_top20_bar
                Reactome_ana5_1 = Reactome_ana5.add_section('Reactome富集分析top20柱状图', link=mv_top(link[0]))
                if len(list(glob(f'{endo_dir}/*通路富集/Reactome/*/Reactome.Bar*.png'))) == 0:
                    Reactome_ana5_plot = Reactome_ana5_1.add_plot('Reactome Level3-example.png',
                                                                  caption='Reactome Level3',
                                                                  description='图片说明：横坐标为每条通路的-log10 p值，纵坐标为不同通路名称，柱子上的数字为注释到该通路的差异代谢物数量。')
                else:
                    Reactome_ana5_plot = Reactome_ana5_1.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*通路富集/Reactome/*/Reactome.Bar*.png')][:10],
                        caption='Reactome Level3',
                        description='图片说明：横坐标为每条通路的-log10 p值，纵坐标为不同通路名称，柱子上的数字为注释到该通路的差异代谢物数量。')
                ##Reactome_top20_point
                Reactome_ana6 = Reactome_ana5.add_section("Reactome富集分析top20气泡图", link=mv_top(link[0]))
                if len(list(glob(f'{endo_dir}/*通路富集/Reactome/*/Reactome.Bubble*.png'))) == 0:
                    Reactome_ana6_plot = Reactome_ana6.add_plot('Reactome.Bubble-example.png',
                                                                caption='Reactome气泡图',
                                                                description='图片说明：图中横坐标Enrichment Score为富集分值，纵坐标为top20的通路信息。气泡越大的通路包含的差异代谢物越多，气泡颜色由蓝-红变化，其富集pvalue值越小，显著程度越大。')
                else:
                    Reactome_ana6_plot = Reactome_ana6.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*通路富集/Reactome/*/Reactome.Bubble*.png')][:10],
                        caption='Reactome气泡图',
                        description='图片说明：图中横坐标Enrichment Score为富集分值，纵坐标为top20的通路信息。气泡越大的通路包含的差异代谢物越多，气泡颜色由蓝-红变化，其富集pvalue值越小，显著程度越大')

                ##Reactome_chord
                Reactome_ana7 = Reactome_ana5.add_section("Reactome富集分析和弦图", link=mv_top(link[0]))
                if len(list(glob(f'{endo_dir}/*通路富集/KEGG/*/KEGG.chord.png'))) == 0:
                    Reactome_ana7_plot = Reactome_ana7.add_plot('Reactome.chord-example.png',
                                                                caption='Reactome富集分析和弦图',
                                                                description='图片说明：左侧为代谢KEGG C number，红色表示上调，蓝色表示下调，右侧为所选Reactome通路。')
                else:
                    Reactome_ana7_plot = Reactome_ana7.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*通路富集/Reactome/*/Reactome.chord.png')][:10],
                        caption='Reactome富集分析和弦图',
                        description='图片说明：左面为代谢KEGG C number，红色表示上调，蓝色表示下调，右面为所选Reactome通路。')

                ##Reactome_up_down
                Reactome_ana9 = Reactome_ana5.add_section("Reactome富集分析上下调通路对比图", link=mv_top(link[0]))
                if len(list(glob(f'{endo_dir}/*通路富集/Reactome/*/Up_Down.Comparison.png'))) == 0:
                    Reactome_ana9_plot = Reactome_ana9.add_plot('Up_Down.Comparison-example.png',
                                                                caption='Reactome富集分析上下调通路对比图',
                                                                description='图片说明：横坐标为该通路前景差异代谢物数量与背景代谢物数量之比（ListHits/TotalHits），纵坐标为Reactome通路名称，有的通路会在上下调均显著富集。')
                else:
                    Reactome_ana9_plot = Reactome_ana9.add_plot(
                        [(f, '') for f in glob(f'{endo_dir}/*通路富集/Reactome/*/Up_Down.Comparison.png')][:10],
                        caption='Reactome富集分析上下调通路对比图',
                        description='图片说明：横坐标为该通路前景差异代谢物数量与背景代谢物数量之比（ListHits/TotalHits），纵坐标为Reactome通路名称，有的通路会在上下调均显著富集。')

                # GSEA
                if list(glob(f'{endo_dir}/*通路富集/GSEA/*/*.gsea.png')):
                    link = list(glob(f'{endo_dir}/*通路富集/*GSEA'))
                    diff_ana3_4 = diff_ana_Term.add_section("GSEA分析", link=mv_top(link[0]))
                    diff_ana3_4_1 = diff_ana3_4.add_section("GSEA富集分析", link=mv_top(link[0]))
                    if len(list(glob(f'{endo_dir}/*通路富集/GSEA/*/*.gsea.png'))) == 0:
                        diff_ana3_4_1_plot = diff_ana3_4_1.add_plot('GSEA-example.png',
                                                                    caption='GSEA')
                    else:
                        diff_ana3_4_1_plot = diff_ana3_4_1.add_plot(
                            [(f, '') for f in glob(f'{endo_dir}/*通路富集/GSEA/*/*gsea.png')][:10],
                            caption='GSEA富集分析',
                            description='图片说明：GSEA分析的图形对应一个数据集的分析结果，图形主要分为4部分，从上至下依次为：① 富集分数（enrichment score， ES）的分布图，绿线为所有差异代谢物的ES分布情况，该曲线在Y轴绝对值最大的位置对应该数据集的富集分数,当ES>0 峰值左侧为核心差异代谢物，ES<0 则峰值右侧为核心差异代谢物；② 数据集代谢物分布图，竖线表示该数据集中代谢物在整个排序中的位置；③ Colorbar，即排序矩阵的颜色映射，值为正值，对应红色，值越大越红，反之，对应蓝色。越趋近与0，越接近白色。④ 排序矩阵分布图，比如差异倍数，signal2noise等数字的分布情况.')

                    diff_ana3_4_2 = diff_ana3_4.add_section("GSEA聚类图", link=mv_top(link[0]))
                    if len(list(glob(f'{endo_dir}/*通路富集/GSEA/*/*.gsea.png'))) == 0:
                        diff_ana3_4_2_plot = diff_ana3_4_2.add_plot('GSEA-heatmap-example.png', caption='GSEA聚类图示例')
                    else:
                        diff_ana3_4_2_plot2 = diff_ana3_4_2.add_plot(
                            [(f, '') for f in glob(f'{endo_dir}/*通路富集/GSEA/*/*heatmap.png')][:10],
                            caption='GSEA聚类图',
                            description='图片说明：图中红色表示相对高表达代谢物，蓝色表示相对低表达代谢物')
        # 入血代谢物与差异代谢物相关性分析
        if list(glob(f'{endo_dir}/*入血代谢物与差异代谢物相关性分析/')):
            rel_diff_blood = Analysis_module.add_section('入血代谢物与差异代谢物相关性分析', link=mv_top(list(glob(f'{endo_dir}/*入血代谢物与差异代谢物相关性分析/'))[0]))
            rel_diff_blood.add_table(list(glob(f'{endo_dir}/*入血代谢物与差异代谢物相关性分析/*.xlsx'))[0],
                                     caption="入血代谢物与差异代谢物相关性数据表",
                                     description="""表头描述如下：
                                     
 [1] Endogenous: 内源代谢物；
 
 [2] Ingredient: 入血成分；
 
 [3] cor: 相关性R值，绝对值越接近1相关性越高；
 
 [4] p: 显著性P值，p<0.05为显著，<0.01为极显著；
 
 [5] class: 相关性类别，正相关或负相关。
                                     
                                     """, show_rows=10, show_search=False)
            rel_diff_blood.add_plot(list(glob(f'{endo_dir}/*入血代谢物与差异代谢物相关性分析/*.jpg'))[0],
                                    caption="入血代谢物与差异代谢物相关性分析网络图",
                                    description='图片说明：方形代表入血代谢物，圆形代表差异代谢物，节点的大小由degree决定。红色连接线代表正相关，蓝色为负相关。相关性节点的筛选标准: P值<=0.01且R2>=0.99。')
        diff_ana8 = report.add_section("参考文献")
    # ---------------------内源代谢物分析 end---------------
    # FAQ & 推广部分
    diff_ana9 = report.add_section("常见问题")
    diff_ana10 = report.add_section("高级分析 & 个性化分析内容展示")
    diff_ana10_1 = diff_ana10.add_section('STEM（时间趋势/浓度趋势）分析')

    diff_ana10_1_plot = diff_ana10_1.add_plot("时间序列.jpg")

    diff_ana10_2 = diff_ana10.add_section('WGCNA分析')
    diff_ana10_2_plot = diff_ana10_2.add_plot("图片4.png")

    diff_ana10_2_1 = diff_ana10_2.add_section('网络构建')
    diff_ana10_2_1_plot = diff_ana10_2_1.add_plot("图片5.png", caption='WGCNA网络构建参数')

    diff_ana10_2_2 = diff_ana10_2.add_section('模块识别及分析')
    diff_ana10_2_2_plot1 = diff_ana10_2_2.add_plot("图片6.png", caption='模块变量聚类数图')
    diff_ana10_2_2_plot2 = diff_ana10_2_2.add_plot("图片7.png", caption='变量TOM聚类热图')

    diff_ana10_2_3 = diff_ana10_2.add_section('模块性状关联分析')
    diff_ana10_2_3_plot = diff_ana10_2_3.add_plot("图片8.png", caption=' 性状模块关联热图')

    diff_ana10_2_4 = diff_ana10_2.add_section('核心变量分析')
    diff_ana10_2_4_plot = diff_ana10_2_4.add_plot("图片9.png",
                                                  description='图片说明：网络图中的线表示变量之间的联系程度，变量与周边点的联系越多，在网络中越处于核心地位')

    diff_ana10_3 = diff_ana10.add_section('大队列样本biomarker筛选')
    diff_ana10_3_1 = diff_ana10_3.add_section('集成机器学习联合策略筛选Biomakers')
    diff_ana10_3_1_plot1 = diff_ana10_3_1.add_plot("图片10.png", caption='集成机器学习分析流程图')
    diff_ana10_3_1_plot2 = diff_ana10_3_1.add_plot("图片11.png", caption='集成机器学习分析结果部分展示')

    diff_ana10_3_2 = diff_ana10_3.add_section('biomarker筛选方法2')

    diff_ana10_3_2_plot = diff_ana10_3_2.add_plot("biomarker2-流程.png", caption='biomarker筛选')
    diff_ana10_3_2_plot = diff_ana10_3_2.add_plot("biomarker2.png", caption='biomarker筛选')

    diff_ana10_4 = diff_ana10.add_section('高维中介分析')
    diff_ana10_4_plot1 = diff_ana10_4.add_plot("图片16.png", caption='高维中介路径')
    diff_ana10_4_plot2 = diff_ana10_4.add_plot("图片17.png", caption='组间差异分析及箱线图分析')

    diff_ana10_5 = diff_ana10.add_section('多组学整合/单组学亚型分析')

    diff_ana10_5_plot = diff_ana10_5.add_plot("多组学整合.png")

    diff_ana10_7 = diff_ana10.add_section('亚细胞定位')
    diff_ana10_7_plot = diff_ana10_7.add_plot("亚细胞定位.png")

    diff_ana10_6 = diff_ana10.add_section('其它个性化分析绘图')

    diff_ana10_6_1 = diff_ana10_6.add_section('分类差异火山图')
    diff_ana10_6_1_plot = diff_ana10_6_1.add_plot("图片22.png", caption='分类差异火山图')

    diff_ana10_6_2 = diff_ana10_6.add_section('分类气泡图')
    diff_ana10_6_2_plot = diff_ana10_6_2.add_plot("图片23.png", caption='分类气泡图',
                                                  description='图片说明:每个气泡代表一个差异脂质代谢物，横坐标是脂质的Sub Class分类，同Class分类下的SubClas分类的脂质代谢物气泡颜色映射同色系，纵坐标是比较组的log2(FC)值，圆点越大，表示P-value值越大。')

    diff_ana10_6_3 = diff_ana10_6.add_section('相关性和弦图')
    diff_ana10_6_3_plot = diff_ana10_6_3.add_plot("图片24.png", caption='相关性和弦图')

    diff_ana10_6_4 = diff_ana10_6.add_section('桑基图')

    diff_ana10_6_4_plot = diff_ana10_6_4.add_plot("图片25.png", caption='桑基图')
    # FAQ & 推广部分

    # 公司简介
    reportinfo.companyinfo(report=report, yamlpath=os.environ['OEBIO'])

    proj_id = config["header_info"]["项目编号"]
    # export
    report.write_to(f"{report_dir}/{proj_id}_report.html", zip_report_name=report_dir+f".zip")

    return report_dir, proj_id
    # 添加与报告到 /public/项目报告检查空间/预报告文件夹/
    # pnum = Path(data).name  # project_info.loc[0, "项目编号"]
    # pname = project_info.loc[0, "项目名称"]
    # # shutil.copy(f"{report_dir}/项目报告.html", f"/public/项目报告检查空间/预报告文件夹/{pnum}-{pname}-LCMS-1.0.0-分析预报告.html")
    # shutil.copy(f"{report_dir}/项目报告.html", f"/data/nas/177/代谢/项目报告检查空间/预报告文件夹/{pnum}-{pname}-LCMS-1.0.0-分析预报告.html")


@add_loginfo
def oecloudProduct(data, version="/data/hstore1/database/report/2023-03-04~2023/src", company="LM",
                 using_product=False, using_herb=False, cloud=True, local=0):

    from lmbio.basic.reportinfo import getreportinfo
    os.environ['OEBIO'] = version
    reportinfo = getreportinfo(sys_logo=company)
    if cloud:
        os.environ['CLOUD'] = "True"
        os.environ['CLOUDUSER'] = "lmbioinformatics@lumingbio.com"
        os.environ['PASSWORD'] = "Lmbio@123"

    src_path = Path("/data/hstore1/database/database/tcm/")

    project_name, project_info = parse_project(data, local)
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
        shutil.copy(str(xls), f"{str(Path(data) / project_name / '实验内容')}/{xls.name}")

    # 拷贝UV data
    uv_path = [path for path in (Path(data) / "UVdata").rglob("*紫外吸收图*")]
    for uv in uv_path:
        shutil.copy(str(uv), f"{str(Path(data) / project_name / '数据分析' / '1.紫外吸收图')}/{uv.name}")

    # 设置项目信息
    config["header_info"]["项目名称"] = '中药入血成分分析'
    # config["header_info"]["客户单位"] = project_info.loc[0, "客户单位"]
    task_id = None
    try:
        task_id = project_info.loc[0, "任务单号"]
    except:
        task_id = project_info.loc[0, "项目编号"]
    config["header_info"]["任务单号"] = task_id if "-b" in task_id else task_id + "-b1"
    config["header_info"]["客户名称"] = project_info.loc[0, "客户名称"] if project_info.loc[0, "客户名称"] != '' and pd.notna(project_info.loc[0, "客户名称"]) else project_info.loc[0, "联系人名称"]
    config["header_info"]["联系人名称"] = project_info.loc[0, "联系人名称"] if project_info.loc[0, "联系人名称"] != '' and pd.notna(project_info.loc[0, "联系人名称"]) else project_info.loc[0, "客户名称"]
    if config["header_info"]["客户名称"] == config["header_info"]["联系人名称"]:
        config["header_info"].pop("客户名称")
    config["header_info"]["项目编号"] = Path(data).cwd().name
    config["header_info"]["样本"] = project_info.loc[0, "样本"]
    # config["header_info"]["完成时间"] = datetime.datetime.now().strftime("%Y-%m-%d")
    # 默认 jiang tao
    # config["header_info"]["执行编码"] = 'LM0460'

    header_info = config.get("header_info")
    # OUT PATH
    report_dir = str(Path(data) / project_name)

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

    proj_id = config["header_info"]["项目编号"]
    # export
    report.write_to(f"{report_dir}/{proj_id}_report.html", zip_report_name=report_dir+f".zip")

    return report_dir, proj_id

    # # 添加与报告到 /public/项目报告检查空间/预报告文件夹/
    # pnum = Path(data).name # project_info.loc[0, "项目编号"]
    # pname = project_info.loc[0, "项目名称"]
    # # shutil.copy(f"{report_dir}/项目报告.html", f"/public/项目报告检查空间/预报告文件夹/{pnum}-{pname}-LCMS-1.0.0-分析预报告.html")
    # shutil.copy(f"{report_dir}/项目报告.html", f"/data/nas/177/代谢/项目报告检查空间/预报告文件夹/{pnum}-{pname}-LCMS-1.0.0-分析预报告.html")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser("oebio report")

    parser.add_argument("-d","--data", type=str, default=".", help="文件夹路径")
    parser.add_argument("-v","--version", type=str, default="/data/hstore1/database/report/2023-03-04~2023/src", help="version")
    parser.add_argument("-s","--company", type=str, default="LM", help="company")
    parser.add_argument("-p","--product", default=False, action='store_true', help="是否分析代谢产物")
    parser.add_argument("-b","--herb", default=False, action='store_true', help="是否分析代谢HERB")
    parser.add_argument("-c","--cloud", default=False, action='store_true', help="是否上传云平台")
    parser.add_argument("-t", "--type", required=True, choices=("tcm", "blood", "product"), help="报告类型")
    parser.add_argument("-l", "--local", type=int, required=True, default=0, help="是否使用本地登记单，默认不使用")
    parser.add_argument("-k", "--skip", type=int, default=0, help="是否跳过 ~do before report~，默认不跳过")

    args = parser.parse_args()
    if args.type == "tcm":
        oecloudTcm(args.data,
                     using_herb=args.herb,
                     cloud=args.cloud,
                     local=args.local)
    elif args.type == "blood":
        oecloudBlood(args.data,
                     using_product=args.product,
                     using_herb=args.herb,
                     cloud=args.cloud,
                     local=args.local,
                     skip=True if args.skip else False)
    elif args.type == "product":
        oecloudProduct(args.data,
                     using_product=args.product,
                     using_herb=args.herb,
                     cloud=args.cloud,
                     local=args.local)