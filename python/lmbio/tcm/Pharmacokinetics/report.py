# File: report
# Author: jiang tao
# Time: 2023/6/27 10:19
import os
import yaml
import argparse
import warnings

from pathlib import Path
from oebio.report import oeweb_register, Report
from oebio.app import *
from lmbio.basic.reportinfo import getreportinfo
from statistics import *

if __name__ == '__main__':
    warnings.filterwarnings(action='ignore')
    # 解析参数
    parser = argparse.ArgumentParser(description = "药代动力学分析流程")
    parser.add_argument('-m', "--mode", default="LM", help='模板')
    parser.add_argument('-p', "--path", default=".", help='项目文件夹的路径')
    args = parser.parse_args()
    # 工作目录
    wkdir = Path(args.path)
    if wkdir.exists():
        os.chdir(wkdir)
    else:
        wkdir.mkdir(parents=True, exist_ok=True)
        os.chdir(wkdir)
    # ------------------------------药代动力学分析------------------------------
    used_cols = ["Sample Name", "Component Name", "Calculated Concentration"]
    # save dir
    f = str(wkdir)
    # input data
    p = wkdir / "data.xlsx"

    nc = pd.read_excel(p, sheet_name='nominal_concentration')
    df1 = pd.read_excel(p, sheet_name='精密度及准确度')[used_cols]
    get_precision_accuracy(df1, nc, f)

    df2 = pd.read_excel(p, sheet_name='提取回收率及基质效应')[used_cols]
    get_recovery_matrix_effect(df2, nc, f)

    df3 = pd.read_excel(p, sheet_name='稳定性')[used_cols]
    get_stability(df3, nc, f)

    df4 = pd.read_excel(p, sheet_name='药代动力学数据')[used_cols]
    get_analyte_cst(df4, f)

    df5 = pd.read_excel(p, sheet_name='药代动力学数据说明')[["样本编号", "时间点"]]
    curve_data = get_curve_data(df5, df4)
    plot_drug_time_curve(curve_data, f)

    desc = pd.read_excel(p, sheet_name="药代动力学数据说明")
    demo = pd.read_excel(p, sheet_name="药代动力学数据")[used_cols]
    calc_kinetics_params(desc, demo, f)
    # ------------------------------药代动力学分析 end---------------------------


    # --------------------------------OEBIO 报告 --------------------------------
    #添加环境变量
    os.environ['OEBIO'] = "/data/hstore1/database/report/2023-03-04~2023/src"

    reportinfo = getreportinfo(sys_logo=args.mode)

    # TODO 动态获取项目信息
    header_info = {
        '项目名称': "中药多组分药代动力学分析",
        '客户单位': "**",
        '客户姓名': "**",
        '样本': "**",
        '项目编号': 'DZLM202304001-T',
        '完成时间': '2023年4月'
    }
    report = Report('中药多组分药代动力学分析报告', title='中药多组分药代动力学分析报告', header_info=header_info, oe_welcome=reportinfo.oe_weclome)

    p = Path("/data/hstore1/database/database/tcm/pharmacokinetics")
    report.add_yaml_config(str(p / 'report.yaml'))
    config = yaml.load(open(str(p / "config.yaml"), encoding="utf8"), Loader=yaml.FullLoader)
    # 1 前言
    intro = report.add_section('前言')

    # 2 材料与方法
    mater_method = report.add_section('材料与方法')
    # 2.1 样本信息
    sample_info = mater_method.add_section('样本信息')
    sample_info.add_table(config['样本信息']['table_path'], caption = config['样本信息']['table_name'], show_rows=11)
    # 2.2 标准品与质控品
    sec22 = mater_method.add_section('标准品与质控品')
    sec22.add_table(config['标准品与质控品']['table_path'], caption = config['标准品与质控品']['table_name'])
    # 2.3 仪器设备
    sec23 = mater_method.add_section('仪器设备')
    sec23.add_table(config['仪器设备']['table_path'], caption = config['仪器设备']['table_name'])
    sec23.add_table(config['仪器设备']['table_path1'], caption = config['仪器设备']['table_name1'])
    # 2.4 LCMS分析方法
    sec24 = mater_method.add_section('LCMS分析方法')
    sec24.add_table(config['LCMS分析方法']['table_path'], caption = config['LCMS分析方法']['table_name'])
    sec24.add_comment(description=config['LCMS分析方法']['table1_desc'])
    sec24.add_table(config['LCMS分析方法']['table_path1'], caption = config['LCMS分析方法']['table_name1'])
    # 2.5 样品前处理
    sec25 = mater_method.add_section('样品前处理')

    # 2.6 方法学验证
    sec26 = mater_method.add_section('方法学验证')

    # 2.6.1 选择性
    sec261 = sec26.add_section('选择性', description=config['选择性']['description'])

    # 2.6.2 线性及灵敏度
    sec262 = sec26.add_section('线性及灵敏度', description=config['线性及灵敏度']['description'])
    # 2.6.3 精密度及准确度
    sec263 = sec26.add_section('精密度及准确度', description=config['精密度及准确度']['description'])
    # 2.6.4 提取回收率及基质效应
    sec264 = sec26.add_section('提取回收率及基质效应')
    # 2.6.5 稳定性
    sec265 = sec26.add_section('稳定性', description=config['稳定性']['description'])

    # 3  数据分析
    sec3 = report.add_section('数据分析')

    # 4 结果与讨论
    sec4 = report.add_section('结果与讨论')
    # 4.1 标准品色谱图及质谱图
    # TODO: GET FROM DIR
    sec41 = sec4.add_section('标准品色谱图及质谱图')
    sec41.add_fig(config['标准品色谱图及质谱图']['fig_path'],caption = config['标准品色谱图及质谱图']['fig_name'])
    sec41.add_fig(config['标准品色谱图及质谱图']['fig_path1'],caption = config['标准品色谱图及质谱图']['fig_name1'])
    sec41.add_fig(config['标准品色谱图及质谱图']['fig_path2'],caption = config['标准品色谱图及质谱图']['fig_name2'])
    # 4.2 方法学验证
    sec42 = sec4.add_section('方法学验证')
    # 4.2.1 选择性
    sec421 = sec42.add_section('选择性')
    sec421.add_fig(config['选择性']['fig_path'],caption = config['选择性']['fig_name'])
    # 4.2.2 线性及灵敏度
    sec422 = sec42.add_section('线性及灵敏度')
    sec422.add_table(config['线性及灵敏度']['table_path'], caption = config['线性及灵敏度']['table_name'])
    # 4.2.3 精密度及准确度
    sec423 = sec42.add_section('精密度及准确度')
    sec423.add_table(config['精密度及准确度']['table_path'], caption = config['精密度及准确度']['table_name'])
    # 4.2.4 回收率及基质效应
    sec424 = sec42.add_section('回收率及基质效应')
    sec424.add_table(config['回收率及基质效应']['table_path'], caption = config['回收率及基质效应']['table_name'])
    # 4.2.5 稳定性
    sec425 = sec42.add_section('稳定性')
    sec425.add_table(config['稳定性']['table_path'], caption = config['稳定性']['table_name'])
    # 4.3 中药多组分含量
    sec43 = sec4.add_section('中药多组分含量')
    sec43.add_table(config['中药多组分含量']['table_path'], caption = config['中药多组分含量']['table_name'])
    # 4.4 药代动力学分析
    sec44 = sec4.add_section('药代动力学分析')
    # 4.4.1 药时曲线
    sec441 = sec44.add_section('药时曲线')
    sec441.add_fig(config['药时曲线']['fig_path'], caption = config['药时曲线']['fig_name'])
    # 4.4.2 药代动力学参数
    sec442 = sec44.add_section('药代动力学参数')
    sec442.add_table(config['药代动力学参数']['table_path'], caption = config['药代动力学参数']['table_name'])
    # 5 参考文献
    sec5 = report.add_section('参考文献')

    #公司简介
    reportinfo.companyinfo(report=report,
                           yamlpath=os.environ['OEBIO'])

    report.write_to("中药多组分药代动力学分析报告.html")
    # --------------------------------OEBIO end --------------------------------




