#!/opt/conda/bin/python
# File: commands
# Author: jiang tao
# Time: 2023/4/10 16:24
import click
import sys


# filter error info: font family not found!
import pandas as pd


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


@click.group()
def cli():
    """\b
        中药成分鉴定和入血分析命令行工具
        命令组成：
            1.`tcm`命令 中药成分鉴定分析流程
            2.`blood`命令 中药入血成分析流程
            3.`blood_check`命令  人工检查删除EIC效果不好的成分
        售后：
            1.`tcm_feedback`命令  中药成分鉴定售后
                a)添加MS2: `lmtcm tcm_feedback -I /data/hstore4/lumingos/project/DZLM2023050885 -t a`
            2.`blood_feedback`命令  中药入血成分析售后
                a)添加MS2: `lmtcm blood_feedback -I /data/hstore4/lumingos/project/DZLM2023050885 -t a`
                b)提供所有候选结果，拓展代谢产物矩阵: `lmtcm blood_feedback -I /data/hstore4/lumingos/project/DZLM2023050885 -t b`
        """
    ...


@cli.command('tcm', short_help="中药成分鉴定分析流程")
@click.option("-I", "--input", 'input_folder', type=click.Path(), required=True, default=".", help="输入文件的文件夹路径，路径中不能包含中文")
@click.option("-B", "--herb", 'using_herb', default=False, is_flag=True, help="是否分析HERB库,布尔值标记[加上参数即开启]")
@click.option("-M", "--water_mask", 'water_mask', default=False, is_flag=True, help="是否添加水印,布尔值标记[加上参数即开启]")
@click.option("-co", "--cloud_off", 'cloud_off', default=True, is_flag=True, help="是否云交付,布尔值标记[默认开启，加上参数关闭]")
@click.option("-oe", "--oereport_only", 'oereport_only', default=False, is_flag=True, help="只执行函数oebio_report")
@click.option("-local", "--local", 'local', default=False, is_flag=True, help="执行本地项目需要加上此参数")
def tcm(input_folder, using_herb=False, water_mask=False, cloud_off=True, oereport_only=False, local=False):
    from TCMreport import Report
    r = Report(data=input_folder, local=local)
    if oereport_only:
        r.oebio_report(using_herb=True, cloud=cloud_off)
    else:
        r.run_pipeline(using_herb=True, water_mask=water_mask, cloud=cloud_off)

    # interested_cpd
    # r.preprocess_qi(using_herb=using_herb, interested_cpd=["Matrine","Ammothamnine"])
    # r.plot_xic(interested_cpd=["Matrine","Ammothamnine"])

        # r.preprocess_qi(using_herb=using_herb)
        # r.plot_chromatograms()
        # r.plot_xic(dpi=600, img_size=250)
        # r.plot_pies()
        # r.oebio_report(using_herb=using_herb, cloud=cloud_off)
        # r.zip_report()
        # r.test()


# 复方中是否包含动物
@cli.command('blood', short_help="中药入血成分析流程")
@click.option("-I", "--input", 'input_folder', type=click.Path(), default=".", required=True, help="输入文件的文件夹路径，路径中不能包含中文")
@click.option("-B", "--herb", 'using_herb', default=False, is_flag=True, help="是否分析HERB库,布尔值标记[加上参数即开启]")
@click.option("-p", "--product", 'using_product', default=False, is_flag=True, help="是否分析代谢产物,布尔值标记[加上参数即开启]")
@click.option("-fc", "--fc", 'FC', default=10, type=float, help="判别入血成分的差异倍数")
@click.option("-co", "--cloud_off", 'cloud_off', default=True, is_flag=True, help="是否云交付,布尔值标记[默认开启，加上参数关闭]")
@click.option("-oe", "--oereport_only", 'oereport_only', default=False, is_flag=True, help="只执行函数oebio_report")
@click.option("-local", "--local", 'local', default=False, is_flag=True, help="执行本地项目需要加上此参数")
def blood(input_folder, using_herb=False, using_product=False, FC=10, cloud_off=True, oereport_only=False, local=False):
    from TCMreport import BloodReport
    br = BloodReport(data=input_folder, project="Blood", local=local)
    if oereport_only:
        br.oebio_report(using_herb=True, using_product=using_product, cloud=cloud_off)
    else:
        br.run_pipeline(using_herb=True, using_product=using_product, FC=FC)

        # # interested_cpd
        # br.preprocess_qi(using_herb=using_herb, using_product=using_product, FC=FC, interested_cpd=["Matrine", "Ammothamnine"])
        # br.plot_xic(interested_cpd=["Matrine", "Ammothamnine"])
        # br.rename_excel(interested_cpd=["Matrine", "Ammothamnine"])

        # br.preprocess_qi(using_herb=using_herb, using_product=using_product, FC=FC)
        # br.plot_chromatograms()
        # if using_product:
        #     br.network_analysis()
    # br.plot_xic()
        # br.rename_excel()
        # br.plot_pies()
        # br.oebio_report(using_herb=using_herb, using_product=using_product)
        # br.zip_report()

        # br.test()
        # br.split_df_only(using_product=using_product)
        # br.add_src_in_prescript()
        # br.plot_tic(absolute_intensity=True)
        # br.plot_tic_split(absolute_intensity=True, add_peak_idx=True)


@cli.command('blood_check', short_help="中药入血成分析流程 - 删除EIC效果不好的入血成分")
@click.option("-I", "--input", 'input_folder', type=click.Path(), default=".", required=True, help="输入文件的文件夹路径，路径中不能包含中文")
@click.option("-S", "--skip", 'skip', default=True, is_flag=True,
              help="是否删除手动检查的入血成分, 默认检查，需填写delCpds.xlsx, 布尔值标记[加上参数--skip跳过检查]")
@click.option("-B", "--herb", 'using_herb', default=False, is_flag=True, help="是否分析HERB库,布尔值标记[加上参数即开启]")
@click.option("-p", "--product", 'using_product', default=True, is_flag=True, help="是否分析代谢产物,布尔值标记[加上参数即开启]")
@click.option("-r", "--raw", 'del_blood', default=True, is_flag=True, help="是否是删除入血代谢产物,布尔值标记[加上参数删除原方cpd]")
@click.option("-co", "--cloud_off", 'cloud_off', default=True, is_flag=True, help="是否云交付,布尔值标记[默认开启，加上参数关闭]")
@click.option("-local", "--local", 'local', default=False, is_flag=True, help="执行本地项目需要加上此参数")
def blood_check(input_folder, skip=False, using_herb=False, using_product=False, del_blood=True, cloud_off=True, local=False):
    from TCMreport import BloodReport
    br = BloodReport(data=input_folder, project="Blood", local=local)
    br.manual_inspect(skip=skip, using_herb=True, using_product=using_product, del_blood=del_blood,
                      cloud=cloud_off)


@cli.command('blood_feedback', short_help="""中药入血成分析流程 - 售后:
                                          a) 扩展代谢网络数据矩阵
                                          b) 提供所有二级""")
@click.option("-I", "--input", 'input_folder', type=click.Path(), default=".", required=True, help="输入文件的文件夹路径，路径中不能包含中文")
@click.option("-t", "--type", 'type', type=click.Choice(['a', 'b']), required=True, help="a)添加MS2 b)提供所有候选结果，拓展代谢产物矩阵 ")
@click.option("-local", "--local", 'local', default=False, is_flag=True, help="执行本地项目需要加上此参数")
def blood_feedback(input_folder, type, local=False):
    from TCMreport import BloodReport
    br = BloodReport(data=input_folder, project="Blood", local=local)
    if type == 'a':
        br.extend_network_df()
    elif type == 'b':
        br.feedback_all_detected_MS2()
    else:
        raise NotImplementedError("未知售后类型")


@cli.command('tcm_feedback', short_help="""中药成分鉴定分析流程 - 售后:
                                          a) 提供所有二级
                                          b) ...""")
@click.option("-I", "--input", 'input_folder', type=click.Path(), default=".", required=True, help="输入文件的文件夹路径，路径中不能包含中文")
@click.option("-t", "--type", 'type', type=click.Choice(['a', 'b']), required=True, help="a)添加MS2")
@click.option("-local", "--local", 'local', default=False, is_flag=True, help="执行本地项目需要加上此参数")
def blood_feedback(input_folder, type, local=False):
    from TCMreport import Report
    r = Report(data=input_folder, local=local)
    if type == 'a':
        r.feedback_all_detected_MS2()
    else:
        raise NotImplementedError("未知售后类型")



@cli.command('product_resolve', short_help="""代谢产物解析""")
@click.option("-I", "--input", 'input_folder', type=click.Path(), default=".", required=True, help="输入文件的文件夹路径，路径中不能包含中文")
@click.option("-B", "--herb", 'using_herb', default=False, is_flag=True, help="是否分析HERB库,布尔值标记[加上参数即开启]")
@click.option("-p", "--product", 'using_product', default=True, is_flag=True, help="是否分析代谢产物,布尔值标记[加上参数即关闭]")
@click.option("-co", "--cloud_off", 'cloud_off', default=True, is_flag=True, help="是否云交付,布尔值标记[默认开启，加上参数关闭]")
@click.option("-local", "--local", 'local', default=False, is_flag=True, help="执行本地项目需要加上此参数")
def blood_feedback(input_folder, using_herb, using_product, cloud_off=True, local=False):
    from TCMreport import ProductResolve
    r = ProductResolve(data=input_folder, project="Blood", local=local)
    r.raw2mzml()
    if not local:
        r.raw2obs()
    # r.preprocess_qi(using_herb=True, using_product=using_product)
    # r.adjustRTByDetected()
    # r.network_analysis()
    # r.plot_xic()
    # r.rename_excel()
    # r.plot_chromatograms()
    # r.oebio_report(using_herb=True, using_product=using_product, cloud=cloud_off)
    # r.zip_report()


@cli.command('after_sales', short_help="""售后模块""")
@click.option("-I", "--input", 'input_folder', type=click.Path(), required=True, default=".",
              help="输入文件的文件夹路径，路径中不能包含中文")
@click.option("-B", "--herb", 'using_herb', default=False, is_flag=True, help="是否分析HERB库,布尔值标记[加上参数即开启]")
@click.option("-p", "--product", 'using_product', default=True, is_flag=True, help="是否分析代谢产物,布尔值标记[加上参数即关闭]")
@click.option("-co", "--cloud_off", 'cloud_off', default=True, is_flag=True, help="是否云交付,布尔值标记[默认开启，加上参数关闭]")
@click.option("-oe", "--oereport_only", 'oereport_only', default=False, is_flag=True, help="只执行函数oebio_report")
@click.option("-local", "--local", 'local', default=False, is_flag=True, help="执行本地项目需要加上此参数")
@click.option("-t", "--type", 'type', type=click.Choice(['TCM', 'Blood']), required=True, help="项目类型： TCM OR Blood")
@click.option("-x", "--xlsx", 'xlsx', type=click.Path(), help="interested_cpd列表： 感兴趣的化合物列表[Metabolites, ID]")
@click.option("-fc", "--fc", 'FC', default=10, type=float, help="判别入血成分的差异倍数")
def dev(input_folder, FC=10, using_herb=True, using_product=False, water_mask=False, cloud_off=True, oereport_only=False, local=False, type="TCM", xlsx=None):
    """lmtcm after_sales -t Blood -x cpds.xlsx --local"""
    peak_id = None
    df = pd.read_excel(xlsx, sheet_name=0)
    interested_cpd = df["Metabolites"].unique().tolist()
    if "ID" in df.columns:
        peak_id = df["ID"].unique().tolist()
    if type == "TCM":
        from TCMreport import Report
        r = Report(data=input_folder, local=local)
        r.preprocess_qi(using_herb=True, interested_cpd=interested_cpd, peak_id=peak_id)
    else:
        from TCMreport import BloodReport
        r = BloodReport(data=input_folder, project="Blood", local=local)
        r.preprocess_qi(FC=FC, using_herb=True, using_product=using_product, interested_cpd=interested_cpd, peak_id=peak_id)
    r.plot_xic(interested_cpd=[1, 2, ])



@cli.command('dev', short_help="""开发测试模块""")
@click.option("-I", "--input", 'input_folder', type=click.Path(), required=True, default=".",
              help="输入文件的文件夹路径，路径中不能包含中文")
@click.option("-B", "--herb", 'using_herb', default=False, is_flag=True, help="是否分析HERB库,布尔值标记[加上参数即开启]")
@click.option("-M", "--water_mask", 'water_mask', default=False, is_flag=True, help="是否添加水印,布尔值标记[加上参数即开启]")
@click.option("-co", "--cloud_off", 'cloud_off', default=True, is_flag=True, help="是否云交付,布尔值标记[默认开启，加上参数关闭]")
@click.option("-oe", "--oereport_only", 'oereport_only', default=False, is_flag=True, help="只执行函数oebio_report")
@click.option("-local", "--local", 'local', default=False, is_flag=True, help="执行本地项目需要加上此参数")
def dev(input_folder, using_herb=False, water_mask=False, cloud_off=True, oereport_only=False, local=False):
    # from TCMreport import Report
    # r = Report(data=input_folder, local=local)
    from TCMreport import BloodReport
    r = BloodReport(data=input_folder, project="Blood", local=local)
    r.relation_analysis()
    # r.oebio_report(using_herb=True, using_product=True, cloud=False, skip=True)
    # r.preprocess_qi(using_herb=True)
    #r.diff_analysis()
    # interested_cpd
    # r.preprocess_qi(using_herb=True, interested_cpd=["Geniposide", ])
    #r.plot_xic(interested_cpd=["Geniposide",])


    # if oereport_only:
    #     r.oebio_report(using_herb=True, cloud=cloud_off)
    # else:
    #     r.run_pipeline(using_herb=True, water_mask=water_mask, cloud=cloud_off)

if __name__ == '__main__':
    cli()
