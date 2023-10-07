#!/opt/conda/bin/python
# File: commands
# Author: jiang tao
# Time: 2023/3/23 15:43
import click

from click import ParamType


class MultiType(ParamType):
    """自定义selectmz的参数类型"""
    name = 'multiType'

    def convert(self, value, param, ctx):
        try:
            value = str(value)
            if "," in value:
                return [v.strip() for v in value.split(",")]
            elif "/" in value:
                with open(value, "r", encoding='utf8') as f:
                    return [mz.strip() for mz in f.read().split("\n") if mz.strip()]
            else:
                return None
        except Exception as e:
            return e


@click.group()
def cli():
    """\b
        空代共定位命令行工具
        操作步骤：
            1.运行`image_align`检查对齐结果，未对齐手动调参直到达到最佳效果
            2.运行`image_correlation`
        """
    ...


@cli.command('image_align', short_help="空代HE染色图与切片图对齐[自动+手动], 需人工检查对齐结果")
@click.option("-O", "--output", 'output_folder', type=click.Path(), required=True, help="输出结果的文件夹路径，不存在会自动创建")
@click.option("-I", "--imzml", 'imzml_path', type=click.Path(exists=True), required=True, help="imzml文件路径")
@click.option("-H", "--he", 'he_path', type=click.Path(exists=True), required=True, help="HE染色图文件路径")
@click.option("-it", "--image_type", 'image_type', type=click.Choice(['png', 'jpg']), default='png', show_default=True, help="图片保存格式")
@click.option("-ps", "--ps", 'ps', type=float, default=5, show_default=True, help="参考点,默认5")
@click.option("-tv", "--tv", 'tv', type=float, default=0.1, show_default=True, help="色差,默认0.1")
@click.option("-ba", "--bg_alpha", 'bg_alpha', type=float, default=0, show_default=True, help="背景透明度, 默认0")
@click.option("-fa", "--fig_alpha", 'fig_alpha', type=float, default=0.2, show_default=True, help="图像透明度,默认0.2")
@click.option("-pt", "--pixel_thresh", 'pixel_thresh', type=float, default=0.85, show_default=True, help="背景识别阈值,默认0.85")
@click.option("-n", "--num", 'num', type=int, default=5, show_default=True, help="切除边界像素数,默认5")
@click.option("-xa", "--xmax_adj", 'xmax_adj', type=int, default=0, show_default=True, help="右边界偏差,默认0")
@click.option("-xm", "--xmin_adj", 'xmin_adj', type=int, default=0, show_default=True, help="左边界偏差,默认0")
@click.option("-ya", "--ymax_adj", 'ymax_adj', type=int, default=0, show_default=True, help="上边界偏差,默认0")
@click.option("-ym", "--ymin_adj", 'ymin_adj', type=int, default=0, show_default=True, help="下边界偏差,默认0")
@click.option("-xa_", "--xmax_adj_", 'xmax_adj_', type=int, default=0, show_default=True, help="旋转后,右边界偏差,默认0")
@click.option("-xm_", "--xmin_adj_", 'xmin_adj_', type=int, default=0, show_default=True, help="旋转后,左边界偏差,默认0")
@click.option("-ya_", "--ymax_adj_", 'ymax_adj_', type=int, default=0, show_default=True, help="旋转后,上边界偏差,默认0")
@click.option("-ym_", "--ymin_adj_", 'ymin_adj_', type=int, default=0, show_default=True, help="旋转后,下边界偏差,默认0")
@click.option("-ita", "--interaction", 'interaction', default=False, is_flag=True, help="是否交互,布尔值标记[加上参数即开启]")
@click.option("-ra", "--rotate_angle", 'rotate_angle', type=int, default=0, show_default=True, help="旋转角度,默认0")
@click.option("-fi", "--flip", 'flip', default=False, is_flag=True, help="是否上下翻转,布尔值标记[加上参数即开启]")
@click.option("-fo", "--flop", 'flop', default=False, is_flag=True, help="是否左右翻转,布尔值标记[加上参数即开启]")
@click.option("-r", "--resize", 'resize', default=False, is_flag=True, help="是否重设图片大小,布尔值标记[加上参数即开启]")
@click.option("-fw", "--fig_width", 'fig_width', type=int, default=1000, show_default=True, help="图片宽度")
@click.option("-fh", "--fig_height", 'fig_height', type=int, default=1000, show_default=True, help="图片高度")
def image_align(output_folder, imzml_path, he_path, image_type, ps, tv, bg_alpha,
                fig_alpha, pixel_thresh, num, xmax_adj, xmin_adj, ymax_adj, ymin_adj,
                xmax_adj_, xmin_adj_, ymax_adj_ , ymin_adj_, interaction, rotate_angle,
                flip, flop, resize, fig_width, fig_height):
    """空代HE染色图与切片图对齐[自动+手动]"""
    from spatialMeta import ColocalizationAnalysis

    ca = ColocalizationAnalysis(output_folder, imzml_path, he_path)
    df = ca.getCoordinates(imzml_path)
    fig_file = ca.heMatch(imzml_path, he_path, df, image_type=image_type)
    ca.heMatchManual(imzml_path, he_path, fig_file,
                     image_type=image_type,
                     ps=ps,
                     tv=tv,
                     bg_alpha=bg_alpha,
                     fig_alpha=fig_alpha,
                     pixel_thresh=pixel_thresh, num=num,
                     ymax_adj=ymax_adj, ymin_adj=ymin_adj,
                     xmax_adj=xmax_adj, xmin_adj=xmin_adj,
                     ymax_adj_=ymax_adj_, ymin_adj_=ymin_adj_,
                     xmax_adj_=xmax_adj_, xmin_adj_=xmin_adj_,
                     interaction=interaction,
                     rotate_angle=rotate_angle,
                     flip=flip,
                     flop=flop,
                     resize=resize,
                     fig_width=fig_width,
                     fig_height=fig_height
                     )


@cli.command('image_correlation', short_help="代谢物成像图 raw特征、resnet特征相关性分析")
@click.option("-O", "--output", 'output_folder', type=click.Path(), required=True, help="输出结果的文件夹路径，不存在会自动创建")
@click.option("-I", "--imzml", 'imzml_path', type=click.Path(exists=True), required=True, help="imzml文件路径")
@click.option("-H", "--he", 'he_path', type=click.Path(exists=True), required=True, help="HE染色图文件路径")
@click.option("-s", "--selectmz", 'selectmz', type=MultiType(),  default=None, help="筛选的mz,可以是英文逗号分割的字符串或者 mz文件(每行一个mz)")
@click.option("-c", "--corfilter", 'corfilter', type=float, default=0.7, show_default=True, help="相关性R2筛选阈值")
@click.option("-ct", "--corfiltertype", 'corfiltertype', type=click.Choice(["+", "-", "+-"]), default="+",
              show_default=True, help="相关性筛选方法,包含+-,+,-")
@click.option("-p", "--pfilter", 'pfilter', type=float, default=0.05, show_default=True, help="显著性P值筛选阈值")
@click.option("-a", "--adjust", 'adjust',
              type=click.Choice(["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"]), default="BH",
              show_default=True, help="显著性校正方法,默认none,包括holm,hochberg,hommel,bonferroni,BH,BY,fdr")
@click.option("-cm", "--cormethod", 'cormethod', type=click.Choice(["pearson", "spearman", "kendall"]),
              default="pearson", show_default=True, help="相关性计算方法,包括pearson, spearman, kendall")
@click.option("-u", "--use", 'use', type=click.Choice(["pairwise", "complete"]), default="pairwise", show_default=True,
              help="相关性计算方法,包括pairwise,complete")
@click.option("-x", "--xname", 'xname', type=str, default="mz-x", show_default=True, help="x数据名称")
@click.option("-y", "--yname", 'yname', type=str, default="mz-y", show_default=True, help="y数据名称")
def image_correlation(output_folder, imzml_path, he_path, selectmz, corfilter, corfiltertype, pfilter, adjust,
                      cormethod, use, xname, yname):
    """代谢物成像图 raw特征、resnet特征相关性分析"""
    from rpy2 import robjects as ro
    from spatialMeta import ColocalizationAnalysis
    from spatialMeta.colocalization.rutils import NULL

    if selectmz is None:
        selectmz = NULL
    elif isinstance(selectmz, list):
        selectmz = ro.StrVector(selectmz)
    else:
        print(str(selectmz))
        return

    ca = ColocalizationAnalysis(output_folder, imzml_path, he_path)
    fig_data = ca.get_fig_data(imzml_path, he_path)
    ca.correlation_raw_metabolite(imzml_path, he_path, fig_data,
                                  selectmz=selectmz,
                                  corfilter=corfilter,
                                  corfiltertype=corfiltertype,
                                  pfilter=pfilter,
                                  adjust=adjust,
                                  cormethod=cormethod,
                                  use=use,
                                  xname=xname,
                                  yname=yname)
    resnet_data = ca.get_resnet_data(imzml_path, he_path, fig_data)
    ca.correlation_resnet(imzml_path, he_path, resnet_data,
                          selectmz=selectmz,
                          corfilter=corfilter,
                          corfiltertype=corfiltertype,
                          pfilter=pfilter,
                          adjust=adjust,
                          cormethod=cormethod,
                          use=use,
                          xname=xname,
                          yname=yname)


if __name__ == '__main__':
    cli()

