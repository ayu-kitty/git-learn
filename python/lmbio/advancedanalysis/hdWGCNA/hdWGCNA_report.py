#! /opt/conda/bin/python
import os
import re
import yaml
import json
import glob
from oebio.report import Report
from collections import defaultdict
import argparse
import shutil

parser = argparse.ArgumentParser(description="hdWGCNA生成report脚本")
parser.add_argument("-n", "--contract", help="合同号。DOE20231111。如果不填写，则默认读取通用的config yaml中的字符串。")
parser.add_argument("-w", "--workdir", help="运行路径，如./ ",default="./")
parser.add_argument("-c", "--pdf2png", help="是否运行pdf向png转换，可选参数[ y/n ]。")
parser.add_argument("-o", "--output_prefix", help="输出html的名字，默认report.html。",default="report")
parser.add_argument("-p", "--parameter", help="分析所用的parameter存储log文件，默认parameter.log。",default="parameter.log")
parser.add_argument("-src", "--srcPath", help="src文件夹的调取路径(不需要带 '/')。默认/data/hstore1/database/report/2023-03-04~2023/hdWGCNA",default="/data/hstore1/database/report/2023-03-04~2023/hdWGCNA")
args = parser.parse_args()

program_path = args.workdir #"./"
pdf2png = args.pdf2png
output_prefix = args.output_prefix
parameter_df = args.parameter
src_dir = args.srcPath
contract = args.contract

# 打开文件并读取内容
with open(parameter_df, 'r') as file:
    lines = file.readlines()
    # 查找包含 "topN1:" 的行
    species_line = next(line for line in lines if "species:" in line)
    species = str(species_line.split(":")[1].split('"')[0].strip())
    hub_n_genes_line = next(line for line in lines if "hub_n_genes:" in line)
    hub_n_genes = int(hub_n_genes_line.split(":")[1].split('"')[0].strip())
    networkType_line = next(line for line in lines if "networkType:" in line)
    networkType = str(networkType_line.split(":")[1].split('"')[0].strip())
    cormethod_line = next(line for line in lines if "cormethod:" in line)
    cormethod = str(cormethod_line.split(":")[1].split('"')[0].strip())


print("spcies:", species)
print("hub_n_genes:", hub_n_genes)
print("networkType:", networkType)
print("cormethod:", cormethod)


# 获取绝对路径，用于后面加入到config.yaml文件
current_working_directory = os.path.abspath('.')

if os.path.exists(program_path+"/src") and os.path.isdir(program_path+"/src"):
    # 删除文件夹及其内容
    shutil.rmtree(program_path+"/src")
    print(f"当前路径存在src文件夹，先进行deleted, 重新copy")
else:
    print(f"当前路径不存在src文件夹,进行copy")
    
shutil.copytree(src_dir + "/src", program_path + "/src")
shutil.copy(src_dir + "/hdwgcna.yaml", "config.yaml")


config = yaml.load(open("config.yaml"),Loader=yaml.FullLoader)
config["html_header"] = current_working_directory + "/src/header/header.jpg"

# 保存更新后的yaml文件，但是变量会自动重新排序......
with open("config.yaml", "w", encoding="utf-8") as file:
    yaml.dump(config, file, default_flow_style=False, allow_unicode=True)


# Define header info
header = dict(config['header'])
header['订单编号'] = contract

def get_file(input):
    file = list(glob.glob(input))
    file = [f.replace('\\', '/') for f in file][0]
    return file


####################################################################################################
## ========================= 一、 1.Create directory  =========================
####################################################################################################
#def biomarker_report(program_path, configfile):
os.chdir(program_path)

if pdf2png == "y":
    os.system('for i in `find ./{}/* -name "*.pdf"`;do new=`echo $i|sed "s/.pdf$//"`;convert $i -background white -alpha remove -flatten -density 300 $new.png;done' .format(output_prefix))

# Define header info
report = Report('hdWGCNA分析', title='hdWGCNA分析', header_info=dict(header),
					oe_welcome='''
    ####  感谢您选择鹿明生物！

    **关于网页报告使用有以下几点提示:**

    1. **报告正文**中的**蓝色字体**均可以**点击**快速导引到感兴趣的章节，便于您快速查看相应内容。
    2. **报告正文**中的**图片**均可以点击后进行**放大查看**，且图片点击后可以包含更多的细节信息，比如左上角会显示具体的**图片数量**，
    右上角有相关的**图片工具**，**键盘上的左右方向键/鼠标滚轮**可以对图片进行**切换**。
    3. **报告正文**中**每个区域**中图片的数量**最多**只显示2张，更多的图片可以在点击之后在弹出界面利用**键盘方向键/鼠标滚轮**查看。
    4. 展示的表格均可以在各列**表头**进行**升序或者降序显示**，当表格多于20行的时候，表格会嵌入到网页之中，可以利用滚动栏进行调整设置，
    同时会在表格上方增加搜索设置，便于快速查询表格信息。
    5. 在报告中浏览，需要返回顶部的时候，可以使用网页**右下角的白色箭头标识**快速返回顶部。
    6. 本提示可以点击右上方的不再显示，则后续打开网页均不会显示本提示。
	''')
report.add_yaml_config(os.path.join(program_path, 'config.yaml'))
#os.environ['OEBIO'] = os.path.join(program_path,"src")
os.environ['OEBIO'] = "%s/src" % program_path
############################################################################
## ========================= 一、总概述 =========================
###############################################################################
project_module = report.add_section('hdWGCNA分析分析概述')
project_module.add_comment('WGCNA是一种从高通量的表达数据中挖掘模块（module）信息的算法，在该方法中module被定义为一组具有类似表达谱的基因，'
                           '并探索module内基因网络与研究者关注的表型之间的关联关系。但WGCNA对于稀疏矩阵的敏感性性能相对较差，'
                           '单细胞数据或者空间组学数据中固有的稀疏性和噪声会导致虚假的基因-基因相关性。基于WGCNA开发的hdWGCNA R包[1]，适用于在单细胞RNA-seq、空间转录组学、空间代谢组学等高维组学。由于单细胞或空间组学数据的相关结构对于不同的子集（细胞类型、细胞状态、解剖区域）有很大差异，所以hdWGCNA流程考虑了这些因素，将高度相似的细胞整合成 "metacells元细胞"或者"metapixels元像素点"，以减少稀疏性，同时保留细胞的异质性。 '
                           '其处理流程可参照下图：'
                           ''
                           '')
project_module.add_plot('%s/src/images/hdwgcna_workflow.png' % program_path,caption ='hdWGCNA分析流程图',content='hdWGCNA分析流程图')
project_module.add_comment('hdWGCNA workflow详细流程解释如下：')
project_module.add_comment('Step1：数据处理'
'在将数据导入到hdWGCNA分析前，一般会对数据做一定处理。可以使用Seurat、Scanpy等R包进行处理。'
'- 质控，比如单细胞转录组数据中，对于线粒体占比高、doublet等细胞的处理；'
'- 数据标准化，在对数据过滤后，一般会进行Normalization操作，如Seurat默认使用通用的“LogNormalize()”方法；'
'其它如workflow中提到的“Feature selection”“Dimensionality redution”、“Batch correction”、“Clustering and annotation”都可以针对具体的项目需求，选择性进行处理。')

project_module.add_comment('Step2：metacells/metapixels的构建'
'这一步也是hdWGCNA中最核心的一步。Metacells定义为：代表不同细胞状态，但feature特征相似的细胞，metapixels类似。'
'基于KNN bagging方法，将相似的cell进行merge。或者将邻近ST spots进行merge。减少了数据的稀疏性，因此也增加了每个“data point”的信息。这样降低了metacells表达谱的噪音。')

project_module.add_comment('Step3：设定表达矩阵'
'基于metacells或者metapixel的定义，得到cell/spot/pixel population矩阵。同时可以过滤掉一些在cells或者spots中一些方差为零的feature，构建矩阵用于下游的Network的构建。')

project_module.add_comment('Step4：Soft power阈值的筛选'
'同WGCNA一样，需要进行软阈值的筛选，软阈值用于减少计算基因-基因相关性的噪音。')

project_module.add_comment('Step5：构建相关性网络'
'在确定软阈值后，进一步进行网络构建、模块相关性分析、TOM矩阵等分析。')

project_module.add_comment('Step6：计算模块特征基因以及模块关联度分析'
'在相关性网络中，一般会重点关注“hub gene“，即和模块高度关联的基因（Module Eigengene）。kME（eigengene connectivity）为特征基因连接度，用于筛选“hub gene”。')

## ========================= 概述  =========================
#res_module = project_module.add_section('hdWGCNA分析结果解析')
project_module.add_comment('hdWGCNA分析结果共分为7个模块，详见下方模块内容解释。' )
project_module.add_comment('本项目富集分析使用**物种缩写编号为：{} .(注：如果为ko，代表使用的为不分物种的总库)。**' .format(species))

## ========================= 1.1 进行metapixel后的降维可视化 =========================
res2_module = project_module.add_section('Metapixel整合及网络构建')
res2_module.add_comment('对数据进行metapixel处理，对处理后的数据，再次进行tSNE和UMAP降维聚类分析，如下：')
res2_module.add_plot('%s/1.WGCNA_setup/Metapixel_tSNE.png' % output_prefix, caption = 'metapixel整合后tSNE降维可视化',content='metapixel整合后tSNE降维可视化',
                     description='图片说明：图中每个点代表Metapixel整合后的点，将所有点进行tSNE降维。')
res2_module.add_plot('%s/1.WGCNA_setup/Metapixel_UMAP.png' % output_prefix, caption = 'metapixel整合后UMAP降维可视化',content='metapixel整合后UMAP降维可视化',
                     description='图片说明：图中每个点代表Metapixel整合后的点，将所有点进行UMAP降维。')
res2_module.add_comment('同时基于metapixel计算后的数据，进行最佳power值的选取分析。'
'基因共表达网络是无尺度（scale-free）的权重基因网络。无尺度特性，又称作无标度特性，是指网络的度分布满足幂律分布。幂律分布这一特性，正说明无尺度网络的度分布是呈集散分布：大部分的节点只有比较少的连接，而少数节点有大量的连接。'
'为了尽量满足无尺度网络分布前提条件，需要选择邻接矩阵权重参数power的取值。设定power值从1-30，分别计算其得到的网络对应的相关系数和网络的平均连接度。其中，相关系数越高（最大为1），网络越接近无尺度网络分布，但同时还需要保证一定的基因连接度，所以这个power值在相关系数足够大的同时保证基因的连接度较大。'
'选取的最佳power值如下图所示，图中黑色圆圈中的数值即为最佳power值：')
res2_module.add_plot('%s/1.WGCNA_setup/SelectPower.png' % output_prefix , caption = 'WGCNA网络构建参数（soft power）',content='WGCNA网络构建参数（soft power）',
                     description='图片说明：上图展示不同power值下的Topology Model Fit，以及对应的网络的最大、最小即平均连接度。黑色圆圈内的数值为选取的最佳power值。')

res2_module.add_table('%s/1.WGCNA_setup/SelectPower.xlsx' % output_prefix , caption = 'WGCNA网络构建参数（soft power）表格',
                     description='Power：power值; new_SFT.R.sq：计算的soft R-squared值。最佳power值的选取主要参照 new_SFT.R.sq 列。'
                     
                     )
res2_module.add_comment("分析结果详见 ：**{link} 文件夹**", link = "%s/1.WGCNA_setup" % output_prefix)

## ========================= 1.2 最佳阈值网络 =========================
###############################################################################
res3_module = project_module.add_section('最佳阈值网络')
res3_module.add_comment('在挑选最佳软阈值后，进行Network构建，并进行Dendrogram以及TOM（即相关性矩阵）的可视化展示分析。')
res3_module.add_comment('网络的构建可以选取singed以及unsiged，本项目采用 {}。hdWGCNA在构建相关性矩阵的时候，一般关系数范围是-1~1，'
                        'WGCNA分析要求转换为0-1范围。有两种转换方式：unsigned:不区分正相关和负相关，signed：区分正负相关。' .format(networkType))
res3_module.add_plot('%s/2.bestPower/Dendrogram.png' % output_prefix, caption = 'Dendrogram图',content='Dendrogram图',description=
                      '图片说明：图中上部分为对加权后的相关系数构建的dissTOM矩阵构建的基因聚类树，图中下部分为每个模块基因的分布情况，'
                      '同一种颜色表示同一个模块，这些模块用于后续的分析。'
                      )
res3_module.add_comment('同时基于TOM矩阵，进行可视化展示：')
res3_module.add_plot('%s/2.bestPower/TOM-NetHeatmap.png' % output_prefix, caption = 'TOM图',content='TOM图')
res3_module.add_comment("分析结果详见 ：**{link} 文件夹**", link = "%s/2.bestPower" % output_prefix)

## ========================= 1.3 模块内Eigen-feature =========================
###############################################################################
res4_module = project_module.add_section('模块内特征feature分析')
res4_module.add_comment('在hdWGCNA中，一般重点关注”Hub-Feature”。 Hub-Feature即和模块有高关联度的这些特征。')
res4_module.add_comment('**Hub-Feature：**即模块特征feature，或者模块的核心feature；')
res4_module.add_comment('')
res4_module.add_comment('**kME (eigengene connectivity)：** 即特征基因连接度，根据kMEs用于筛选Hubgene；')
res4_module.add_comment('\n')
res4_module.add_comment('模块特征可以通过以下3个维度来评估，并用于模块计算，比如模块间相关性分析、模块和性状相关性分析等：')
res4_module.add_comment('-  **MEs（Module Eigenes）**：模块特征基因(MEs)，是一种常用的度量，用来总结整个共表达模块的基因表达谱。')
res4_module.add_comment('-  **average expression**：模块内基因平均表达。')
res4_module.add_comment('-  **module score**：模块内基因评分，用于计算每个模块给定数量的基因的基因分数。基因评分是通过计算模块特征基因来总结模块表达的另一种方法。')
res4_module.add_comment('在此分析模块，会基于以下步骤进行分析结果展示：')
res4_module.add_comment('（1）基于kMEs，找到 Hub-Feature ，并进行kMEs分析；')
#res4_module.add_comment('(2) 模块特征分析：以umap降维投影方式，分别通过average expression（模块feature平均表达）、MEs、hMEs、Module Score进行展示；')
res4_module.add_comment('（2）模块相关性分析：以模块间相关性结果展示，'
                             '分别通过average expression（模块feature平均表达）、MEs、Module Score进行相关性计算；')
res4_module.add_comment('（3）模块趋势分析：以小提琴图、气泡图、拟合线图形式，分别通过average expression（模块feature平均表达）、MEs、Module Score进行展示。')
res4_module.add_comment('本项目挑选各模块 **top {}** 的基因为Hubgene。' .format(hub_n_genes))


## ========================= 3.1 模块kMEs图 =========================
res41_module = res4_module.add_section('特征feature连接度分析')
res41_module.add_comment('kME (eigengene connectivity) 为特征基因连接度，用于筛选Hubgene。对各模块Hub-Feature的kMEs进行可视化展示如下:')
res41_module.add_plot('%s/3.ModuleEigengenes/*kME_allFeatures.png' % output_prefix, caption = '特征基因连接度（kMEs）图',description=
                      '图片说明：图中为各模块的特征基因（feature）连接度，x轴为各模块对应的基因，y轴为模块内基因的kMEs数值，同时挑选了各模块top特征feature列出。'
                      )
res41_module.add_comment('**对各模块内的总的基因列表进行展示，详见 excel表格 : {link}**', link = "%s/3.ModuleEigengenes/*allFeatures_inModules.xlsx" % output_prefix)
res41_module.add_comment('**同时对各模块内的Hubgene列表进行展示，详见 excel表格：{link}**',link = "%s/3.ModuleEigengenes/1_top25_HubFeatures_inModules.xlsx" % output_prefix)
res41_module.add_comment("**分析结果详见 ：{link} 文件夹**。", link = "%s/3.ModuleEigengenes" % output_prefix)


## ========================= 1.4.2 模块内特征feature MEs降维 umap图 =========================

# res42_module = res4_module.add_section('模块特征降维 umap图')
# res42_module.add_comment('对各模块Hub-Feature，从4个维度（average expression、hMEs、MEs、module score），进行降维UMAP可视化展示如下:')
# res42_module.add_plot('%s/3.ModuleEigengenes/*_umap.png' % output_prefix, caption = '模块特征降维UMAP图',description=
                      # '图中为各模块的特征基因（feature）连接度的UMAP降维。不同颜色代表不同的模块，颜色的深浅代表数值大小。'
                      # )
# res42_module.add_comment("分析结果详见 ：{link} 文件夹， 对应文件夹中 2.1~2.4部分 ", link = "%s/3.ModuleEigengenes" % output_prefix)

                      

## ========================= 3.2 模块相关性热图 =========================
res43_module = res4_module.add_section('模块相关性')
res43_module.add_comment('分别从3个维度（average expression、MEs、module score），计算每个共表达模块间的相关性，并进行可视化展示：')
res43_module.add_plot('%s/3.ModuleEigengenes/*Module_Corplot_*.png' % output_prefix, caption = '模块间相关性热图',description=
                      '图片说明：上图为模块间相关性热图，粉红色代表正相关，蓝色代表负相关。颜色越深同时形状越趋于一条直线，代表相关性数值越高；'
                      '颜色越浅同时形状越趋于圆形，代表相关性数值越低。如果相关性p值不显著，则在图中标注一个黑色的叉。'
                      )
res43_module.add_comment("**分析结果详见 ：{link} 文件夹，对应文件夹中 3.2部分。**", link = "%s/3.ModuleEigengenes" % output_prefix)


## ========================= 3.3 模块特征feature丰度趋势 =========================
res44_module = res4_module.add_section('模块特征feature趋势')
res44_module.add_comment('提取各模块Hub Feature的特征（average expression、hMEs、MEs、module score），'
                         '并基于特征值进行小提琴图和气泡图的绘制。')
res44_module.add_plot('%s/3.ModuleEigengenes/*_violinAll*.png' % output_prefix, caption = '模块Hub Feature趋势小提琴图',description=
                      '图片说明：小提琴图中箱线显示了数据的中位值、上四分位、下四分位等信息。'
                      )
res44_module.add_comment('同时进行气泡图的绘制：')
res44_module.add_plot('%s/3.ModuleEigengenes/*Module2Cluster_DotPlot*.png' % output_prefix , caption = '各模块Hub Feature趋势气泡图'
                      ,description='图片说明：上图为各模块Hub Feature趋势气泡图，x轴为分组，y轴为不同的模块，气泡大小对应值的高低。'
                      )

res44_module.add_comment('同时基于模块内的average expresssion，即Hub Feature的丰度平均值，进行拟合线图的绘制如下：')
res44_module.add_plot('%s/3.ModuleEigengenes/*hubFitline_all_Expr.png' % output_prefix, caption = '模块特征feature表达趋势拟合图',content='模块特征feature表达趋势拟合图',
                      description='图片说明：上图为表达趋势拟合图，采用loess拟合方法。'
                      )
res44_module.add_comment("**分析结果详见 ：{link} 文件夹， 对应文件夹中 3.3 部分**", link = "%s/3.ModuleEigengenes" % output_prefix)


## ========================= 1.5 各模块特征feature网络 =========================
###############################################################################
res5_module = project_module.add_section('各模块特征feature网络')
res5_module.add_comment('WGCNA构建的网络，是无尺度网络。正如前面所述，无尺度网络图中，绝大多数节点（基因）的连接度较低，只有少数节点连接度较高。一个模块内的核心基因，即该模块中连接度最高的一群基因。'
'本次模块核心基因分析，分别分析了每个模块中连接度最高的top特征feature，展示这些特征feature之间的关系；具体分析结果见：')
res5_module.add_plot('%s/4.ModuleNetwork/HubGeneNetwork.png' % output_prefix, caption = '模块特征feature总网络图',content='模块特征feature总网络图',description=
                      '图片说明：网络图中的线表示基因与基因之间的联系程度，基因与周边点的联系越多，在网络中越处于核心地位。'
                      )
res5_module.add_plot('%s/4.ModuleNetwork/*_top*.png' % output_prefix, caption = '各模块top特征feature网络图',description=
                      '图片说明：分别对各模块进行网络图绘制，网络图中的线表示基因与基因之间的联系程度，基因与周边点的联系越多，在网络中越处于核心地位。'
                      )
res5_module.add_comment("**分析结果详见 ：{link} 文件夹。** ", link = "%s/4.ModuleNetwork" % output_prefix)

## ========================= 1.6 模块间差异分析 =========================
###############################################################################
res6_module = project_module.add_section('模块间差异分析')
res6_module.add_comment('基于各模块内的Hub-Feature，分别在各模块内，计算模块内各不同分组间的差异。'
                        '差异计算的方法为：计算某1组和其他所有组的差异。'
                        '详细结果如下：'
)
res6_module.add_plot('%s/5.diffGroup_in_Module/*VolcanoPlot.png' % output_prefix, caption = '模块差异火山图',description=
                      '图片说明：上图为模块差异火山图，x轴为差异倍数，y轴为 −log10(Adj.P−value)'
                      )
res6_module.add_plot('%s/5.diffGroup_in_Module/*LollipopPlot.png' % output_prefix, caption = '模块差异棒棒糖图',description=
                      '图片说明：上图为模块差异棒棒糖图，x轴为Avg. log2 (Fold Change)。'
                      )
res6_module.add_comment("**分析结果详见 ：{link} 文件夹。**", link = "%s/5.diffGroup_in_Module" % output_prefix)

## ========================= 1.7 模块和性状相关性 =========================
###############################################################################
res7_module = project_module.add_section('模块和性状相关性')
res7_module.add_comment('基于模块内Feature平均丰度，计算模块和各性状的相关性, 并绘制相关性heatmap如下:'

)
res7_module.add_plot('%s/6.ModuleTraitCorrelation/ModuleTrais-cor-p-heatmap.png' % output_prefix, caption = '模块性状相关性heatmap',content='模块性状相关性heatmap',
                       description=
                      '图片说明：上图为模块性状相关性heatmap，红色代表正相关，蓝色代表负相关，颜色越深，相关性系数越大。'
                      )


## ========================= 1.8 模块富集分析 =========================
###############################################################################
res8_module = project_module.add_section('模块富集分析')
res8_module.add_comment('KEGG( https://www.kegg.jp/ )数据库的构建旨在了解生物系统（如细胞，组织等）中基因、蛋白及代谢物的功能及相互作用关系。'
'可以查询到与代谢物相关的代谢通路、人类疾病及药物研发等信息。该数据库的代谢物及代谢通路涉及两大类：真核生物（动物、植物、真菌及原生生物）和原核生物（细菌、古细菌）。'
'通过对差异代谢物进行通路富集分析，有助于理解在差异样品中代谢途径变化机制。基于KEGG数据库对差异代谢物进行代谢通路富集分析。'
)
res8_module.add_comment('利用差异代谢物的KEGG ID进行通路富集分析，获得代谢通路富集结果。'
'应用超几何检验，找出与整个背景相比，在显著性差异表达代谢物中显著富集的pathway条目，其计算公式为：'
)
res8_module.add_comment('**注：富集分析部分采用各模块总的features进行富集分析分析**。'
)
res8_module.add_plot('%s/src/images/超几何分布.png' % program_path
                      )
res8_module.add_comment('N为代谢物总数；n为N中差异表达代谢物的数目；M为注释为某特定pathway的代谢物数目；m为注释为某特定pathway的差异代谢物数目。'
'以p-value≤0.05为阈值，满足此条件的pathway为在差异代谢物中显著富集的pathway。p-value越小，则该代谢通路的差异性越显著。'
)
res8_module.add_comment('**注：对各个模块进行富集分析，富集分析所用基因为模块内所有基因。**'
)
res8_module.add_plot('%s/7.Enrichment/*M1_KEGGlist/代谢通路富集图-top20.png' % output_prefix, caption = 'TOP-20代谢通路富集图',content='TOP-20代谢通路富集图',description=
                      '图片说明：代谢通路中p-value为该代谢通路富集的显著性。红线示意p值为0.01，蓝线示意p值为0.05，条柱的顶端高于蓝线时，其所代表的信号通路具有显著性。'
                      )
res8_module.add_plot('%s/7.Enrichment/*M1_KEGGlist/代谢通路富集图-p小于0.05.png' % output_prefix, caption = 'p<0.05代谢通路富集图',content='p<0.05代谢通路富集图',description=
                      '图片说明：代谢通路中p-value为该代谢通路富集的显著性。红线示意p值为0.01，蓝线示意p值为0.05，条柱的顶端高于蓝线时，其所代表的信号通路具有显著性。'
                      )
res8_module.add_plot('%s/7.Enrichment/*M1_KEGGlist/气泡图-top20.png' % output_prefix, caption = 'TOP-20气泡图',content='TOP-20气泡图',description=
                      '图片说明：代谢通路中p-value为该代谢通路富集的显著性，选择显著性富集pathway进行气泡图绘制。纵坐标为代谢通路名称；'
                      '横坐标为富集因子（Rich factor，Rich factor=显著差异代谢物个数/该pathway中的总代谢物个数），'
                      'Rich factor越大，则说明富集程度越大；颜色由绿到红表示p-value依次降低；点越大，说明富集到该pathway上的代谢物数目越多。'
                      )
res8_module.add_plot('%s/7.Enrichment/*M1_KEGGlist/气泡图-p小于0.05.png' % output_prefix, caption = 'p<0.05气泡图',content='p<0.05气泡图',description=
                      '图片说明：代谢通路中p-value为该代谢通路富集的显著性，选择显著性富集pathway进行气泡图绘制。纵坐标为代谢通路名称；'
                      '横坐标为富集因子（Rich factor，Rich factor=显著差异代谢物个数/该pathway中的总代谢物个数），'
                      'Rich factor越大，则说明富集程度越大；颜色由绿到红表示p-value依次降低；点越大，说明富集到该pathway上的代谢物数目越多。'
                      )
res8_module.add_comment("**分析结果详见 ：{link} 文件夹。**", link = "%s/7.Enrichment" % output_prefix)


## ========================= 二 、文献软件及数据库信息  =========================
lit_module = report.add_section('文献参考')
lit_module.add_comment('[1] Morabito, Samuel, et al. "hdWGCNA identifies co-expression networks in high-dimensional transcriptomics data." Cell Reports Methods (2023).')
## ========================= 2.6 SECTION 6 申明 =========================
shenming_module = report.add_section('申明')
shenming_module.add_comment('  本项目报告由上海鹿明生物科技有限公司提供给项目相关客户。本公司承诺：未经客户同意，不向第三方泄露数据及数据'
                            '分析内容，不将客户数据用于任何商业行为（遵循合同保密协议）。客户未经本公司同意，不得以任何目的向第三方出示项目报告。'
                            '本报告的最终解释权归上海鹿明生物科技有限公司。')
## =========================   3.  写入html     =========================

# Generate Report HTML
xx = output_prefix + ".html"
report.write_to(xx)

#if __name__ == "__main__":
#	scrna_report()
