#! /opt/conda/bin/python
# -*- coding: utf-8 -*-
from lmbio.basic.system import system
import argparse

parser = argparse.ArgumentParser(description="hdWGCNA pipeline\n\n usage: hdWGCNA_pip.py -i file1_df.xls -m file2_group.xls -a annotation.xlsx")
parser.add_argument("-i", "--input", help="输入文件1：表达量矩阵。")
parser.add_argument("-g", "--group", help="输入文件2：Group文件, 第一列为样本列, 第二列为分组列，列名为Group。")
parser.add_argument("-a", "--annotation", help="注释文件, 如果输入矩阵第一列为 mz,或者gene ID等, 需要根据注释文件将mz或者gene ID转换为富集分析用的 C00001 号,或者gene Name, 便于提供富集分析用列表。")
parser.add_argument("-l", "--orders", help="输入既可以是一个只包含一列的文件,也可以是[]内的字符: [ N-epi,OLK-epi,OSCC-nepi,OSCC-ca2 ].指定Group列元素的顺序,绘图则也会按照这个顺序展示. 比如: c('N-epi','OLK-epi','OSCC-nepi','OSCC-ca2')。",default="None")
parser.add_argument("-s", "--species", help="第二种方法选取的biomarker个数。默认top [5] 。",default="ko")
parser.add_argument("-hubN", "--hubgeneN", help="Top N Hubgene,默认为25 。",default=25)
parser.add_argument("-minM", "--minModuleSize", help="构建Module 和 TOM的时候, 允许Module的最小基因数目. 默认为20。",default=20)
parser.add_argument("-deepS", "--deepSplit", help="0~4之间的数字,数字越大,模块个数越多,模块内的基因数目越少.默认为4 。",default=4)
parser.add_argument("-net", "--networkType", help="signed代表构建Network的时候会考虑方向. unsigned/signed hybrid/signed。",default="signed")

parser.add_argument("-cor", "--cormethod", help="相关性算法，可选方法： pearson / spearman / kendall。",default="pearson")
parser.add_argument("-trais", "--traiscol", help="计算相关性时候所用到的性状的列名。默认为all，即分组文件中的所有列。否则输入列名，以逗号分隔，比如 [ p1,p2,p3,p4,p5]。。",default="all")
parser.add_argument("-k", "--kvalue", help="nearest-neighbors parameter。默认为20。注：当数据集比较小的时候，默认的kvalue、min_cells、max_shared这三个参数，需要调整数值更小，否则会报错。。",default=20)
parser.add_argument("-minc", "--min_cells", help="Metapixel中最小的cells数目，需要大于kvalue。默认为25。",default=25)
parser.add_argument("-maxs", "--max_shared", help="maximum number of shared cells between two metacells。默认为10。",default=10)
parser.add_argument("-name", "--wgcnaname", help="WGCNA对象的命名, 比如如果分析的为肿瘤不同恶性程度的Group数据, 可以命名为Malignancy, 会影响部分绘图的命名。",default="WGCNA")


parser.add_argument("-w", "--workdir", help="report部分：运行路径，如./ 。",default="./")
parser.add_argument("-n", "--contract", help="合同号。DOE20231111。如果不填写，则默认读取通用的config yaml中的字符串。")
parser.add_argument("-c", "--pdf2png", help="report部分：是否运行pdf向png转换，可选参数[ y/n ]。",default="y")
parser.add_argument("-o", "--output_prefix", help="report部分：输出html的名字，默认report.html。",default="hdWGCNA_report")
parser.add_argument("-p", "--parameter", help="report部分：分析所用的parameter存储log文件，默认。",default="parameter.log")

args = parser.parse_args()
inputfile = args.input
groupfile = args.group
annotation = args.annotation
orders = args.orders
species = args.species
hubgeneN= args.hubgeneN
minModuleSize= args.minModuleSize
deepSplit= args.deepSplit
networkType= args.networkType
cormethod= args.cormethod
traiscol= args.traiscol
kvalue= args.kvalue
min_cells= args.min_cells
max_shared= args.max_shared
wgcnaname= args.wgcnaname

workdir = args.workdir
pdf2png = args.pdf2png
output_prefix = args.output_prefix
parameter = args.parameter
contract = args.contract

## run biomarker分析
print(">>>> step1 script")
command="flow_hdWGCNA_r  -i %s -g %s -a %s -l %s -s %s -n %s -z %s -d %s -t %s -c %s -r %s -k %s -b %s -x %s -e %s -o %s" %(inputfile,groupfile,annotation,orders,species,hubgeneN,minModuleSize,deepSplit,networkType,cormethod,traiscol,kvalue,min_cells,max_shared,wgcnaname,output_prefix)
print(command)
print(">>>>>>>>>>>>>>>>>")
system(command)

with open("run.sh", "w") as file:
    file.write("# step1: Run hdWGCNA analysis \n"+command+"\n")

## 生成html报告
print(">>>> step2 script")
command="flow_hdWGCNA_report -w %s -c %s -o %s -p %s -n %s" %(workdir,pdf2png,output_prefix,parameter,contract)
print(command)
print(">>>>>>>>>>>>>>>>>")
system(command)

with open("run.sh", "a") as file:
    file.write("# step2: run html report \n"+command+"\n")


