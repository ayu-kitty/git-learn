#!/opt/conda/bin/python
#coding=utf-8
from lmbio.basic.runsnakemake import runsnakemake
from lmbio.basic.printfilter import FilteredPrinter,filtered_print,filtered_print2
from lmbio.basic import path
import pandas as pd
import argparse
import snakemake
import shutil
import glob
import sys
import re
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-ip", "--inputpath",default = "rawdata/",help = "输入文件路径，默认rawdata/")
    parser.add_argument("-s","--savepath",default = "KSEA_Activity/",help = "结果保存路径，KSEA_Activity/")
    parser.add_argument("-or", "--org",default = "hsa",help = "物种kegg缩写，默认hsa")
    parser.add_argument("-pt","--diffilephos",default = "差异表达矩阵(未筛选)-ptm.xlsx", help = "磷酸化表格名称，默认为差异表达矩阵(未筛选)-ptm.xlsx")
    parser.add_argument("-da","--inputdata",default = "KSEA Kinase Scores.xlsx", help = "输入绘图文件，默认KSEA Kinase Scores.xlsx")
    parser.add_argument("-font","--fontfamily",default = "sans", help = "富集绘图字体，默认sans")
    parser.add_argument("-fc","--foldchange",default = 1.2,type = float, help = "差异筛选FC值，默认1.2")
    parser.add_argument("-c","--cores",default = 30, type= int,help = "线程数，默认30")
    parser.add_argument("-sm","--snakefile",default = "snakemake/Diffanalysis/diffgetdata/KSEA_plot.smk",help = "snakefile路径")
    
    args = parser.parse_args()
    
    def collect_reads(path,diffilephos):
        group = pd.ExcelFile(path + diffilephos).sheet_names
        return group

    groups = collect_reads(args.inputpath,args.diffilephos)
    
    config = {"path":args.inputpath,
              "ptmdata":args.diffilephos,
              "savepath":args.savepath,
              "drawdata":args.inputdata,
              "org":args.org,
              "compare":groups,
              "fontfamily" : args.fontfamily}
    
    runsnakemake(snakefilepath = args.snakefile,
                 config = config,
                 configfiles =[],
                 cores = args.cores,
                 targets = ["Ksea"])

    os.system("find -name .snakemake_timestamp -print0  | xargs -0  rm -rf")
    shutil.rmtree(".snakemake")

