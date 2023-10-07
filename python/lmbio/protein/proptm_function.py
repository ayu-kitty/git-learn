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
    parser.add_argument("-t", "--enrichtype",default = "all",help = "富集数据库类型，默认全选，单项示例G/K/I")
    parser.add_argument("-ph","--diffile",default = "去除本底差异数据.xlsx", help = "需要富集的名称，默认为去除本底差异数据.xlsx")
    parser.add_argument("-p","--inputpath",default = "./", help = "差异表格路径，默认当前路径")
    parser.add_argument("-d","--database",default = "background/", help = "富集背景文件路径，默认background/")
    parser.add_argument("-di","--diffout",default = "temp/", help = "diff表格保存路径，默认temp/")
    parser.add_argument("-kd","--keggdb",default = "/data/hstore3/database/kegg/", help="keggmap图所有物种存放位置，默认为/data/hstore3/database/kegg/")
    parser.add_argument("-s","--savepath",default = "Enrichment",help = "富集结果保存路径，默认Enrichment")
    #parser.add_argument("-pp","--ppipath",default = "./",help = "富集结果保存路径，默认./")
    parser.add_argument("-n","--number",default = 25, help = "ppi绘图蛋白top数量，默认25")
    parser.add_argument("-o","--org",default = "hsa", help = "kegg物种拉丁名缩写，默认hsa")
    parser.add_argument("-font","--fontfamily",default = "sans", help = "富集绘图字体，默认sans")
    parser.add_argument("-sm","--snakefile",default = "snakemake/Diffanalysis/diffgetdata/diff_enrich.smk",help = "snakefile路径")
    parser.add_argument("-c","--cores",default = 30, type= int,help = "线程数，默认30")
    
    args = parser.parse_args()
    os.system("makediff -p {path} -f {file} -o {diffout} -om P".format(path=args.inputpath,file=args.diffile,diffout=args.diffout))

    def collect_reads():
      files = glob.glob("temp/*-diff-*.xls")
      files.sort()
      reads = [file.split('/')[-1].split('-diff-')[0] for file in files]
      other = files[0].split('/')[-1].split('-diff-')[1]
      return reads,other
  
    compare = list(set(collect_reads()[0]))

    if args.enrichtype=="all":
      args.enrichtype="G,K,ppi"
    config = {"path":{"projectpath":args.savepath,
                      "projectppipath":args.savepath,
                      "diff_ana_rds_savepath":"",
                      "diff_filter_rds_savepath":"",
                      "enrich_xls_savepath":"temp"},
              "params":{"omic":"P",
                        "all_compare":compare,
                        "enrichtype":args.enrichtype,
                        "enrich":{"org":args.org,
                                  "databasefrom":"/data/hstore3/database/background/",
                                  "database":args.database,
                                  "keggdb":args.keggdb,
                                  "ppitopn":args.number},
                        "plot":{"family":args.fontfamily}}}
    runsnakemake(snakefilepath = args.snakefile,
                 config = config,
                 configfiles =[],
                 cores = args.cores,
                 targets = ["enrich_result"])
    os.system("find -name .snakemake_timestamp -print0  | xargs -0  rm -rf")
    shutil.rmtree(".snakemake",ignore_errors=True)
    shutil.rmtree("temp",ignore_errors=True)

