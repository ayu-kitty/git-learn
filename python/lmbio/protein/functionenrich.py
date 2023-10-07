#!/opt/conda/bin/python
from lmbio.basic.runsnakemake import runsnakemake
from lmbio.basic.printfilter import FilteredPrinter,filtered_print,filtered_print2
from lmbio.basic import path
from lmbio.basic.system import system
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
    parser.add_argument("-om", "--omic",default = "P",help = "组学分类，蛋白为P，代谢为M，默认P")
    parser.add_argument("-t", "--enrichtype",default = "all",help = "富集数据库类型，默认全选，单项示例G,K,W,R,I,ppi")
    parser.add_argument("-f","--diffile",default = "差异表达矩阵.xlsx", help = "需要富集的表格名称，默认为差异表达矩阵.xlsx")
    parser.add_argument("-p","--inputpath",default = "./", help = "差异表格路径，默认当前路径")
    parser.add_argument("-d","--database",default = "background/", help = "富集背景文件路径，默认background/")
    parser.add_argument("-kd","--keggdb",default = "", help="keggmap图所有物种存放位置")
    parser.add_argument("-s","--savepath",default = "Enrichment",help = "富集结果保存路径，默认Enrichment")
    parser.add_argument("-pp","--ppipath",default = "PPI_network",help = "PPI结果保存路径，默认PPI_network")
    parser.add_argument("-o","--org",default = "hsa", help = "kegg物种拉丁名缩写，默认hsa")
    parser.add_argument("-n","--number",default = 25, help = "ppi绘图蛋白top数量，默认25")
    parser.add_argument("-bf", "--backfile", help="background file",type = str)
    parser.add_argument("-y", "--ypt",default = "F", help="是否用于云平台",type = str)
    parser.add_argument("-font","--fontfamily",default = "sans", help = "富集绘图字体，默认sans")
    parser.add_argument("-c","--cores",default = 5, type= int,help = "线程数，默认5")
    parser.add_argument("-sm","--snakefile",default = "snakemake/Diffanalysis/diffgetdata/diff_enrich.smk",help = "snakefile路径")
    
    args = parser.parse_args()
    # protn = "富集蛋白列表.xlsx"
    # system("cp {inputfile} {file}".format(inputfile=args.diffile,file=protn))
    # if len(glob.glob("*-diff-*.xls")) == 0:
    
    if args.org == "updata":
      os.system(f"cp {args.backfile} ./")
      if os.path.exists("background.oecloud"):
        os.mkdir("./background/")
        os.system(f"tar -xf background.oecloud -C ./background")
      elif os.path.exists("data.oecloud"):
        os.system(f"unzip -q -o data.oecloud -d ./background/")
        os.system(f"cp -r ./background/data/* ./background/")
      args.keggdb="background"
      args.database="background"
      htmlfile=glob.glob("background/*.html")[0]
      htmlfile=htmlfile.rsplit("/",1)[1]
      args.org=re.match("^[a-z]{2,5}",htmlfile).group()
    else:
      if args.keggdb=="":
        if args.omic == "M":
          args.keggdb="/data/hstore4/database/kegg/kegg/"
        else:
          args.keggdb="/data/hstore3/database/kegg/"
    
    if args.enrichtype=="all":
      args.enrichtype="G,K,W,R,I,ppi"
    
    if args.omic == "M" and re.match(".*R.*",args.enrichtype) and args.org in["hsa","mmu","rno"]:
      # print("makediff -p {path} -f {file} -om {omic} -cb -o temp".format(path=args.inputpath,file=args.diffile,omic=args.omic))
      system("makediff -p {path} -f {file} -om {omic} -cb -o temp".format(path=args.inputpath,file=args.diffile,omic=args.omic))
    
    # print("makediff -p {path} -f {file} -om {omic} -o temp".format(path=args.inputpath,file=args.diffile,omic=args.omic))
    system("makediff -p {path} -f {file} -om {omic} -o temp".format(path=args.inputpath,file=args.diffile,omic=args.omic))
    
    def collect_reads():
      files = glob.glob("temp/*-diff-*.xls")
      files.sort()
      reads = [file.split('/')[-1].split('-diff-')[0] for file in files]
      other = files[0].split('/')[-1].split('-diff-')[1]
      return reads,other
  
    compare = list(set(collect_reads()[0]))
    
    # print(compare)
    sys.stderr = FilteredPrinter(filtered_print, sys.stderr)
    sys.stdout = FilteredPrinter(filtered_print2, sys.stdout)
    
    config = {"path":{"projectpath":args.savepath,
                      "projectppipath":args.ppipath,
                      "diff_ana_rds_savepath":"",
                      "diff_filter_rds_savepath":"",
                      "enrich_xls_savepath":"temp"},
              "params":{"omic":args.omic,
                        "all_compare":compare,
                        "enrichtype":args.enrichtype,
                        "enrich":{"org":args.org,
                                  "databasefrom":"/data/hstore3/database/background/" if args.omic=="P" else "/data/hstore4/database/kegg/kegg/",
                                  "database":args.database,
                                  "keggdb":args.keggdb,
                                  "ppitopn":args.number},
                        "plot":{"family":args.fontfamily}}}
    
    runsnakemake(snakefilepath = args.snakefile,
                 config = config,
                 configfiles =[],
                 cores = args.cores,
                 targets = ["enrich_result"])
    
    system("find -name .snakemake_timestamp -print0  | xargs -0  rm -rf")
    shutil.rmtree(".snakemake",ignore_errors=True)
    shutil.rmtree("temp",ignore_errors=True)
    
    if args.ypt=="T":
        os.system("chmod -R 777 ./* >/dev/null 2>&1")
        system("rm -rf background")
        system("rm -rf *.oecloud")
