#!/opt/conda/bin/python
from lmbio.advancedanalysis.organizedata import organizedata
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
    parser.add_argument("-rf","--rawfile",default = "数据矩阵.xlsx", help = "数据矩阵")
    parser.add_argument("-rdf","--rawdatafile",default = [""], help = "表达数据",nargs = "+")
    parser.add_argument("-rcf","--rawclassfile",default =  [""], help = "分组数据",nargs = "+")
    parser.add_argument("-rcpf","--rawcomparefile",default =  [""], help = "比较组数据",nargs = "+")
    parser.add_argument("-tfp","--tempfilepath",default = "oecloud", help = "数据保存路径")
    parser.add_argument("-min","--gseamin",default = 0, type= int, help = "基因集最小数量，代谢默认5，蛋白默认15")
    parser.add_argument("-graph","--gseagraph",default = 20, type= int, help = "绘制激活+抑制总图数，默认20")
    parser.add_argument("-lt","--listfile",default = "", help = "指定绘图的通路，传入文件")
    
    parser.add_argument("-om", "--omic",default = "P",help = "组学分类，蛋白为P，代谢为M，默认P")
    parser.add_argument("-t", "--gseatype",default = "all",help = "富集数据库类型，默认全选，单项示例G,K")
    parser.add_argument("-f","--file",default = "数据矩阵.xlsx", help = "需要做gsea表格名称，默认为数据矩阵.xlsx")
    parser.add_argument("-d","--database",default = "background/", help = "富集背景文件路径，默认background/")
    parser.add_argument("-s","--savepath",default = "result",help = "富集结果保存路径，默认result")
    parser.add_argument("-o","--org",default = "hsa", help = "kegg物种拉丁名缩写，默认hsa")
    parser.add_argument("-bf", "--backfile", help="background file",type = str)
    parser.add_argument("-y", "--ypt",default = "F", help="是否用于云平台",type = str)
    parser.add_argument("-font","--fontfamily",default = "sans", help = "富集绘图字体，默认sans")
    parser.add_argument("-c","--cores",default = 5, type= int,help = "线程数，默认5")
    parser.add_argument("-sm","--snakefile",default = "snakemake/Diffanalysis/diffgetdata/diff_gsea.smk",help = "snakefile路径")
    
    args = parser.parse_args()
    
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
      if args.omic == "M" or args.omic == "ML":
        args.keggdb="/data/hstore4/database/kegg/kegg/"
      else:
        args.keggdb="/data/hstore3/database/kegg/"
    
    if args.gseatype=="all":
      args.gseatype="G,K"
    if args.listfile!="":
      args.listfile = " -lt "+args.listfile

    args.inputfile = args.tempfilepath+"/rawdata/inputfile.txt"
    args.datafile = args.tempfilepath+"/rawdata/datafile.txt"
    args.infofile = args.tempfilepath+"/rawdata/infofile.txt"
    args.classfile = args.tempfilepath+"/rawdata/classfile.yaml"
    args.classtypefile = args.tempfilepath+"/rawdata/classtype.xlsx"
    args.comparefile = args.tempfilepath+"/rawdata/compare.yaml"
    
    organizedata(rawfile = [args.rawfile],
                 rawdatafile = args.rawdatafile,
                 rawclassfile = args.rawclassfile,
                 rawcomparefile = args.rawcomparefile,
                 datafile = args.datafile,
                 infofile = args.infofile,
                 classfile = args.classfile,
                 classtypefile = args.classtypefile,
                 comparefile = args.comparefile)
    if args.omic=="M" or args.omic=="ML":
        getinfo = args.infofile
    else:
        getinfo = "F"
    # print(compare)
    sys.stderr = FilteredPrinter(filtered_print, sys.stderr)
    sys.stdout = FilteredPrinter(filtered_print2, sys.stdout)
    config = {"path":{"projectpath":args.savepath,
                      "datafile":args.datafile,
                      "infofile":args.infofile,
                      "classfile":args.classfile,
                      "classtypefile":args.classtypefile,
                      "comparefile":args.comparefile,
                      "diff_ana_rds_savepath":"",
                      "diff_filter_rds_savepath":"",
                      "enrich_xls_savepath":"temp"},
              "params":{"omic":args.omic,
                        "gseatype":args.gseatype,
                        "enrich":{"org":args.org,
                                  "databasefrom":"/data/hstore3/database/background/" if args.omic=="P" else "/data/hstore4/database/kegg/kegg/",
                                  "database":args.database,
                                  "keggdb":args.keggdb,
                                  "infofile":getinfo,
                                  "gseamin":args.gseamin,
                                  "graph":args.gseagraph,
                                  "termlist":args.listfile,
                                  "gseadatapath":args.inputfile},
                        "plot":{"family":args.fontfamily}}}
    
    runsnakemake(snakefilepath = args.snakefile,
                 config = config,
                 configfiles =[args.comparefile],
                 cores = args.cores,
                 targets = ["gsea_result"])
    
    system("find -name .snakemake_timestamp -print0  | xargs -0  rm -rf")
    shutil.rmtree(".snakemake",ignore_errors=True)
    shutil.rmtree("temp",ignore_errors=True)
    
    if args.ypt=="T":
        os.system("chmod -R 777 ./* >/dev/null 2>&1")
        system("rm -rf background")
        shutil.rmtree(args.tempfilepath,ignore_errors=True)
