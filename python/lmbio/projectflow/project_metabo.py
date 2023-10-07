#!/opt/conda/bin/python

import rpy2.robjects as robjects
from lmbio.basic.printfilter import FilteredPrinter,filtered_print,filtered_print2
from lmbio.basic.runsnakemake import runsnakemake
from lmbio.basic.system import system
from lmbio.advancedanalysis.fetch import metafetch
import argparse
import tarfile
import shutil
import wget
import time
import glob
import sys
import os
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #一键式流程相关参数
    parser.add_argument("-om","--omic",default = "M", help = "项目类型,M/P")
    parser.add_argument("-or","--org",default = "", help = "物种简写")
    parser.add_argument("-rd","--rawdatapath",default = "projectdata/", help = "初始数据保存路径")
    parser.add_argument("-tfp","--tempfilepath",default = "rawdata/oecloud", help = "中间数据保存路径")
    parser.add_argument("-fp","--fillpalette",default =  "procol", help = "颜色")
    parser.add_argument("-f", "--force",default="F",help="强制重新分析")
                        
    # 差异分析参数
    parser.add_argument("-pm","--p_method",default = "t.test", help = "双组p值计算方法",choices = ["t.test","f.test","ft.test","wilcox.test","oneway.test","kruskal.test"])
    parser.add_argument("-mpm","--mul_p_method",default = "oneway.test", help = "三组及以上p值计算方法",choices = ["oneway.test","kruskal.test"])
    parser.add_argument("-pam","--p_adjust_method",default = "BH", help = "校正p值计算方法",choices = ["holm","hochberg","hommel","bonferroni","BH","BY","fdr"])
    parser.add_argument("-pd","--paired",default = "F", help = "是否配对")
    parser.add_argument("-l","--log",default = False, help = "数据是否经过log标准化",action='store_true',dest = "log")
    parser.add_argument("-ff","--fcfilter",default = [0], help = "FC筛选标准",type = float,nargs = "+")
    parser.add_argument("-fft","--fcfiltertype",default = "+-", help = "FC筛选标准")
    parser.add_argument("-eff","--errorptofcfilter",default = 2, help = "FC筛选标准",type = float)
    parser.add_argument("-vf","--vipfilter",default = -1, help = "vip筛选标准",type = float)
    parser.add_argument("-pf","--pfilter",default = -1, help = "p-value筛选标准",type = float)
    parser.add_argument("-apf","--adjpfilter",default = -1, help = "adjp-value筛选标准",type = float)
    
    # 多元统计分析参数
    parser.add_argument("-m", "--mulstatistics_mode",default = ["PCA","OPLS-DA"],nargs = "+",
                        choices = ["PCA","PLS-DA","OPLS-DA"],help = "多元统计分析模式")
    parser.add_argument("-lo","--log10L",default = "auto", help = "多元统计是否log处理")
    parser.add_argument("-pp","--plspermI",default = 0, type= int,help = "PLS-DA响应检测排序次数")
    parser.add_argument("-op","--oplspermI",default = 200, type= int,help = "OPLS-DA响应检测排序次数")
    parser.add_argument("-pcs","--pcascaleC",default = "standard",type=str, help = "PCA多元统计归一化模式",
                        choices = ["standard","pareto","none","center"])
    parser.add_argument("-pls","--plsscaleC",default = "pareto",type=str, help = "PLS-DA多元统计归一化模式",
                        choices = ["standard","pareto","none","center"])
    parser.add_argument("-ops","--oplsscaleC",default = "pareto",type=str, help = "OPLS-DA多元统计归一化模式",
                        choices = ["standard","pareto","none","center"])
  
    # 绘图参数
    parser.add_argument("-it","--imagetype",default = ["png","pdf"],help = "图片格式",nargs = "+")
    parser.add_argument("-fa","--family",default = "sans",type=str, help = "字体")
    parser.add_argument("-wi","--width",default = 0, type= float,help = "图片宽度")
    parser.add_argument("-he","--height",default = 0, type= float,help = "图片高度")
    
    # 运算相关参数
    parser.add_argument("-pr","--projectpath",default = "", type=str,help = "保存目录")
    parser.add_argument("-z","--zip",default = False, help = "是否进行压缩",action='store_true')
    parser.add_argument("-sm","--snakefile",default = "snakemake/MetaboAnalysis-2023/metaboanalysis.smk",help = "snakefile路径")
    parser.add_argument("-c","--cores",default = 5, type= int,help = "线程数")
    parser.add_argument("-t","--targets",default = ["metabo_result"],nargs="*",help = "运行规则",
                        choices = ["metabo_result",
                                   "metabo_result_pre",# 预处理给给实验
                                   "metabo_result_qc",# 质控
                                   "metabo_result_mul",# 多元统计
                                   "metabo_result_oecloud",# oecloud文件
                                   "metabo_result_lipid",# 差异脂质绘图
                                   "metabo_result_diff_xlsx",# 差异表格
                                   "metabo_result_diff_map",#差异绘图
                                   "metabo_result_diff_function",# 差异功能（gsea，富集）
                                   "metabo_result_diff_enrich",# 差异富集
                                   "metabo_result_diff_gsea",# 差异gsea
                                   "metabo_result_diff",# 整套差异
                                   "diff_getdata_result"]) # 云平台结果
    # 生成报告参数
    parser.add_argument("-sh","--sh", default = False, action = "store_true", help = "是否为售后，默认F")
    parser.add_argument("-r", "--report", help="是否出网页报告，默认T", nargs="?", type=str, default="T")
    #parser.add_argument("-ch","--check",default = F, help = "检查")
    parser.add_argument("-cl","--CLOUD",default = "T", help = "云交付")
    parser.add_argument("-ob","--moveobs",default = "T", help = "上传obs")

    args = parser.parse_args()
    
    # args.tempfilepath = "raw/oecloud"
    
    if args.force == "T":
      for file in ["raw.RData","classtype.xlsx","classfile.yaml"]:
        if os.path.exists(file):
          os.remove(file)
      shutil.rmtree(".snakemake",ignore_errors=True)
      shutil.rmtree("rawdata",ignore_errors=True)
      shutil.rmtree(args.tempfilepath,ignore_errors=True)
    
    pretargets = ["metabo_result","diff_getdata_result","diff_getdata_result_test"]
    judgetargets = False
    for targets in pretargets:
      if targets in args.targets:
        judgetargets = True
        break
    
    if len(args.fcfilter) > 1 and judgetargets:
      targets = ["diff_result"]
    else:
      targets = args.targets
    
    args.datafile = args.tempfilepath+"/rawdata/datafile.txt"
    args.infofile = args.tempfilepath+"/rawdata/infofile.txt"
    args.classfile = "classfile.yaml"
    args.classtypefile = "classtype.xlsx"
    args.comparefile = args.tempfilepath+"/rawdata/compare.yaml"
    
    projecttype = metafetch(savepath = args.tempfilepath,
                            datafile = args.datafile,
                            infofile = args.infofile,
                            classfile = args.classfile,
                            classtypefile = args.classtypefile,
                            comparefile = args.comparefile,
                            fillpalette = args.fillpalette,
                            p_method = args.p_method,
                            mul_p_method = args.mul_p_method,
                            p_adjust_method = args.p_adjust_method,
                            paired = args.paired)
    
    if args.log10L == "auto":
      if re.match("^Untarget.*Gcms",str(projecttype.rx2("class")[0])) or re.match(".*GCMS",str(projecttype.rx2("pro_type")[0])) or re.match(".*蜡质",str(projecttype.rx2("pro_type")[0])) or re.match(".*顶空",str(projecttype.rx2("pro_type")[0])):
        args.log10L = True
      else:
        args.log10L = False
    elif args.log10L == "T":
      args.log10L = True
    elif args.log10L == "F":
      args.log10L = False
    else:
      raise ValueError("--log10L请输入auto、T、F")
    
    if str(projecttype.rx2("class")[0]) == "Target":
      if args.adjpfilter < 0:
        args.adjpfilter = 0
      if args.vipfilter < 0:
        args.vipfilter = 0
      if (str(projecttype.rx2("pro_type")[0]) not in ["有机酸","氨基酸","黄酮酚类"]):
        if args.pfilter < 0:
          args.pfilter = 0
      else:
        args.pfilter = 0.05
    else:
      if args.vipfilter < 0:
        args.vipfilter = 1
      if args.pfilter < 0:
        args.pfilter = 0.05
      if args.adjpfilter < 0:
        args.adjpfilter = 0
    
    # print(projecttype)
    sys.stderr = FilteredPrinter(filtered_print, sys.stderr)
    sys.stdout = FilteredPrinter(filtered_print2, sys.stdout)
    
    # if args.vipfilter > 0 :
    #   args.diff_ana_rds_savepath = args.tempfilepath+"/diff_ana_vip"
    # else:
    #   args.diff_ana_rds_savepath = args.tempfilepath+"/diff_ana"
    args.diff_ana_rds_savepath = args.tempfilepath+"/diff_ana_vip"
    
    args.diff_filter_rds_savepath = args.tempfilepath+"/diff_filter"
    args.config_savepath = args.tempfilepath+"/config"  
    args.mulstatistics_rds_savepath = args.tempfilepath+"/mulstatisticsanalyst" 
    args.enrich_xls_savepath = args.tempfilepath+"/enrich" 
    
    if args.projectpath == "":
      if os.path.exists(args.tempfilepath+"/project_title.txt"):
        f = open(args.tempfilepath+"/project_title.txt")
        project_title = f.readline()
        project_title = project_title.strip()
        args.projectpath = project_title
        f.close()
      else:
        project_title = "report"
    
    if args.org == "":
      if os.path.exists(args.tempfilepath+"/org.txt"):
        f = open(args.tempfilepath+"/org.txt")
        org = f.readline()
        org = org.strip()
        args.org = org
        f.close()
      else:
        raise ValueError("请输入物种")
    
    if args.force == "T":
      print(args.projectpath)
      shutil.rmtree(args.projectpath,ignore_errors=True)
    
    config = {"path":{"tempfilepath":args.tempfilepath,
                      "rawdatapath":"rawdata/",
                      "rawdatafile":args.tempfilepath+"/rawdata/rawdatafile.txt",
                      "datafile":args.datafile,
                      "infofile":args.infofile,
                      "classfile":args.classfile,
                      "classtypefile":args.classtypefile,
                      "comparefile":args.comparefile,
                      "projectpath":args.projectpath,
                      "diff_ana_rds_savepath":args.diff_ana_rds_savepath,
                      "diff_filter_rds_savepath":args.diff_filter_rds_savepath,
                      "diff_xls":args.tempfilepath+"/diff",
                      "config_savepath":args.config_savepath,
                      "mulstatistics_rds_savepath":args.mulstatistics_rds_savepath,
                      "enrich_xls_savepath":args.enrich_xls_savepath},
              "params":{"projecttype":{"class":str(projecttype.rx2("class")[0]),
                                       "pro_type":str(projecttype.rx2("pro_type")[0]),
                                       "qc":bool(projecttype.rx2("qc")[0]),
                                       "compare":bool(projecttype.rx2("compare")[0])},
                        "omic":args.omic,
                        "predeal":{"processsavepath":args.tempfilepath+"/predeal/",
                                   "missvarremovepercent":0.5,
                                   "missvarremovebygroup":True,
                                   "missvarfillmethod":"halfmin",
                                   "rowNorm":"NULL",
                                   "transNorm":"NULL",
                                   "scaleNorm":"NULL",
                                   "ref":"NULL",
                                   "filter":"mean",
                                   "remainnum":100000,
                                   "qcFilter":True,
                                   "rsd":30},
                        "enrich":{"org":args.org,
                                  "databasefrom":"/data/hstore3/database/background/" if  args.omic=="P" or args.omic=="PF" else "/data/hstore4/database/kegg/kegg/",
                                  "database":"background/",
                                  "keggdb":"/data/hstore3/database/kegg/" if  args.omic=="P" or args.omic=="PF" else "/data/hstore4/database/kegg/kegg/"},
                        "diff_filter":{"log":args.log,
                                       "fcfilter":args.fcfilter,
                                       "fcfiltertype":args.fcfiltertype,
                                       "errorptofcfilter":args.errorptofcfilter,
                                       "vipfilter":args.vipfilter,
                                       "pfilter":args.pfilter,
                                       "adjpfilter":args.adjpfilter},
                        "mulstatistics_mode":args.mulstatistics_mode,
                        "PCA":{"log10L":args.log10L,
                               "scaleC":args.pcascaleC},
                        "PLS-DA":{"log10L":args.log10L,
                                  "permI":args.plspermI,
                                  "scaleC":args.plsscaleC},
                        "OPLS-DA":{"log10L":args.log10L,
                                   "permI":args.oplspermI,
                                   "scaleC":args.oplsscaleC},
                        "imagetype":args.imagetype,
                        "plot":{"family":args.family,
                                "width":args.width,
                                "height":args.height}}}
                                   
    # print(config)
    
    if (not os.path.exists(args.diff_ana_rds_savepath+"/fitfc.txt")) or len(config["params"]["diff_filter"]["fcfilter"]) == 1 or not judgetargets:
      runsnakemake(snakefilepath = args.snakefile,
                   config = config,
                   configfiles =[args.comparefile],
                   cores = args.cores,
                   targets = targets)
    
    if len(config["params"]["diff_filter"]["fcfilter"]) > 1 and judgetargets:
      f = open(args.diff_ana_rds_savepath+"/fitfc.txt")
      fcfilter = f.readline()
      f.close()
      fcfilter = float(fcfilter)
      config["params"]["diff_filter"]["fcfilter"] = [fcfilter]
      runsnakemake(snakefilepath = args.snakefile,
                   config = config,
                   configfiles =[args.comparefile],
                   cores = args.cores,
                   targets = args.targets)
    
    projectpath = args.projectpath
    if glob.glob("./"+projectpath+"/*/KEGG_map") :
      keggfile=glob.glob("./"+projectpath+"/*/KEGG_map")[0]
      system("zip -qr '"+keggfile+".zip' '"+keggfile+"'")
      system("rm -rf '"+keggfile+"'")
    if args.sh :
      args.moveobs = "F"
      args.CLOUD = "F"
    if args.report == "T" and args.targets[0] != "metabo_result_pre":
      os.system("find -name .snakemake_timestamp -print0  | xargs -0  rm -rf")
      shutil.rmtree(".snakemake")
      system("mkreport_oebio2 -cl " + args.CLOUD)
    os.system("chmod -R 777 ./* >/dev/null 2>&1")
    
    if args.moveobs == "T" and args.targets[0] != "metabo_result_pre" and args.CLOUD == "T":
      system("movealldata")

    if args.zip :
      project_title = args.projectpath
      with tarfile.open(project_title+".tar", "w") as tar:
        tar.add(args.projectpath, arcname=os.path.basename(args.projectpath))
      shutil.rmtree(".snakemake",ignore_errors=True)
      shutil.rmtree("rawdata",ignore_errors=True)
      shutil.rmtree("background",ignore_errors=True)
      shutil.rmtree("background_out",ignore_errors=True)
      shutil.rmtree(args.tempfilepath,ignore_errors=True)
      shutil.rmtree(args.projectpath,ignore_errors=True)
      
    print("oeInfo:运行完成")
