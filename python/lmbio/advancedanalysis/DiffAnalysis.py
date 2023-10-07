#!/opt/conda/bin/python
import rpy2.robjects as robjects
from lmbio.advancedanalysis.organizedata import organizedata
from lmbio.advancedanalysis.fetch import getpredealparams
from lmbio.basic.printfilter import FilteredPrinter,filtered_print,filtered_print2
from lmbio.basic.runsnakemake import runsnakemake
from lmbio.basic.system import system
import argparse
import tarfile
import shutil
import glob
import sys
import os
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-om","--omic",default = "M", help = "项目类型,M(代谢),ML(脂质),P(蛋白),PF(工业蛋白)",required = True)
    parser.add_argument("-mt","--MStype",default = "auto", help = "项目子类别")
    parser.add_argument("-or","--org",default = "hsa", help = "物种简写")
    
    parser.add_argument("-rf","--rawfile",default = "数据矩阵.xlsx", help = "数据矩阵")
    parser.add_argument("-rdf","--rawdatafile",default = [""], help = "表达数据",nargs = "+")
    parser.add_argument("-rcf","--rawclassfile",default =  [""], help = "分组数据",nargs = "+")
    parser.add_argument("-rcpf","--rawcomparefile",default =  [""], help = "比较组数据",nargs = "+")
    parser.add_argument("-fp","--fillpalette",default =  "procol", help = "颜色")
    
    parser.add_argument("-tfp","--tempfilepath",default = "oecloud", help = "数据保存路径")
    # parser.add_argument("-df","--datafile",default = "oecloud/rawdata/datafile.txt", help = "数据保存路径")
    # parser.add_argument("-if","--infofile",default = "oecloud/rawdata/infofile.txt", help = "信息保存路径")
    # parser.add_argument("-cf","--classfile",default = "classfile.yaml", help = "分组保存路径")
    # parser.add_argument("-ctf","--classtypefile",default = "classtype.xlsx", help = "分组绘图保存路径")
    # parser.add_argument("-cpf","--comparefile",default = "compare.yaml", help = "比较组保存路径")
    
    parser.add_argument("-oe","--oecloud",default =  False, help = "是否已存在oecloud",action = 'store_true')
    
    # 数据预处理
    parser.add_argument("-p","--predeal",default = "T", help = "是否进行数据预处理")
    parser.add_argument("-mrp","--missvarremovepercent",default = -1, type= float,help = "缺失值筛选比例,默认0.5")
    parser.add_argument("-mrg","--missvarremovebygroup",default = "auto", help = "是否按组进行缺失值比例计算,使用参数将不按组进行计算")
    parser.add_argument("-mvf","--missvarfillmethod",default = "auto", 
                        help = "缺失值填充方式,包含none,mean_half,na_half,halfmin,valuemin,min,mean,median,knn_var,knn_smp,ppca,bpca,svdlmpute",
                        choices = ["none","mean_half","na_half","halfmin","valuemin","min","mean","median","knn_var","knn_smp","ppca","bpca","svdlmpute"])
    parser.add_argument("-rn","--rowNorm",default = "auto", 
                        help = "归一化方式,包含NULL,SumNorm,MedianMeanNorm,MedianNorm,QuantileNorm,SamplePQN(ref需提供样本名),GroupPQN(ref需提供组名),CompNorm(ref需提供蛋白或者代谢名),SpecNorm(ref需提供样本名)",
                        choices = ["NULL","SumNorm","MedianMeanNorm","MedianNorm","QuantileNorm","SamplePQN","GroupPQN","CompNorm"])
    parser.add_argument("-tn","--transNorm",default = "auto", 
                        help = "转化方式,包含NULL,LogNorm,Log2minNorm,SrNorm,CrNorm",
                        choices=["NULL","LogNorm","Log2minNorm","SrNorm","CrNorm"])
    parser.add_argument("-sn","--scaleNorm",default = "auto", 
                        help = "标准化方式,包含NULL,MeanCenter,AutoNorm,ParetoNorm,RangeNorm",
                        choices = ["NULL","MeanCenter","AutoNorm","ParetoNorm","RangeNorm"])
    parser.add_argument("-r","--ref",default = "auto", help = "归一化方式提供参考样本或者特征的参数")
    parser.add_argument("-fi","--filter",default = "auto", 
                        help = "数量筛选方式,默认为mean,包含rsd,nrsd,mean,sd,iqr,mad,median,none",
                        choices = ["rsd","nrsd","mean","sd","iqr","mad","median","none"])
    parser.add_argument("-re","--remainnum",default = -1, type= int,help = "特征筛选上限,默认为100000,输入0为自动筛选数量")
    parser.add_argument("-qf","--qcFilter",default = "auto",help = "是否进行QC样本的rsd筛选")
    parser.add_argument("-rs","--rsd",default = -1,type= int,help = "QC样本rsd筛选范围,默认30")
    
    # 多元统计分析参数
    parser.add_argument("-m", "--mulstatistics_mode",default = ["PCA","PLS-DA","OPLS-DA"],nargs = "+",
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
    
    # 差异分析参数
    parser.add_argument("-pm","--p_method",default = "t.test", help = "双组p值计算方法",choices = ["t.test","f.test","ft.test","wilcox.test","oneway.test","kruskal.test"])
    parser.add_argument("-mpm","--mul_p_method",default = "oneway.test", help = "三组及以上p值计算方法",choices = ["oneway.test","kruskal.test"])
    parser.add_argument("-pam","--p_adjust_method",default = "BH", help = "校正p值计算方法",choices = ["holm","hochberg","hommel","bonferroni","BH","BY","fdr"])
    parser.add_argument("-pd","--paired",default = "F", help = "是否配对")
    parser.add_argument("-l","--log",default = False, help = "数据是否经过log标准化",action='store_true')
    parser.add_argument("-ff","--fcfilter",default = [1.2,1.5,2], help = "FC筛选标准",type = float,nargs = "+")
    parser.add_argument("-fft","--fcfiltertype",default = "+-", help = "FC筛选标准")
    parser.add_argument("-eff","--errorptofcfilter",default = 2, help = "FC筛选标准",type = float)
    parser.add_argument("-vf","--vipfilter",default = 0, help = "vip筛选标准",type = float)
    parser.add_argument("-pf","--pfilter",default = 0.05, help = "p-value筛选标准",type = float)
    parser.add_argument("-apf","--adjpfilter",default = 0, help = "adjp-value筛选标准",type = float)
  
    # 绘图参数
    parser.add_argument("-it","--imagetype",default = ["png","pdf"],help = "图片格式",nargs = "+")
    parser.add_argument("-fa","--family",default = "sans",type=str, help = "字体")
    parser.add_argument("-wi","--width",default = 0, type= float,help = "图片宽度")
    parser.add_argument("-he","--height",default = 0, type= float,help = "图片高度")
    
    parser.add_argument("-bf","--backfile",default = "", type=str,help = "oecloud背景文件")
    
    # 运算相关参数
    parser.add_argument("-rp", "--report", help="是否出网页报告，默认T", nargs="?", type=str, default="T")
    parser.add_argument("-pr","--projectpath",default = "差异分析结果", type=str,help = "保存目录")
    parser.add_argument("-z","--zip",default = False, help = "是否进行压缩",action='store_true')
    parser.add_argument("-sm","--snakefile",default = "snakemake/Diffanalysis/whole/diff_analysis.smk",help = "snakefile路径")
    parser.add_argument("-c","--cores",default = 5, type= int,help = "线程数")
    parser.add_argument("-t","--targets",default = ["diff_getdata_result"],nargs="*",help = "运行规则",
                        choices = ["diff_getdata_result",
                                   "diff_getdata_result_test",
                                   "diff_result",
                                   "mul_result_all",
                                   "mul_result_score",
                                   "mul_result_score2"])

    args = parser.parse_args()
    args2 = vars(args)
    
    pretargets = ["protein_result","diff_getdata_result","diff_getdata_result_test"]
    judgetargets = False
    for targets in pretargets:
      if targets in args.targets:
        judgetargets = True
        break
    
    if len(args.fcfilter) > 1 and judgetargets:
      targets = ["diff_result"]
    else:
      targets = args.targets
    
    # print(targets)
    # print(args)
    args.datafile = args.tempfilepath+"/rawdata/datafile.txt"
    args.infofile = args.tempfilepath+"/rawdata/infofile.txt"
    args.classfile = args.tempfilepath+"/rawdata/classfile.yaml"
    args.classtypefile = args.tempfilepath+"/rawdata/classtype.xlsx"
    args.comparefile = args.tempfilepath+"/rawdata/compare.yaml"
    
    sys.stderr = FilteredPrinter(filtered_print, sys.stderr)
    sys.stdout = FilteredPrinter(filtered_print2, sys.stdout)
    
    args.log10L = getpredealparams(omic = args.omic,
                                   MStype = args.MStype,
                                   name = "log10L",
                                   value = args.log10L)
                                   
    if args.log10L == "T":
      if args.log :
        args.log10L = False
      else:  
        args.log10L = True
    else:
      args.log10L = False
    
    if args.predeal == "T" :
      datafile = args.tempfilepath+"/rawdata/rawdatafile.txt"
      for i in ["missvarremovepercent","missvarremovebygroup","missvarfillmethod",
                "rowNorm","transNorm","scaleNorm",
                "ref","filter","remainnum","qcFilter","rsd"]:
        args2[i] = getpredealparams(omic = args.omic,
                                    MStype = args.MStype,
                                    name = i,
                                    value = args2[i])
      
      if args.transNorm != "NULL" or args.scaleNorm != "NULL":
        args.log=True
        args.log10L=False
    else:
      datafile = args.datafile
    
    # print(args)
    
    if not args.oecloud :
      organizedata(rawfile = [args.rawfile],
                   rawdatafile = args.rawdatafile,
                   rawclassfile = args.rawclassfile,
                   rawcomparefile = args.rawcomparefile,
                   datafile = datafile,
                   infofile = args.infofile,
                   classfile = args.classfile,
                   classtypefile = args.classtypefile,
                   comparefile = args.comparefile,
                   fillpalette = args.fillpalette,
                   p_method = args.p_method,
                   mul_p_method = args.mul_p_method,
                   p_adjust_method = args.p_adjust_method,
                   paired = args.paired)
    if args.org == "updata":
      os.system(f"cp {args.backfile} ./")
      if os.path.exists("background.oecloud"):
        os.mkdir("./background/")
        os.system(f"tar -xf background.oecloud -C ./background")
      elif os.path.exists("data.oecloud"):
        os.system(f"unzip -q -o data.oecloud -d ./background/")
        os.system(f"cp -r ./background/data/* ./background/")
    
    if args.vipfilter > 0 :
      args.diff_ana_rds_savepath = args.tempfilepath+"/diff_ana_vip"
    else:
      args.diff_ana_rds_savepath = args.tempfilepath+"/diff_ana"

    args.diff_filter_rds_savepath = args.tempfilepath+"/diff_filter"
    args.config_savepath = args.tempfilepath+"/config"
    args.mulstatistics_rds_savepath = args.tempfilepath+"/mulstatisticsanalyst"
    args.enrich_xls_savepath = args.tempfilepath+"/enrich"

    config = {"path":{"tempfilepath":args.tempfilepath,
                      "rawfile":args.rawfile,
                      "rawdatafile":datafile,
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
              "params":{"omic":args.omic,
                        "MStype":args.MStype,
                        "predeal":{"processsavepath":args.tempfilepath+"/predeal/",
                                   "missvarremovepercent":args.missvarremovepercent,
                                   "missvarremovebygroup":args.missvarremovebygroup,
                                   "missvarfillmethod":args.missvarfillmethod,
                                   "rowNorm":args.rowNorm,
                                   "transNorm":args.transNorm,
                                   "scaleNorm":args.scaleNorm,
                                   "ref":args.ref,
                                   "filter":args.filter,
                                   "remainnum":args.remainnum,
                                   "qcFilter":args.qcFilter,
                                   "rsd":args.rsd},
                        "enrich":{"org":args.org,
                                  "databasefrom":"/data/hstore3/database/background/" if args.omic=="P" or args.omic=="PF" else "/data/hstore4/database/kegg/kegg/",
                                  "database":"background/",
                                  "keggdb":"/data/hstore3/database/kegg/" if args.omic=="P" or args.omic=="PF" else "/data/hstore4/database/kegg/kegg/"},
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

    if "OPLS-DA" in args.mulstatistics_mode and args.oplspermI > 0:
      config["params"].update({"plot_permutation":{"permutation_mode":"OPLS-DA"}})
    elif "PLS-DA" in args.mulstatistics_mode and args.plspermI > 0:
      config["params"].update({"plot_permutation":{"permutation_mode":"PLS-DA"}})
    else:
      config["params"].update({"plot_permutation":{"permutation_mode":"none"}})

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
      fcfilter = float(fcfilter)
      config["params"]["diff_filter"]["fcfilter"] = [fcfilter]
      runsnakemake(snakefilepath = args.snakefile,
                   config = config,
                   configfiles =[args.comparefile],
                   cores = args.cores,
                   targets = args.targets)
                   
    if args.report == "T":
      if args.omic == "P" or args.omic == "PF":
        system("Prorecloud")
        os.remove("analysis.yaml")
        os.remove("Report.py")
        shutil.rmtree("src",ignore_errors=True)
      elif args.omic == "M" or args.omic == "ML":
        if glob.glob("./"+args.projectpath+"/*/KEGG_map") :
          keggfile=glob.glob("./"+args.projectpath+"/*/KEGG_map")[0]
          system("zip -qr '"+keggfile+".zip' '"+keggfile+"'")
          system("rm -rf '"+keggfile+"'")
        system("tool_Report --type 全谱代谢-云平台 --savepath '"+args.projectpath+"'")
    
    os.system("chmod -R 777 ./* >/dev/null 2>&1")
    os.system("find -name .snakemake_timestamp -print0  | xargs -0  rm -rf")
    if args.zip :
      with tarfile.open(args.projectpath+".tar", "w") as tar:
        tar.add(args.projectpath, arcname=os.path.basename(args.projectpath))
      shutil.rmtree(".snakemake",ignore_errors=True)
      shutil.rmtree("background",ignore_errors=True)
      shutil.rmtree("background_out",ignore_errors=True)
      shutil.rmtree(args.tempfilepath,ignore_errors=True)
      shutil.rmtree(args.projectpath,ignore_errors=True)
      # os.remove("log.txt")

