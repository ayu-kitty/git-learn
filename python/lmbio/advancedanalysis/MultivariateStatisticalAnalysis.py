#!/opt/conda/bin/python
import rpy2.robjects as robjects
from lmbio.advancedanalysis.organizedata import organizedata
from lmbio.basic.printfilter import FilteredPrinter,filtered_print,filtered_print2
from lmbio.basic.runsnakemake import runsnakemake
import lmbio
import argparse
import os
import tarfile
import shutil
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-rf","--rawfile",default = "数据矩阵.xlsx", help = "数据矩阵")
    parser.add_argument("-df","--rawdatafile",default = [""], help = "表达数据",nargs = "+")
    parser.add_argument("-cf","--rawclassfile",default =  [""], help = "分组数据",nargs = "+")
    parser.add_argument("-cpff","--rawcomparefile",default =  [""], help = "比较组数据",nargs = "+")
    parser.add_argument("-fp","--fillpalette",default =  "procol",
                        choices = list(robjects.r("names(lmbio::palette_color)")), 
                        help = "颜色")
    
    # 数据预处理
    parser.add_argument("-p","--predeal",default = False, 
                        help = "是否进行数据预处理",action='store_true')
    parser.add_argument("-mrp","--missvarremovepercent",default = 0.5, type= float,
                        help = "缺失值筛选比例,默认0.5")
    parser.add_argument("-nmrg","--nmissvarremovebygroup",default = True, 
                        help = "是否按组进行缺失值比例计算,使用参数将不按组进行计算",action='store_false',
                        dest = "missvarremovebygroup")
    parser.add_argument("-mvf","--missvarfillmethod",default = "halfmin", 
                        help = "缺失值填充方式,包含none,mean_half,halfmin,valuemin,min,mean,median,knn_var,knn_smp,ppca,bpca,svdlmpute",
                        choices = ["none","mean_half","halfmin","valuemin","min","mean","median","knn_var","knn_smp","ppca","bpca","svdlmpute"])
    parser.add_argument("-rn","--rowNorm",default = "NULL", 
                        help = "归一化方式,包含NULL,MedianMeanNorm,SumNorm,MedianNorm,QuantileNorm,SamplePQN(ref需提供样本名),GroupPQN(ref需提供组名),CompNorm(ref需提供蛋白或者代谢名),SpecNorm(ref需提供样本名)",
                        choices = ["NULL","SumNorm","MedianMeanNorm","MedianNorm","QuantileNorm","SamplePQN","GroupPQN","CompNorm"])
    parser.add_argument("-tn","--transNorm",default = "NULL", 
                        help = "转化方式,包含LogNorm,SrNorm,CrNorm",
                        choices=["NULL","LogNorm","SrNorm","CrNorm"])
    parser.add_argument("-sn","--scaleNorm",default = "NULL", 
                        help = "标准化方式,包含MeanCenter,AutoNorm,ParetoNorm,RangeNorm",
                        choices = ["MeanCenter","AutoNorm","ParetoNorm","RangeNorm"])
    parser.add_argument("-r","--ref",default = "NULL", 
                        help = "归一化方式提供参考样本或者特征的参数")
    parser.add_argument("-fi","--filter",default = "mean", 
                        help = "数量筛选方式,默认为mean,包含rsd,nrsd,mean,sd,iqr,mad,median,none",
                        choices = ["rsd","nrsd","mean","sd","iqr","mad","median","none"])
    parser.add_argument("-re","--remainnum",default = 100000, type= int,
                        help = "特征筛选上限,默认为100000,输入0为自动筛选数量")
    parser.add_argument("-nqf","--nqcFilter",default = True,action = "store_false", 
                        help = "是否进行QC样本的rsd筛选",dest = "qcFilter")
    parser.add_argument("-rs","--rsd",default = 30,type= int, 
                        help = "QC样本rsd筛选范围,默认30")
    
     # 多元统计分析参数
    parser.add_argument("-m", "--mulstatistics_mode",default = ["PCA","PLS-DA","OPLS-DA"],nargs = "+",
                        choices = ["PCA","PLS-DA","OPLS-DA"],help = "多元统计分析模式")
    parser.add_argument("-ma", "--mulstatistics_all",default = False,
                        help = "是否进行所有样本分析",action='store_true')
    parser.add_argument("-mq", "--mulstatistics_all_qc",default = False,
                        help = "是否进行除QC外所有样本分析",action='store_true')
    parser.add_argument("-l","--log10L",default = False, 
                        help = "多元统计是否log处理",action='store_true')
    parser.add_argument("-pp","--plspermI",default = 0, type= int,help = "PLS-DA响应检测排序次数")
    parser.add_argument("-op","--oplspermI",default = 200, type= int,help = "OPLS-DA响应检测排序次数")
    parser.add_argument("-pcs","--pcascaleC",default = "standard",type=str, help = "PCA多元统计归一化模式",
                        choices = ["standard","pareto","none","center"])
    parser.add_argument("-pls","--plsscaleC",default = "pareto",type=str, help = "PLS-DA多元统计归一化模式",
                        choices = ["standard","pareto","none","center"])
    parser.add_argument("-ops","--oplsscaleC",default = "pareto",type=str, help = "OPLS-DA多元统计归一化模式",
                        choices = ["standard","pareto","none","center"])
                        
    # 绘图参数
    parser.add_argument("-fa","--family",default = "sans",type=str, help = "字体")
    parser.add_argument("-wi","--width",default = 0, type= float,help = "图片宽度")
    parser.add_argument("-he","--height",default = 0, type= float,help = "图片高度")
    
    # 运算相关参数
    parser.add_argument("-pr","--projectpath",default = "多元统计分析", type=str,help = "保存目录")
    parser.add_argument("-z","--zip",default = False, help = "是否进行压缩",action='store_true')
    parser.add_argument("-sm","--snakefile",default = "snakemake/MultivariateStatisticalAnalysis-2023/whole/multivariate.smk",help = "snakefile路径")
    parser.add_argument("-c","--cores",default = 5, type= int,help = "线程数")
    parser.add_argument("-t","--targets",default = ["mul_result_all"],nargs="*",help = "运行规则",
                        choices = ["mul_result_all",
                                   "mul_result_score",
                                   "mul_result_score2"])

    args = parser.parse_args()
    
    # print(args)
    sys.stderr = FilteredPrinter(filtered_print, sys.stderr)
    sys.stdout = FilteredPrinter(filtered_print2, sys.stdout)
    
    if args.predeal :
      datafile = "oecloud/rawdata/rawdatafile.txt"
    else:
      datafile = "oecloud/rawdata/datafile.txt"
    
    args.tempfilepath="oecloud"
    args.datafile = args.tempfilepath+"/rawdata/datafile.txt"
    args.infofile = args.tempfilepath+"/rawdata/infofile.txt"
    args.classfile = args.tempfilepath+"/rawdata/classfile.yaml"
    args.classtypefile = args.tempfilepath+"/rawdata/classtype.xlsx"
    args.comparefile = args.tempfilepath+"/rawdata/compare.yaml"
    
    organizedata(rawfile = [args.rawfile],
                 rawdatafile = args.rawdatafile,
                 rawclassfile = args.rawclassfile,
                 rawcomparefile = args.rawcomparefile,
                 fillpalette = args.fillpalette,
                 datafile = datafile,
                 infofile = args.infofile,
                 classfile = args.classfile,
                 classtypefile = args.classtypefile,
                 comparefile = args.comparefile)
    
    config = {"path":{"rawfile":args.rawfile,
                      "rawdatafile":datafile,
                      "datafile":args.datafile,
                      "infofile":args.infofile,
                      "classfile":args.classfile,
                      "classtypefile":args.classtypefile,
                      "comparefile":args.comparefile,
                      "projectpath":args.projectpath},
              "params":{"predeal":{"processsavepath":"oecloud/predeal/",
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
                        "mulstatistics_all":args.mulstatistics_all,
                        "mulstatistics_all-qc":args.mulstatistics_all_qc,
                        "mulstatistics_mode":args.mulstatistics_mode,
                        "PCA":{"log10L":args.log10L,
                               "scaleC":args.pcascaleC},
                        "PLS-DA":{"log10L":args.log10L,
                                  "permI":args.plspermI,
                                  "scaleC":args.plsscaleC},
                        "OPLS-DA":{"log10L":args.log10L,
                                   "permI":args.oplspermI,
                                   "scaleC":args.oplsscaleC},
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
    
    runsnakemake(snakefilepath = args.snakefile,
                 config = config,
                 configfiles =["oecloud/rawdata/compare.yaml"],
                 cores = args.cores,
                 targets = args.targets)
    
    os.system("chmod -R 777 ./*")
    if args.zip :
      with tarfile.open(args.projectpath+".tar", "w") as tar:
        tar.add(args.projectpath, arcname=os.path.basename(args.projectpath))
      # os.remove("log.txt")
      shutil.rmtree(".snakemake")
      shutil.rmtree("oecloud")
      shutil.rmtree(args.projectpath)
