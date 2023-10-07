#!/opt/conda/bin/python
import rpy2.robjects as robjects
from lmbio.basic.printfilter import FilteredPrinter,filtered_print,filtered_print2
from lmbio.basic.runsnakemake import runsnakemake
from lmbio.basic.system import system
from lmbio.advancedanalysis.fetch import profetch
from lmbio.advancedanalysis.fetch import getpredealparams
import argparse
import wget
import os
import time
import glob
import tarfile
import zipfile
import shutil
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #一键式流程相关参数
    parser.add_argument("-l", "--link", help="分析单下载链接", nargs="?", type=str, default='None')
    parser.add_argument("-nd", "--nasdir", help="NAS盘rawdata路径，没有下载链接必须填写", nargs="?", type=str, default='./')
    parser.add_argument("-fd", "--fdir", help="运行流程路径，如/data/hstore3/Projects/202306，没有下载链接必须填写", nargs="?", type=str, default='./')
    parser.add_argument("-i", "--proid", help="项目编号，没有下载链接必须填写", nargs="?", type=str, default="")
    parser.add_argument("-b", "--background", help="背景文件夹路径，如/data/hstore3/shouhou/2023/202303/DLM202211261-售后1/background", nargs="?", type=str, default='None')
    parser.add_argument("-r", "--report", help="是否出网页报告，默认T", nargs="?", type=str, default="T")
    parser.add_argument("-m", "--mode", help="网页报告模板类型，默认自动选择", nargs="?", type=str, default="auto")
    parser.add_argument("-cl", "--cloud", help="是否云交付，默认F", nargs="?", type=str, default="F")
    parser.add_argument("-om","--omic",default = "P", help = "项目类型,P、PF")
    parser.add_argument("-or","--org",default = "hsa", help = "物种简写",required = True)
    parser.add_argument("-fet","--fetch",default = "T", help = "原始数据是否做下机数据标准提取，售后如剔除样本选F")
    parser.add_argument("-rd","--rawdatapath",default = "projectdata/", help = "初始数据保存路径")
    # parser.add_argument("-tfp","--tempfilepath",default = "rawdata/oecloud", help = "中间数据保存路径")
    parser.add_argument("-fp","--fillpalette",default =  "procol", help = "颜色")
                        
    # 差异分析参数
    parser.add_argument("-pm","--p_method",default = "t.test", help = "双组p值计算方法",choices = ["t.test","f.test","ft.test","wilcox.test","oneway.test","kruskal.test"])
    parser.add_argument("-mpm","--mul_p_method",default = "oneway.test", help = "三组及以上p值计算方法",choices = ["oneway.test","kruskal.test"])
    parser.add_argument("-pam","--p_adjust_method",default = "BH", help = "校正p值计算方法",choices = ["holm","hochberg","hommel","bonferroni","BH","BY","fdr"])
    parser.add_argument("-pd","--paired",default = "F", help = "是否配对")
    parser.add_argument("-lo","--log",default = False, help = "数据是否经过log标准化",action='store_true')
    parser.add_argument("-ff","--fcfilter",default = [1.2,1.5,2], help = "FC筛选标准",type = float,nargs = "+")
    parser.add_argument("-fft","--fcfiltertype",default = "+-", help = "FC筛选标准")
    parser.add_argument("-eff","--errorptofcfilter",default = 2, help = "FC筛选标准",type = float)
    parser.add_argument("-vf","--vipfilter",default = 0, help = "vip筛选标准",type = float)
    parser.add_argument("-pf","--pfilter",default = 0.05, help = "p-value筛选标准",type = float)
    parser.add_argument("-apf","--adjpfilter",default = 0, help = "adjp-value筛选标准",type = float)
  
    # 数据预处理
    parser.add_argument("-mrp","--missvarremovepercent",default = -1, type= float,help = "缺失值筛选比例,默认0.5")
    parser.add_argument("-mrg","--missvarremovebygroup",default = "auto", help = "是否按组进行缺失值比例计算,使用参数将不按组进行计算")
    parser.add_argument("-mvf","--missvarfillmethod",default = "auto", 
                        help = "缺失值填充方式,包含none,mean_half,na_half,halfmin,valuemin,min,mean,median,knn_var,knn_smp,ppca,bpca,svdlmpute",
                        choices = ["none","mean_half","na_half","halfmin","valuemin","min","mean","median","knn_var","knn_smp","ppca","bpca","svdlmpute"])
    parser.add_argument("-rn","--rowNorm",default = "auto", 
                        help = "归一化方式,包含NULL,SumNorm,MedianMeanNorm,MedianNorm,QuantileNorm,SamplePQN(ref需提供样本名),GroupPQN(ref需提供组名),CompNorm(ref需提供蛋白或者代谢名),SpecNorm(ref需提供样本名)",
                        choices = ["NULL","SumNorm","MedianMeanNorm","MedianNorm","QuantileNorm","SamplePQN","GroupPQN","CompNorm"])
    parser.add_argument("-tn","--transNorm",default = "auto", 
                        help = "转化方式,包含LogNorm,Log2minNorm,SrNorm,CrNorm",
                        choices=["NULL","LogNorm","Log2minNorm","SrNorm","CrNorm"])
    parser.add_argument("-sn","--scaleNorm",default = "auto", 
                        help = "标准化方式,包含NULL,MeanCenter,AutoNorm,ParetoNorm,RangeNorm",
                        choices = ["NULL","MeanCenter","AutoNorm","ParetoNorm","RangeNorm"])
    parser.add_argument("-rf","--ref",default = "auto", help = "归一化方式提供参考样本或者特征的参数")
    parser.add_argument("-fi","--filter",default = "auto", 
                        help = "数量筛选方式,默认为mean,包含rsd,nrsd,mean,sd,iqr,mad,median,none",
                        choices = ["rsd","nrsd","mean","sd","iqr","mad","median","none"])
    parser.add_argument("-re","--remainnum",default = -1, type= int,help = "特征筛选上限,默认为100000,输入0为自动筛选数量")
    parser.add_argument("-qf","--qcFilter",default = "auto",help = "是否进行QC样本的rsd筛选")
    parser.add_argument("-rs","--rsd",default = -1,type= int,help = "QC样本rsd筛选范围,默认30")
  
    # 绘图参数
    parser.add_argument("-it","--imagetype",default = ["png","pdf"],help = "图片格式",nargs = "+")
    parser.add_argument("-fa","--family",default = "sans",type=str, help = "字体")
    parser.add_argument("-wi","--width",default = 0, type= float,help = "图片宽度")
    parser.add_argument("-he","--height",default = 0, type= float,help = "图片高度")
    parser.add_argument("-voln","--volcanonum",default = 5,type= int,help = "火山图标注数量")
    
    # 运算相关参数
    parser.add_argument("-pr","--projectpath",default = "report", type=str,help = "保存目录")
    parser.add_argument("-z","--zip",default = False, help = "是否进行压缩",action='store_true')
    parser.add_argument("-sm","--snakefile",default = "snakemake/ProteinAnalysis-2023/proteinanalysis.smk",help = "snakefile路径")
    parser.add_argument("-c","--cores",default = 5, type= int,help = "线程数")
    parser.add_argument("-t","--targets",default = ["protein_result"],nargs="*",help = "运行规则",
                        choices = ["protein_result",
                                   "protein_result_pca",# pca分析
                                   "protein_result_oecloud",# oecloud文件
                                   "protein_result_qc",# 质控与统计
                                   "protein_result_diff_cor",# 差异相关性分析
                                   "protein_result_diff_motif",# 差异motif分析
                                   "protein_result_diff_function",# 差异功能（ppi，gsea，富集）
                                   "protein_result_diff_enrich",# 差异富集
                                   "protein_result_diff_gsea",# 差异gsea
                                   "protein_result_diff_ppi",# 差异ppi
                                   "protein_result_diff_xlsx",# 差异表格
                                   "protein_result_diff_map",#差异绘图
                                   "protein_result_diff",# 整套差异
                                   "diff_getdata_result"])

    args = parser.parse_args()
    args2 = vars(args)
    
    args.tempfilepath = "rawdata/oecloud"
    if args.link != "None":
        loctime = time.strftime("%Y%m", time.localtime())
        prodir="/data/hstore3/localPro/1.Projects/"+ loctime #NAS路径
        rundir="/data/hstore3/Projects/" + loctime #服务器执行路径
        experdir="/data/nas/175/蛋白实验室-蛋白实验部-cm/1.完整项目报告/1.正式项目/" + loctime#实验报告路径
        ##1、NAS盘数据处理
        os.chdir(prodir)
        # 下载分析单
        file_name = wget.download(args.link)
        print(file_name) 
        Proid = file_name.split("_")[0]
        MSfile = glob.glob("/data/hstore3/localPro/6.Temp/{fn}*".format(fn=Proid))[0]
        if not os.path.exists(MSfile+"/pic/"):
            RPfile = glob.glob(experdir+"/{fn}*".format(fn=Proid))[0]
            system("cp -r '"+RPfile+"'/* '"+MSfile+"'")
        rawfilname = os.path.basename(MSfile)
        #移动下机数据
        if not os.path.exists(prodir+"/"+Proid):
            os.mkdir(prodir+"/"+Proid)
            temd = shutil.move(MSfile,prodir+"/"+Proid)
        #移动确认单
        nasrawd = prodir+"/"+Proid+"/"+rawfilname
        shutil.move(prodir+"/"+file_name,prodir+"/"+Proid+"/"+rawfilname)
    else:
        nasrawd = args.nasdir
        rundir= args.fdir
        Proid = args.proid
    if nasrawd != "./":
        os.chdir(nasrawd)
        system("rm -rf @eaDir")
    ##2、服务器数据分析
        if not os.path.exists(rundir+"/"+Proid):
            os.mkdir(rundir+"/"+Proid)
        os.chdir(rundir+"/"+Proid)
        system("cp -r '{naspd}' projectdata".format(naspd=nasrawd))
    if not os.path.exists("background/"):
        os.mkdir("background/")
    if args.background != "None":
        system("cp -r {bd}/* background/".format(bd=args.background))
    else:
        system("cp -r /data/hstore3/database/background/{org}/* background/".format(org=args.org))
        system("cp -r /data/hstore3/database/kegg/{org}/* background/".format(org=args.org))
    
    pretargets = ["protein_result","diff_getdata_result",
                  "protein_result_diff_cor","protein_result_diff_motif",
                  "protein_result_diff_function","protein_result_diff_enrich","protein_result_diff_gsea","protein_result_diff_ppi",
                  "protein_result_diff_xlsx","protein_result_diff_map","protein_result_diff"]
    judgetargets = False
    for targets in pretargets:
      if targets in args.targets:
        judgetargets = True
        break
    
    if len(args.fcfilter) > 1 and judgetargets:
      targets = ["diff_result"]
    else:
      targets = args.targets
    
    if not os.path.exists(args.tempfilepath):
      if not os.path.isdir(args.rawdatapath):
        with zipfile.ZipFile(args.rawdatapath, 'r') as zip_ref: 
          for info in zip_ref.infolist(): 
            filename = info.filename.encode('cp437').decode('gbk') 
            zip_ref.extract(info, path='.', pwd=None) 
            os.rename(os.path.join('.', info.filename), os.path.join('.', filename))
          for info in zip_ref.infolist():
            filename = info.filename.encode('cp437').decode('gbk')
            os.rename(os.path.join('.', filename), os.path.join('.', "projectdata"))
            shutil.rmtree(os.path.join('.', info.filename),ignore_errors=True)
            break
        args.rawdatapath="projectdata/"
      # system("profetch -ip '"+args.rawdatapath+"' -sp rawdata/ -f "+args.fetch)
      proteintype=profetch(inputpath=args.rawdatapath,
                           savepath="./rawdata/",fetch=args.fetch,
                           missvarremovepercent=args.missvarremovepercent,
                           missvarremovebygroup=args.missvarremovebygroup,
                           missvarfillmethod=args.missvarfillmethod,
                           rowNorm=args.rowNorm,
                           transNorm=args.transNorm,
                           scaleNorm=args.scaleNorm,
                           ref=args.ref,
                           filter=args.filter,
                           remainnum=args.remainnum,
                           qcFilter=args.qcFilter,
                           rsd=args.rsd,
                           fillpalette = args.fillpalette,
                           p_method = args.p_method,
                           mul_p_method = args.mul_p_method,
                           p_adjust_method = args.p_adjust_method,
                           paired = args.paired)
      for i in ["missvarremovepercent","missvarremovebygroup","missvarfillmethod",
                "rowNorm","transNorm","scaleNorm",
                "ref","filter","remainnum","qcFilter","rsd"]:
        args2[i] = getpredealparams(omic = args.omic,
                                    MStype = proteintype,
                                    name = i,
                                    value = args2[i])
    else:
      f = open("rawdata/protein_type.txt")
      proteintype = f.readline()
      proteintype = proteintype.strip()
      f.close()

    print(proteintype)
    path0 = os.getcwd()
    print("~项目路径：",path0)
    # if os.path.exists("rawdata/project_title.txt") and args.projectpath == "项目报告":
    #   f = open("rawdata/project_title.txt")
    #   project_title = f.readline()
    #   project_title = project_title.strip()
    #   f.close()
    #   args.projectpath = project_title
    
    # print(targets)
    # print(args)
    args.rawdatafile = args.tempfilepath+"/rawdata/rawdatafile.txt"
    args.datafile = args.tempfilepath+"/rawdata/datafile.txt"
    args.infofile = args.tempfilepath+"/rawdata/infofile.txt"
    args.classfile = "classfile.yaml"
    args.classtypefile = "classtype.xlsx"
    args.comparefile = args.tempfilepath+"/rawdata/compare.yaml"
    
    sys.stderr = FilteredPrinter(filtered_print, sys.stderr)
    sys.stdout = FilteredPrinter(filtered_print2, sys.stdout)
    
    if args.transNorm != "NULL" or args.scaleNorm != "NULL":
      args.log=True
    
    if args.vipfilter > 0 :
      args.diff_ana_rds_savepath = args.tempfilepath+"/diff_ana_vip"
    else:
      args.diff_ana_rds_savepath = args.tempfilepath+"/diff_ana"
    
    args.diff_filter_rds_savepath = args.tempfilepath+"/diff_filter"
    args.config_savepath = args.tempfilepath+"/config"  
    args.mulstatistics_rds_savepath = args.tempfilepath+"/mulstatisticsanalyst" 
    args.enrich_xls_savepath = args.tempfilepath+"/enrich" 
    
    config = {"path":{"tempfilepath":args.tempfilepath,
                      "rawdatapath":"rawdata/",
                      "rawdatafile":args.rawdatafile,
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
                                  "keggdb":"/data/hstore3/database/kegg/" if  args.omic=="P" or args.omic=="PF" else "/data/hstore4/database/kegg/kegg/"},
                        "diff_filter":{"log":args.log,
                                       "fcfilter":args.fcfilter,
                                       "fcfiltertype":args.fcfiltertype,
                                       "errorptofcfilter":args.errorptofcfilter,
                                       "vipfilter":args.vipfilter,
                                       "pfilter":args.pfilter,
                                       "adjpfilter":args.adjpfilter},
                        "imagetype":args.imagetype,
                        "MStype":proteintype,
                        "plot":{"family":args.family,
                                "width":args.width,
                                "height":args.height,
                                "volcanonum":args.volcanonum}}}
                                   
    print(config)
    
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
    
    os.system("find -name .snakemake_timestamp -print0  | xargs -0  rm -rf")
    
    if args.report == "T":
      system("Proreport -c "+args.cloud+" -m "+args.mode)
        
    shutil.rmtree(".snakemake",ignore_errors=True)
    os.system("chmod -R 777 ./* >/dev/null 2>&1")
        
    if args.zip :
      if args.report != "T":
        if os.path.exists("rawdata/project_title.txt"):
          f = open("rawdata/project_title.txt")
          project_title = f.readline()
          project_title = project_title.strip()
          f.close()
        else:
          project_title = args.projectpath
        with tarfile.open(project_title+".tar", "w") as tar:
          tar.add(args.projectpath, arcname=os.path.basename(args.projectpath))
      shutil.rmtree(".snakemake",ignore_errors=True)
      shutil.rmtree("rawdata",ignore_errors=True)
      shutil.rmtree("background",ignore_errors=True)
      shutil.rmtree("background_out",ignore_errors=True)
      shutil.rmtree(args.tempfilepath,ignore_errors=True)
      shutil.rmtree(args.rawdatapath,ignore_errors=True)
      shutil.rmtree(args.projectpath,ignore_errors=True)
      os.remove("classfile.yaml")
      os.remove("classtype.xlsx")
