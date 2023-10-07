#!/opt/conda/bin/python
from lmbio.basic import path
import pandas as pd
import argparse
import snakemake
import os
import yaml
import tarfile
import shutil
import sys

def cytofana(config,cores = 5,targets = None):
    snakefile = path.packagepath(path = "snakemake/CytofAnalysis-2022/cytofanalysis.smk")
    snakemake.snakemake(snakefile,config = config,configfiles = None,
                        cores = cores,targets = targets,
                        quiet = True,keepgoing = True,verbose = False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mode",default = "LM",help = "公司模板")
    parser.add_argument("-pr","--projectname",default = None, help = "保存目录")
    parser.add_argument("-ad","--analysis_id",default = None, help = "分析编号")
    parser.add_argument("-t","--targets",default = None,nargs="*",help = "运行规则",
                        choices = [None,"metadata_panel"])
    parser.add_argument("-z","--zip",default = False, help = "是否进行压缩",action='store_true')
    parser.add_argument("-c","--cores",default = 5, type= int,help = "线程数")
    
    args = parser.parse_args()
    
    # print(args)
    
    projectname = args.projectname
    
    if projectname is None:
      if not os.path.exists('项目登记单.xlsx'):
        if args.analysis_id is None:
          raise ValueError('未找到项目登记单,请使用-ad输入分析编号')
        else:
          os.system("cytof_1.1_GetAnalystInfo -ad "+args.analysis_id)
      if not os.path.exists('项目登记单.xlsx'):
        raise ValueError('分析编号异常,未查到线上登记单')
      project = pd.read_excel('项目登记单.xlsx')
      projectdata = dict(zip(project['key'], project['value']))
      if projectdata['客户名称'] == projectdata['联系人']:
        projectname = projectdata['项目编号'] + "-" + projectdata['客户名称'] + "-" + projectdata['项目类型'] + "-" + '项目报告'
      else:
        projectname = projectdata['项目编号'] + "-" + projectdata['客户名称'] + "-" + projectdata['联系人'] + "-" + projectdata['项目类型'] + "-" + '项目报告'

    
    config = {"mode":args.mode,
              "projectname":projectname,
              "analysis_id":args.analysis_id}
                                   
    # print(config)
    
    cytofana(config = config,
             cores = args.cores,
             targets = args.targets)
    
    if args.zip :
      with tarfile.open(projectname+".tar", "w") as tar:
        tar.add(projectname, arcname=os.path.basename(projectname))
      os.remove("log.txt")
      shutil.rmtree(".snakemake")
      shutil.rmtree("oecloud")
      shutil.rmtree(args.projectpath)
