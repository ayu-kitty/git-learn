#!/opt/conda/bin/python
import rpy2.robjects as robjects
from lmbio.basic.printfilter import FilteredPrinter,filtered_print,filtered_print2
from lmbio.basic.runsnakemake import runsnakemake
import lmbio
import pandas as pd
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
    parser.add_argument("-pt", "--ptmfile", help="修饰项目文件路径，如：/data/hstore3/Projects/202306/ZLM2023010890", nargs="?", type=str, default='None')
    parser.add_argument("-pr", "--profile", help="蛋白项目文件路径，如：/data/hstore3/Projects/202306/ZLM2023010891", nargs="?", type=str, default='None')
    parser.add_argument("-b", "--background", help="背景文件夹路径，如：/data/hstore3/Projects/202306/ZLM2023010891/background", nargs="?", type=str, default='None')
    parser.add_argument("-fd", "--fdir", help="运行流程路径，如：/data/hstore3/Personal_dir/pj/ksea-2023", nargs="?", type=str, default='./')
    parser.add_argument("-i", "--proptm", help="项目编号", nargs="?", type=str, default="")
    parser.add_argument("-ip", "--inputpath", help="输入文件路径", nargs="?", type=str, default="./rawdata/")
    parser.add_argument("-r", "--report", help="是否出网页报告，默认T", nargs="?", type=str, default="T")
    parser.add_argument("-m", "--mode", help="网页报告模板类型，默认LM", nargs="?", type=str, default="LM")
    parser.add_argument("-cl", "--cloud", help="是否云交付，默认F", nargs="?", type=str, default="F")
    parser.add_argument("-om","--omic",default = "P", help = "项目类型,P、PF")
    parser.add_argument("-or","--org",default = "hsa", help = "物种简写",required = True)
    parser.add_argument("-rd","--rawdatapath",default = "projectdata/", help = "初始数据保存路径")
        
    args = parser.parse_args()
    
    fdir = args.fdir
    proptm = args.proptm
    if not os.path.exists(fdir+'/'+proptm):
        os.mkdir(fdir+'/'+proptm)
    os.chdir(fdir+'/'+proptm)
    if not os.path.exists('background/'):
        os.system('cp -r {back} ./'.format(back=args.background))
    if not os.path.exists('rawdata/'):
        os.system("proptm_copy_files -pt '"+args.ptmfile+"' -pr '"+args.profile+"' -i '"+args.proptm+"'")
    os.system("proptm_Venn -ip '"+args.inputpath+"' -s ./report/result/2.Related_protein/Venn/")
    os.system("proptm_table -ip '"+args.inputpath+"' -s ./report/result/2.Related_protein/")
    os.system("proptm_Scatterplot -ip '"+args.inputpath+"' -sp ./report/result/3.Different_protein/Scatterplot/")
    os.mkdir('./report/result/1.Project_information/')
    os.system("cp  '{path}/ratio.xlsx' ./report/result/1.Project_information/".format(path=args.inputpath))
    #os.system("cp  '{path}/*表达矩阵*' ./report/result/2.Related_protein/".format(path=args.inputpath))
    #os.system("cp  '{path}/' ./report/result/2.Related_protein/".format(path=args.inputpath))
    os.system("proptm_Heatmap_FC -ip ./ -s ./report/result/3.Different_protein/Heatmap_FC/")
    os.system("proptm_Heatmap_Expression -ip ./ -sp ./report/result/3.Different_protein/Heatmap_Expression/ -rp '"+args.inputpath+"'")
    os.system("cp '去除本底差异数据.xlsx'  ./report/result/3.Different_protein/")
    os.system("proptm_function -o '"+args.org+"'")
    os.system("mv Enrichment/PPI_network/ report/result/5.PPI_network/")
    os.system("mv Enrichment/ report/result/4.Enrichment/")
    os.system("rm '去除本底差异数据.xlsx'")
    org_type=['hsa','mmu','rno']
    if pd.read_table('rawdata/ptm_type.txt',header=None).iloc[1,0]=='磷酸化':
        if args.org in org_type:
            os.system("KSEA_plot -s 6.KSEA/KSEA_Activity/ -o '"+args.org+"'")
            os.system("KSEA_Heatmap -s 6.KSEA/KSEA_Heatmap/ -ip 6.KSEA/KSEA_Activity/")
            os.system("rm 'Rplots.pdf'")
            os.system("KSEA_Substrate -s 6.KSEA/Network/ -o '"+args.org+"'")
            os.system("mv 6.KSEA report/result/6.KSEA")
    os.system("ptmreport -pt '"+args.proptm+"'")






