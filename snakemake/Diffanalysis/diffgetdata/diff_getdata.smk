import pandas as pd
import yaml
import os
from snakemake.utils import validate
from itertools import zip_longest

configfile: "config_diff_getdata.yaml"

if "params" not in config:
    config["params"]={}

if "plot_permutation" not in config["params"]:   
    config["params"]["plot_permutation"]={}

if "diff_filter" not in config["params"]:   
    config["params"]["diff_filter"]={}
if "enrich" not in config["params"]:   
    config["params"]["enrich"]={}

if "path" not in config:
    config["path"]={}

if isinstance(config["compare"],list):
	config["compare"]=dict(zip_longest(config["compare"],[{}],fillvalue={}))
	for compare in config["compare"]:
	  config["compare"][compare]={}

for compare in config["compare"]:
	if not isinstance(config["compare"][compare],dict):
		config["compare"][compare]={}
	if "rawname" not in config["compare"][compare]:
		config["compare"][compare].update({"rawname":compare.replace("-vs-","/")})
	if "group" not in config["compare"][compare]:
		config["compare"][compare].update({"group":compare.split("-vs-")})

# 验证config
validate(config, "config_diff_getdata.schema.yaml")

# if not os.path.exists(config["path"]["config_savepath"]): 
#     os.makedirs(config["path"]["config_savepath"])
# with open("{path}/config_diff_get.yaml".format(path=config["path"]["config_savepath"]),"w",encoding="utf-8") as f:
#     yaml.dump(config,f,allow_unicode=True,sort_keys=False)

rule getdata_all:
    input:
        datafile=lambda wildcards: config["path"]["datafile"],
        infofile=lambda wildcards: config["path"]["infofile"]
    output:
        "{savepath}/表达矩阵.xls"
    threads: 1
    shell:
        '''
        /opt/conda/bin/Rscript -e "lmbio::get_alldata_file( \
        datafile = '{input.datafile}', \
        infofile = '{input.infofile}', \
        savepath = '{output}')"
        '''

rule getdata_all_cloud:
    input:
        datafile=config["path"]["datafile"],
        infofile=config["path"]["infofile"],
        classfile=config["path"]["classfile"],
        comparefile=config["path"]["comparefile"]
    output:
        "{savepath}/表达矩阵.xlsx"
    threads: 1
    shell:
        '''
        /opt/conda/bin/Rscript -e "lmbio::get_alldata_file_cloud( \
        datafile = '{input.datafile}', \
        infofile = '{input.infofile}', \
        classfile = '{input.classfile}',\
        comparefile = '{input.comparefile}',\
        savepath = '{output}')"
        '''

if config["params"]["omic"] == "M":
  rule getdata_all_cloud2:
      input:
          datafile=config["path"]["datafile"],
          infofile=config["path"]["infofile"],
          classfile=config["path"]["classfile"],
          comparefile=config["path"]["comparefile"]
      output:
          "{savepath}/data_for_oecloud.xlsx"
      threads: 1
      shell:
          '''
          /opt/conda/bin/Rscript -e "lmbio::get_alldata_file_cloud( \
          datafile = '{input.datafile}', \
          infofile = '{input.infofile}', \
          classfile = '{input.classfile}',\
          comparefile = '{input.comparefile}',\
          savepath = '{output}')"
          '''
else:
    rule getdata_all_cloud2:
      input:
          datafile=config["path"]["rawdatafile"],
          infofile=config["path"]["infofile"],
          classfile=config["path"]["classfile"],
          comparefile=config["path"]["comparefile"]
      output:
          "{savepath}/data_for_oecloud.xlsx"
      threads: 1
      shell:
          '''
          /opt/conda/bin/Rscript -e "lmbio::get_alldata_file_cloud( \
          datafile = '{input.datafile}', \
          infofile = '{input.infofile}', \
          classfile = '{input.classfile}',\
          comparefile = '{input.comparefile}',\
          savepath = '{output}')"
          '''

rule getdata_all_cloud3:
    input:
        datafile=config["path"]["datafile"],
        infofile=config["path"]["infofile"],
        classfile=config["path"]["classfile"],
        comparefile=config["path"]["comparefile"]
    output:
        "{savepath}/数据矩阵.xlsx"
    params:
        rawdatafile=config["path"]["tempfilepath"]+"/rawdata/zerodatafile.txt"
    threads: 1
    shell:
        '''
        /opt/conda/bin/Rscript -e "lmbio::get_alldata_file_cloud2( \
        rawdatafile = '{params.rawdatafile}', \
        datafile = '{input.datafile}', \
        infofile = '{input.infofile}', \
        classfile = '{input.classfile}',\
        comparefile = '{input.comparefile}',\
        savepath = '{output}')"
        '''

rule diff_getdata_all:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{compare}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"],
               compare=config["params"]["all_compare"])
    output:
        "{savepath}/差异表达矩阵(未筛选).xlsx"
    params:
        rdspath = config["path"]["diff_filter_rds_savepath"]
    threads: 1
    shell:
        '''
        flow_diff_getalldata \
        --rdspath '{params.rdspath}' \
        --savefilename '{output}'
        '''

rule diff_getdata_filter:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{compare}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"],
               compare=config["params"]["all_compare"])
    output:
        "{savepath}/差异表达矩阵.xlsx"
    params:
        rdspath = config["path"]["diff_filter_rds_savepath"]
    threads: 1
    shell:
        '''
        flow_diff_getfilterdata \
        --rdspath '{params.rdspath}' \
        --savefilename '{output}'
        '''
        
rule diff_getdata_filter2:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{compare}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"],
               compare=config["params"]["all_compare"])
    output:
        directory("{savepath}/diff")
    threads: 1
    priority: 10
    shell:
        '''
        flow_diff_getfilterdata2 \
        --rdsname {input} \
        --savepath '{output}'
        flow_diff_getfilterdata2 \
        --rdsname {input} \
        --savepath '{output}' \
        --datafrom 'diffdata'
        '''

rule diff_getdata_filter3:
    input:
        config["path"]["diff_xls"]
    output:
        "{savepath}/diff-data.oecloud"
    threads: 1
    params:
        path = config["path"]["diff_xls"]
    priority: 10
    shell:
        '''
        echo "压缩中" > 文件压缩中.log
        tar -cf '{output}' --directory='{params.path}' .
        rm 文件压缩中.log
        '''

rule diff_plot_num:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{compare}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"],
               compare=config["params"]["all_compare"])
    output:
        directory("{savepath}/Foldchange/")
    params:
        rdspath = config["path"]["diff_filter_rds_savepath"],
        family = config["params"]["plot"]["family"],
        filename = lambda wildcards, output, input: "'"+"' '".join(input)+"'"
    threads: 1
    shell:
        '''
        prodiff_fc \
        --inputfile {params.filename} \
        --savepath '{output}' \
        --fontfamily '{params.family}'
        '''

rule diff_getdata_boxplot_top50:
    input:
        rdspath=config["path"]["diff_filter_rds_savepath"]+"/diff_filter-{compare}.rds"
    output:
        "{savepath}/boxplot_top50-{compare}.xls"
    threads: 1
    shell:
        '''
        flow_diff_getheatmapdata \
        --rdspath '{input.rdspath}' \
        --savepath '{output}' \
        --top 50
        '''

# 质谱图
module mzmap:
    snakefile: "mzmap.smk"
    config: config

use rule * from mzmap

# qcmap
module qcmap:
    snakefile: "qcmap.smk"
    config: config

use rule * from qcmap

# 热图
module heatmap:
    snakefile: "diff_heatmap.smk"
    config: config

use rule * from heatmap

# 棒棒糖图
module lolipopmap:
    snakefile: "diff_lolipopmap.smk"
    config: config

use rule * from lolipopmap

# 箱线图
module boxplot:
    snakefile: "diff_boxplot.smk"
    config: config

use rule * from boxplot

# 相关性图
module corrplot:
    snakefile: "diff_corr.smk"
    config: config

use rule * from corrplot

# 火山图
module volcano:
    snakefile: "diff_volcano.smk"
    config: config

use rule * from volcano

# 韦恩图
module venn:
    snakefile: "diff_venn.smk"
    config: config

use rule * from venn

# 热图
module zscore:
    snakefile: "diff_zscoremap.smk"
    config: config

use rule * from zscore

# 饼图
module pie_graph:
    snakefile: "pie_graph.smk"
    config: config

use rule * from pie_graph

# 脂质
module lipid:
    snakefile: "diff_lipid.smk"
    config: config

use rule * from lipid

# 通路富集
module enrich:
    snakefile: "diff_enrich.smk"
    config: config

use rule * from enrich

# gsea
module gsea:
    snakefile: "diff_gsea.smk"
    config: config

use rule * from gsea
