import pandas as pd
import yaml
import os
from snakemake.utils import validate
from itertools import zip_longest

configfile: "config_3Dscoremap.yaml"

if isinstance(config["compare"],list):
	config["compare"]=dict(zip_longest(config["compare"],[{}],fillvalue={}))
	for compare in config["compare"]:
	  config["compare"][compare]={}

mode = "plot_3Dscoremap" 
for compare in config["compare"]:
	if not isinstance(config["compare"][compare],dict):
		config["compare"][compare]={}
	if mode not in config["compare"][compare]:
		config["compare"][compare][mode]=dict()
	if "params" in config and isinstance(config["params"],dict) and mode in config["params"] and isinstance(config["params"][mode],dict):
		config["compare"][compare][mode].update(config["params"][mode])
	if "params" in config and isinstance(config["params"],dict) and "plot" in config["params"]:
		config["compare"][compare][mode].update(config["params"]["plot"])

# 验证config
validate(config, "config_3Dscoremap.schema.yaml")

# if not os.path.exists(config["path"]["config_savepath"]): 
# 	os.makedirs(config["path"]["config_savepath"])
# with open("{path}/config_mulstatistic_plot_3Dscoremap.yaml".format(path=config["path"]["config_savepath"]),"w",encoding="utf-8") as f:
# 	yaml.dump(config,f,allow_unicode=True,sort_keys=False)

rule scoremap3D_all:
    input:
        expand("{projectpath}/{mode}-3Dscore-{compare}.{imagetype}",
               projectpath=["多元统计分析"],
               mode=["PCA","PLS-DA","OPLS-DA"],
               compare=config["compare"],
               imagetype=["jpg","pdf"])

rule scoremap3D:
    input:
        "{savepath}/{mode}-score-{compare}.xls"
    output:
        "{savepath,((?!其他图像)(?!label).)*}/{mode,(OPLS-DA)|(PLS-DA)|(PCA)}-3Dscore-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["family"],
        width = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["width"],
        height = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["height"],
        classfile = config["path"]["classtypefile"]
    shell:
        '''
        map_mulstatistics_3dscoremap \
        --filename '{input}' \
        --savepath '{wildcards.savepath}' \
        --mapname '{wildcards.mode}-3Dscore-{wildcards.compare}' \
        --mode '{wildcards.mode}' \
        --imagetype '{wildcards.imagetype}' \
        --classfile {params.classfile} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''

rule scoremap3D1_2:
    input:
        "{savepath}/{mode}-score-{compare}.xls"
    output:
        "{savepath}/none_label/{mode,(OPLS-DA)|(PLS-DA)|(PCA)}-3Dscore-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["family"],
        width = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["width"],
        height = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["height"],
        classfile = config["path"]["classtypefile"]
    shell:
        '''
        map_mulstatistics_3dscoremap \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/none_label/' \
        --mapname '{wildcards.mode}-3Dscore-{wildcards.compare}' \
        --mode '{wildcards.mode}' \
        --imagetype '{wildcards.imagetype}' \
        --classfile {params.classfile} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''
        
rule scoremap3D1_3:
    input:
        "{savepath}/{mode}-score-{compare}.xls"
    output:
        "{savepath}/label/{mode,(OPLS-DA)|(PLS-DA)|(PCA)}-3Dscore-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["family"],
        width = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["width"],
        height = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["height"],
        classfile = config["path"]["classtypefile"]
    shell:
        '''
        map_mulstatistics_3dscoremap \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/label/' \
        --mapname '{wildcards.mode}-3Dscore-{wildcards.compare}' \
        --mode '{wildcards.mode}' \
        --imagetype '{wildcards.imagetype}' \
        --classfile {params.classfile} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''