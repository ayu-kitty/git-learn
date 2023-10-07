import pandas as pd
import yaml
import os
from snakemake.utils import validate
from itertools import zip_longest

configfile: "config_splotmap2.yaml"

if isinstance(config["compare"],list):
	config["compare"]=dict(zip_longest(config["compare"],[{}],fillvalue={}))
	for compare in config["compare"]:
	  config["compare"][compare]={}

mode = "plot_splotmap2" 
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
validate(config, "config_splotmap2.schema.yaml")

# if not os.path.exists(config["path"]["config_savepath"]): 
# 	os.makedirs(config["path"]["config_savepath"])
# with open("{path}/config_mulstatistic_plot_splotmap2.yaml".format(path=config["path"]["config_savepath"]),"w",encoding="utf-8") as f:
# 	yaml.dump(config,f,allow_unicode=True,sort_keys=False)

rule splotmap2_all:
    input:
        expand("{projectpath}/{mode}-loading2-{compare}.{imagetype}",
               projectpath=["多元统计分析"],
               mode=["PCA","PLS-DA","OPLS-DA"],
               compare=config["compare"],
               imagetype=["jpg","pdf"])

rule splotmap2:
    input:
        "{savepath}/{mode}-splot-{compare}.xls"
    output:
        "{savepath}/{mode,(OPLS-DA)|(PLS-DA)|(PCA)}-splot2-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["family"],
        width = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["width"],
        height = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["height"]
    shell:
        '''
        map_mulstatistics_splotmap2 \
        --filename '{input}' \
        --savepath '{wildcards.savepath}' \
        --mapname '{wildcards.mode}-splot2-{wildcards.compare}' \
        --imagetype '{wildcards.imagetype}' \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''
        
rule splotmap2_2:
    input:
        "{savepath}/{mode}-splot-{compare}.xls"
    output:
        "{savepath}/其他图像/{mode,(OPLS-DA)|(PLS-DA)|(PCA)}-splot-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["family"],
        width = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["width"],
        height = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["height"]
    shell:
        '''
        map_mulstatistics_splotmap2 \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/其他图像/' \
        --mapname '{wildcards.mode}-splot-{wildcards.compare}' \
        --imagetype '{wildcards.imagetype}' \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''

