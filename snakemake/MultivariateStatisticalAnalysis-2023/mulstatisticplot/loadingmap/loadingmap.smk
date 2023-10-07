import pandas as pd
import yaml
import os
from snakemake.utils import validate
from itertools import zip_longest

configfile: "config_loadingmap.yaml"

if isinstance(config["compare"],list):
	config["compare"]=dict(zip_longest(config["compare"],[{}],fillvalue={}))
	for compare in config["compare"]:
	  config["compare"][compare]={}

mode = "plot_loadingmap" 
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
validate(config, "config_loadingmap.schema.yaml")

# if not os.path.exists(config["path"]["config_savepath"]): 
# 	os.makedirs(config["path"]["config_savepath"])
# with open("{path}/config_mulstatistic_plot_loadingmap.yaml".format(path=config["path"]["config_savepath"]),"w",encoding="utf-8") as f:
# 	yaml.dump(config,f,allow_unicode=True,sort_keys=False)

rule loadingmap_all:
    input:
        expand("{projectpath}/{mode}-loading-{compare}.{imagetype}",
               projectpath=["多元统计分析"],
               mode=["PCA","PLS-DA","OPLS-DA"],
               compare=config["compare"],
               imagetype=["jpg","pdf"])

rule loadingmap:
    input:
        "{savepath}/{mode}-loading-{compare}.xls"
    output:
        "{savepath,((?!其他图像).)*}/{mode,(OPLS-DA)|(PLS-DA)|(PCA)}-loading-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["family"],
        width = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["width"],
        height = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["height"]
    shell:
        '''
        map_mulstatistics_loadingmap \
        --filename '{input}' \
        --savepath '{wildcards.savepath}' \
        --mapname '{wildcards.mode}-loading-{wildcards.compare}' \
        --imagetype '{wildcards.imagetype}' \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''
