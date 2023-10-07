import pandas as pd
import yaml
import os
from snakemake.utils import validate
from itertools import zip_longest

configfile: "config_permutation.yaml"

if isinstance(config["compare"],list):
	config["compare"]=dict(zip_longest(config["compare"],[{}],fillvalue={}))
	for compare in config["compare"]:
	  config["compare"][compare]={}

mode = "plot_permutation" 
for compare in config["compare"]:
	if not isinstance(config["compare"][compare],dict):
		config["compare"][compare]={}
	if mode not in config["compare"][compare]:
		config["compare"][compare][mode]=dict()
	if "PLS-DA" in config["compare"][compare] and "permI" in config["compare"][compare]["PLS-DA"] and config["compare"][compare]["PLS-DA"]["permI"] > 0:
	  config["compare"][compare][mode].update({"permutation_mode":"PLS-DA"})
	if "OPLS-DA" in config["compare"][compare] and "permI" in config["compare"][compare]["OPLS-DA"] and config["compare"][compare]["OPLS-DA"]["permI"] > 0:
	  config["compare"][compare][mode].update({"permutation_mode":"OPLS-DA"})
	if "params" in config and isinstance(config["params"],dict) and mode in config["params"] and isinstance(config["params"][mode],dict):
		config["compare"][compare][mode].update(config["params"][mode])
	if "params" in config and isinstance(config["params"],dict) and "plot" in config["params"]:
		config["compare"][compare][mode].update(config["params"]["plot"])


# 验证config
validate(config, "config_permutation.schema.yaml")

# if not os.path.exists(config["path"]["config_savepath"]): 
# 	os.makedirs(config["path"]["config_savepath"])
# with open("{path}/config_mulstatistic_plot_permutation.yaml".format(path=config["path"]["config_savepath"]),"w",encoding="utf-8") as f:
# 	yaml.dump(config,f,allow_unicode=True,sort_keys=False)

rule permutation_all:
    input:
        expand("{projectpath}/permutation-{compare}.{imagetype}",
               projectpath=["多元统计分析"],
               compare=config["compare"],
               imagetype=["jpg","pdf"])

rule permutation:
    input:
        lambda wildcards: config["path"]["mulstatistics_rds_savepath"]+"/mulstatistic_"+config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["permutation_mode"]+"-"+wildcards.compare+".rds"
    output:
        "{savepath}/permutation-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["family"],
        width = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["width"],
        height = lambda wildcards: config["compare"][wildcards.compare.replace("~LCMS","").replace("~GCMS","")][mode]["height"]
    shell:
        '''
        map_mulstatistics_permutation \
        --filename '{input}' \
        --savepath '{wildcards.savepath}' \
        --mapname 'permutation-{wildcards.compare}' \
        --imagetype '{wildcards.imagetype}' \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''
