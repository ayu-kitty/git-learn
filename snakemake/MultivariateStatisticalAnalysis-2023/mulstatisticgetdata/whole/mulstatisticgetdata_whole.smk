import pandas as pd
import yaml
import os
from snakemake.utils import validate
from itertools import zip_longest

configfile: "config_mulstatisticgetdata_whole.yaml"

if isinstance(config["compare"],list):
	config["compare"]=dict(zip_longest(config["compare"],[{}],fillvalue={}))
	for compare in config["compare"]:
	  config["compare"][compare]={}

# 验证config
validate(config, "config_mulstatisticgetdata_whole.schema.yaml")

rule getdata_whole:
    input:
        expand("{projectpath}/{compare}/{mode}-score-{compare}.xls",
               projectpath=config["path"]["projectpath"],
               mode=config["params"]["mulstatistics_mode"],
               compare={i for i in config["compare"]  if config["compare"][i]["groupnum"]==2}),
        expand("{projectpath}/{compare}/{mode}-score-{compare}.xls",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i != "OPLS-DA"],
               compare={i for i in config["compare"]  if config["compare"][i]["groupnum"]>2}),
        expand("{projectpath}/{compare}/{mode}-loading-{compare}.xls",
               projectpath=config["path"]["projectpath"],
               mode=config["params"]["mulstatistics_mode"],
               compare={i for i in config["compare"]  if config["compare"][i]["groupnum"]==2}),
        expand("{projectpath}/{compare}/{mode}-loading-{compare}.xls",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i != "OPLS-DA"],
               compare={i for i in config["compare"]  if config["compare"][i]["groupnum"]>2}),
        expand("{projectpath}/{compare}/{mode}-splot-{compare}.xls",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i != "PCA"],
               compare={i for i in config["compare"]  if config["compare"][i]["groupnum"]==2}),
        expand("{projectpath}/{compare}/{mode}-splot-{compare}.xls",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i == "PLS-DA"],
               compare={i for i in config["compare"]  if config["compare"][i]["groupnum"]>2}),
        expand("{projectpath}/{compare}/summarydata-{compare}.xls",
               projectpath=config["path"]["projectpath"],
               compare=config["compare"]),
        expand("{projectpath}/summarydata.xls",
               projectpath=config["path"]["projectpath"])

module mul_getscore:
    snakefile: "../getscore/getscore.smk"
    config: config
    skip_validation: False

use rule * from mul_getscore as mul_*

module mul_getloading:
    snakefile: "../getloading/getloading.smk"
    config: config
    skip_validation: False

use rule * from mul_getloading as mul_*

module mul_getsplot:
    snakefile: "../getsplot/getsplot.smk"
    config: config
    skip_validation: False

use rule * from mul_getsplot as mul_*

module mul_getsummary:
    snakefile: "../getsummary/getsummary.smk"
    config: config
    skip_validation: False

use rule * from mul_getsummary as mul_*

# if not os.path.exists(config["path"]["config_savepath"]): 
# 	os.makedirs(config["path"]["config_savepath"])
# with open("{path}/config_mulstatistic_getdata_whole.yaml".format(path=config["path"]["config_savepath"]),"w",encoding="utf-8") as f:
# 	yaml.dump(config,f,allow_unicode=True,sort_keys=False)
