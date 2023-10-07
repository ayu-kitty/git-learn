from snakemake.utils import validate
from itertools import zip_longest

configfile: "config_getsplot.yaml"

if isinstance(config["compare"],list):
	config["compare"]=dict(zip_longest(config["compare"],[{}],fillvalue={}))
	for compare in config["compare"]:
	  config["compare"][compare]={}

validate(config, "config_getsplot.schema.yaml")

rule getsplot_all:
    input:
        expand("{projectpath}/{mode}-splot-{compare}.xls",
               projectpath=config["path"]["projectpath"],
               mode=["PLS-DA","OPLS-DA"],
               compare=config["compare"])

rule getsplot:
    input:
        rdspath=config["path"]["mulstatistics_rds_savepath"]+"/mulstatistic_{mode}-{compare}.rds"
    output:
        "{datasavepath}/{mode,(OPLS-DA)|(PLS-DA)}-splot-{compare}.xls"
    script:
        "getsplot.R"

