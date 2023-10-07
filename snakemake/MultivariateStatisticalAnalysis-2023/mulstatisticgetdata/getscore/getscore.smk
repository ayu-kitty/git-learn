from snakemake.utils import validate
from itertools import zip_longest

configfile: "config_getscore.yaml"

if isinstance(config["compare"],list):
	config["compare"]=dict(zip_longest(config["compare"],[{}],fillvalue={}))
	for compare in config["compare"]:
	  config["compare"][compare]={}

validate(config, "config_getscore.schema.yaml")

rule getscore_all:
    input:
        expand("{projectpath}/{mode}-score-{compare}.xls",
               projectpath=config["path"]["projectpath"],
               mode=["PCA","PLS-DA","OPLS-DA"],
               compare=config["compare"])

rule getscore:
    input:
        rdspath=config["path"]["mulstatistics_rds_savepath"]+"/mulstatistic_{mode}-{compare}.rds"
    output:
        "{datasavepath}/{mode,(OPLS-DA)|(PLS-DA)|(PCA)}-score-{compare}.xls"
    script:
        "getscore.R"
