from snakemake.utils import validate
from itertools import zip_longest

configfile: "config_getsummary.yaml"

if isinstance(config["compare"],list):
	config["compare"]=dict(zip_longest(config["compare"],[{}],fillvalue={}))
	for compare in config["compare"]:
	  config["compare"][compare]={}

if "params" not in config:
  config["params"]={}

validate(config, "config_getsummary.schema.yaml")

rule getsummary_all:
    input:
        expand("{projectpath}/summarydata-all.xls",
               projectpath=config["path"]["projectpath"],
               compare=config["compare"]),
        expand("{projectpath}/{compare}/summarydata-{compare}.xls",
               projectpath=config["path"]["projectpath"],
               compare=config["compare"]),

rule getsummary:
    input:
        expand("{mulstatistics_rds_savepath}/mulstatistic_{mode}-{{compare}}.rds",
               mulstatistics_rds_savepath=config["path"]["mulstatistics_rds_savepath"],
               mode=config["params"]["mulstatistics_mode"])
    output:
        "{datasavepath}/summarydata-{compare}.xls"
    script:
        "getsummary.R"

rule getsummary_allgroup:
    input:
        expand("{mulstatistics_rds_savepath}/mulstatistic_{mode}-{compare}.rds",
               mulstatistics_rds_savepath=config["path"]["mulstatistics_rds_savepath"],
               mode=config["params"]["mulstatistics_mode"],
               compare=config["params"]["all_compare"])
    output:
        "{datasavepath}/summarydata.xls"
    script:
        "getsummary.R"
