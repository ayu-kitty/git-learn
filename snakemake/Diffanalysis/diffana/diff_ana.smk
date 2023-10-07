import pandas as pd
import yaml
import os
from snakemake.utils import validate
from itertools import zip_longest

configfile: "config_diff_ana.yaml"

if "params" not in config:
    config["params"]={}

if "plot_permutation" not in config["params"]:   
    config["params"]["plot_permutation"]={}

if "diff_filter" not in config["params"]:   
    config["params"]["diff_filter"]={}

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
validate(config, "config_diff_ana.schema.yaml")

rule diff_result:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{compare}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"],
               compare=config["params"]["all_compare"])

rule diff_ana:
    input:
        datafile=lambda wildcards: config["path"]["datafile"],
        classfile=lambda wildcards: config["path"]["classfile"],
        infofile=lambda wildcards: config["path"]["infofile"],
        comparefile=lambda wildcards: config["path"]["comparefile"]
    output:
        "{savepath}/diff_ana/diff_ana-{compare}.rds"
    params:
        log = config["params"]["diff_filter"]["log"]
    threads: 1
    shell:
        '''
        flow_diff_ana \
        --comparefile '{input.comparefile}' \
        --datafile '{input.datafile}' \
        --infofile '{input.infofile}' \
        --classfile '{input.classfile}' \
        --mulstatisticsrds 'nonefile' \
        --group '{wildcards.compare}' \
        --savediffrds '{output}' \
        --log {params.log}
        '''

rule diff_ana_vip:
    input:
        datafile=lambda wildcards: config["path"]["datafile"],
        classfile=lambda wildcards: config["path"]["classfile"],
        infofile=lambda wildcards: config["path"]["infofile"],
        comparefile=lambda wildcards: config["path"]["comparefile"],
        mulstatisticsrds = expand("{mulstatistics_rds_savepath}/mulstatistic_OPLS-DA-{{compare}}.rds",
                                  mulstatistics_rds_savepath=config["path"]["mulstatistics_rds_savepath"])
    output:
        "{savepath}/diff_ana_vip/diff_ana-{compare}.rds"
    params:
        log = config["params"]["diff_filter"]["log"]
    threads: 1
    shell:
        '''
        flow_diff_ana \
        --comparefile '{input.comparefile}' \
        --datafile '{input.datafile}' \
        --infofile '{input.infofile}' \
        --classfile '{input.classfile}' \
        --mulstatisticsrds '{input.mulstatisticsrds}' \
        --group '{wildcards.compare}' \
        --savediffrds '{output}' \
        --log {params.log}
        '''
        
rule diff_filter:
    input:
        expand("{diff_ana_rds_savepath}/diff_ana-{compare}.rds",
               diff_ana_rds_savepath=config["path"]["diff_ana_rds_savepath"],
               compare=config["params"]["all_compare"])
    output:
        expand("{diff_filter_rds_savepath}/diff_filter-{compare}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"],
               compare=config["params"]["all_compare"])
    params:
        rdspath = config["path"]["diff_ana_rds_savepath"],
        savepath = config["path"]["diff_filter_rds_savepath"],
        fcfilter = config["params"]["diff_filter"]["fcfilter"],
        fcfiltertype = config["params"]["diff_filter"]["fcfiltertype"],
        errorptofcfilter = config["params"]["diff_filter"]["errorptofcfilter"],
        vipfilter = config["params"]["diff_filter"]["vipfilter"],
        pfilter = config["params"]["diff_filter"]["pfilter"],
        adjpfilter = config["params"]["diff_filter"]["adjpfilter"]
    threads: 1
    shell:
        '''
        flow_diff_filter \
        --rdspath '{params.rdspath}' \
        --savepath '{params.savepath}' \
        --fcfilter {params.fcfilter} \
        --fcfiltertype '{params.fcfiltertype}' \
        --errorptofcfilter {params.errorptofcfilter} \
        --vipfilter {params.vipfilter} \
        --pfilter {params.pfilter} \
        --adjpfilter {params.adjpfilter}
        '''

# 多元统计运算
module mul_all:
    snakefile: "../../MultivariateStatisticalAnalysis-2023/whole/multivariate.smk"
    config: config
    skip_validation: False

use rule * from mul_all
