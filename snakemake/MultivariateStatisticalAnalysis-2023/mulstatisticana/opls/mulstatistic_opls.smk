import pandas as pd
import yaml
import os
from snakemake.utils import validate
from itertools import zip_longest

configfile: "config_mulstatistic_opls.yaml"

if isinstance(config["compare"],list):
	config["compare"]=dict(zip_longest(config["compare"],[{}],fillvalue={}))
	for compare in config["compare"]:
	  config["compare"][compare]={}

mode="OPLS-DA"
for compare in config["compare"]:
	if not isinstance(config["compare"][compare],dict):
		config["compare"][compare]={}
	if "rawname" not in config["compare"][compare]:
		config["compare"][compare].update({"rawname":compare.replace("-vs-","/")})
	if "group" not in config["compare"][compare]:
		config["compare"][compare].update({"group":compare.split("-vs-")})
	if mode not in config["compare"][compare]:
		config["compare"][compare][mode]=dict()
	if "params" in config and isinstance(config["params"],dict) and mode in config["params"] and isinstance(config["params"][mode],dict):
		config["compare"][compare][mode].update(config["params"][mode])

# 验证config
validate(config, "config_mulstatistic_opls.schema.yaml")

# if not os.path.exists(config["path"]["config_savepath"]): 
# 	os.makedirs(config["path"]["config_savepath"])
# with open("{path}/config_mulstatistic_ana_opls.yaml".format(path=config["path"]["config_savepath"]),"w",encoding="utf-8") as f:
# 	yaml.dump(config,f,allow_unicode=True,sort_keys=False)

rule opls_all:
    input:
        expand("{mulstatistics_rds_savepath}/mulstatistic_OPLS-DA-{compare}.rds",
               mulstatistics_rds_savepath=config["path"]["mulstatistics_rds_savepath"],
               compare=config["compare"])
               
if config["params"]["mulboth"]:
  rule opls:
      input:
          datafile=lambda wildcards: config["path"]["datafile"],
          classfile=lambda wildcards: config["path"]["classfile"]
      output:
          rdssavepath=config["path"]["mulstatistics_rds_savepath"]+"/mulstatistic_OPLS-DA-{compare2,([A-Za-z0-9._+-]*-vs-[A-Za-z0-9._+-]*)|(All)|(Allsample)}.rds",
          rdssavepathlc=config["path"]["mulstatistics_rds_savepath"]+"/mulstatistic_OPLS-DA-{compare2,([A-Za-z0-9._+-]*-vs-[A-Za-z0-9._+-]*)|(All)|(Allsample)}~LCMS.rds",
          rdssavepathgc=config["path"]["mulstatistics_rds_savepath"]+"/mulstatistic_OPLS-DA-{compare2,([A-Za-z0-9._+-]*-vs-[A-Za-z0-9._+-]*)|(All)|(Allsample)}~GCMS.rds"
      params:
          group=lambda wildcards: config["compare"][wildcards.compare2]["group"],
          name=lambda wildcards: config["compare"][wildcards.compare2]["rawname"],
          amount=lambda wildcards: config["compare"][wildcards.compare2][mode]["amount"],
          log10L=lambda wildcards: config["compare"][wildcards.compare2][mode]["log10L"],
          orthoI=lambda wildcards: config["compare"][wildcards.compare2][mode]["orthoI"],
          permI=lambda wildcards: config["compare"][wildcards.compare2][mode]["permI"],
          predI=lambda wildcards: config["compare"][wildcards.compare2][mode]["predI"],
          scaleC=lambda wildcards: config["compare"][wildcards.compare2][mode]["scaleC"],
          mode=mode,
          both="T"
      threads: 1
      script:
          "../mulstatistic.R"
else:
  rule opls:
      input:
          datafile=lambda wildcards: config["path"]["datafile"],
          classfile=lambda wildcards: config["path"]["classfile"]
      output:
          rdssavepath=config["path"]["mulstatistics_rds_savepath"]+"/mulstatistic_OPLS-DA-{compare,([A-Za-z0-9._+-]*-vs-[A-Za-z0-9._+-]*)|(All)|(Allsample)}.rds"
      params:
          group=lambda wildcards: config["compare"][wildcards.compare]["group"],
          name=lambda wildcards: config["compare"][wildcards.compare]["rawname"],
          amount=lambda wildcards: config["compare"][wildcards.compare][mode]["amount"],
          log10L=lambda wildcards: config["compare"][wildcards.compare][mode]["log10L"],
          orthoI=lambda wildcards: config["compare"][wildcards.compare][mode]["orthoI"],
          permI=lambda wildcards: config["compare"][wildcards.compare][mode]["permI"],
          predI=lambda wildcards: config["compare"][wildcards.compare][mode]["predI"],
          scaleC=lambda wildcards: config["compare"][wildcards.compare][mode]["scaleC"],
          mode=mode,
          both="F"
      threads: 1
      script:
          "../mulstatistic.R"
