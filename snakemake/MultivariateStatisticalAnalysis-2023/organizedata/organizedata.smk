import pandas as pd
import yaml
import os
from snakemake.utils import validate
from copy import deepcopy

configfile: "config_organizeddata.yaml"

if "compare" not in config:
  compare=pd.read_excel(config["path"]["rawfile"],"比较",index_col=0).T
  compare.columns=[i.replace("/","-vs-") for i in compare.columns]
  config["compare"]=compare.to_dict()

if isinstance(config["compare"],list):
	config["compare"]={}.fromkeys(config["compare"])

# 参数列表
if "params" not in config:
  config["params"]={}

if "sample" not in config["params"]:
  sample=pd.read_excel(config["path"]["rawfile"],"分组")
  sample=sample.iloc[:,0].values.tolist()
  config["params"]["sample"]=sample

if "all_group" not in config["params"]:
  config["params"]["all_group"]={}
  group=pd.read_excel(config["path"]["rawfile"],"分组")
  listdata=list()
  for (columnName, columnData) in group.iloc[:,1:].iteritems():
    listdata2=columnData.values.tolist()
    listdata2=[i for i in listdata2 if i != '' and not pd.isna(i)]
    if "QC" in listdata2:
      listdata.append("QC")
    for i in listdata2:
      if i not in listdata:
        listdata.append(i)

if "all_compare" not in config["params"]:
  config["params"]["all_compare"]=deepcopy(config["compare"])
  for compare in config["params"]["all_compare"]:
    config["params"]["all_compare"][compare].update({"group":compare.split("-vs-")})
    config["params"]["all_compare"][compare].update({"groupnum":len(config["params"]["all_compare"][compare]["group"])})

if "main_group" not in config["params"]:
  group=pd.read_excel(config["path"]["rawfile"],"分组")
  group=group.iloc[:,1].values.tolist()
  group=[i for i in group if i != '' and not pd.isna(i)]
  listdata=[]
  if "QC" in group:
    listdata.append("QC")
  for i in group:
    if i not in listdata:
      listdata.append(i)
  config["params"]["main_group"]=listdata
  config["compare"].update({"All":{"name":"-vs-".join(config["params"]["main_group"])}})

if "main_group-qc" not in config["params"]:
  config["params"]["main_group-qc"]=[i for i in config["params"]["main_group"] if i != '' and not pd.isna(i) and i != "QC"]
  config["compare"].update({"Allsample":{"name":"-vs-".join(config["params"]["main_group-qc"])}})

for compare in config["compare"]:
	if not isinstance(config["compare"][compare],dict):
		config["compare"][compare]={}
	if "name" not in config["compare"][compare]:
		config["compare"][compare].update({"name":compare})
	if "rawname" not in config["compare"][compare]:
		config["compare"][compare].update({"rawname":config["compare"][compare]["name"].replace("-vs-","/")})
	if "group" not in config["compare"][compare]:
		config["compare"][compare].update({"group":config["compare"][compare]["rawname"].split("/")})
		config["compare"][compare].update({"groupnum":len(config["compare"][compare]["group"])})
	if "paired" not in config["compare"][compare]:
	  config["compare"][compare].update({"paired":False})
	if "p-value_method" not in config["compare"][compare]:
	  if config["compare"][compare]["groupnum"] > 2:
	    config["compare"][compare].update({"p-value_method":"oneway.test"})
	  else:
	    config["compare"][compare].update({"p-value_method":"t.test"})
	if "adjust_p-value_method" not in config["compare"][compare]:
	  config["compare"][compare].update({"adjust_p-value_method":"BH"})

if "main_compare" not in config["params"]:
  config["params"]["main_compare"] = [i for i in config["compare"] if all(elements in config["params"]["main_group"] for elements in config["compare"][i]["group"]) and i != "All" and i != "Allsample"]

validate(config, "config_organizeddata.schema.yaml")

if not os.path.exists(config["path"]["config_savepath"]): 
	os.makedirs(config["path"]["config_savepath"])
with open("{path}/config_organizedata.yaml".format(path=config["path"]["config_savepath"]),"w",encoding="utf-8") as f:
	yaml.dump(config,f,allow_unicode=True,sort_keys=False)

rule organizedata:
    input:
        rawfile=config["path"]["rawfile"]
    output:
        datafile=config["path"]["datafile"],
        infofile=config["path"]["infofile"],
        classfile=config["path"]["classfile"],
        classtypefile=config["path"]["classtypefile"]
    script:
        "organizedata.R"
