import yaml
import os
from snakemake.utils import validate

configfile: "config_mulstatisticana_whole.yaml"
# 验证config
validate(config, "config_mulstatisticana_whole.schema.yaml")

rule mul_whole:
    input:
        expand("{mulstatistics_rds_savepath}/mulstatistic_{mode}-{compare}.rds",
               mulstatistics_rds_savepath=config["path"]["mulstatistics_rds_savepath"],
               mode=["PCA","PLS-DA","OPLS-DA"],
               compare=config["compare"])
               
if "params" not in config :
  config["params"] = {}
  config["params"]["mulboth"] = False
else:
  if "mulboth" not in config["params"]:
    config["params"]["mulboth"] = False

module pca:
    snakefile: "../pca/mulstatistic_pca.smk"
    config: config
    skip_validation: False

use rule * from pca as mul_ana_*

module pls:
    snakefile: "../pls/mulstatistic_pls.smk"
    config: config
    skip_validation: False

use rule * from pls as mul_ana_*

module opls:
    snakefile: "../opls/mulstatistic_opls.smk"
    config: config
    skip_validation: False

use rule * from opls as mul_ana_*

# if not os.path.exists(config["path"]["config_savepath"]): 
# 	os.makedirs(config["path"]["config_savepath"])
# with open("{path}/config_mulstatistic_ana_whole.yaml".format(path=config["path"]["config_savepath"]),"w",encoding="utf-8") as f:
# 	yaml.dump(config,f,allow_unicode=True,sort_keys=False)
