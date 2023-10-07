import pandas as pd
import yaml
import os
from snakemake.utils import validate
from itertools import zip_longest

# configfile: "config_multivariate.yaml"

if "params" not in config:
    config["params"]={}

validate(config, "config_multivariate.schema.yaml")

# 数据整理
# module organizedata:
#     snakefile: "../organizedata/organizedata.smk"
#     config: config
#     skip_validation: False

# use rule * from organizedata

# 多元统计运算
module mul_ana:
    snakefile: "../mulstatisticana/whole/mulstatisticana_whole.smk"
    config: config
    skip_validation: False

use rule * from mul_ana

# 多元统计数据提取
module mul_getdata:
    snakefile: "../mulstatisticgetdata/whole/mulstatisticgetdata_whole.smk"
    config: config
    skip_validation: False

use rule * from mul_getdata

# 多元统计绘图
module mul_plot:
    snakefile: "../mulstatisticplot/whole/mulstatisticplot_whole.smk"
    config: config
    skip_validation: False

use rule * from mul_plot

compare=config["params"]["all_compare"]
compare_all=dict({"All":config["compare"]["All"]})
compare_allsample=dict({"Allsample":config["compare"]["Allsample"]})

rule predeal_data:
    input:
        rawdatafile = config["path"]["rawdatafile"],
        classfile = config["path"]["classfile"]
    output:
        resultfile = config["path"]["datafile"]
    params:
        processsavepath = config["params"]["predeal"]["processsavepath"],
        missvarremovepercent = config["params"]["predeal"]["missvarremovepercent"],
        missvarremovebygroup = config["params"]["predeal"]["missvarremovebygroup"],
        missvarfillmethod = config["params"]["predeal"]["missvarfillmethod"],
        rowNorm = config["params"]["predeal"]["rowNorm"],
        transNorm = config["params"]["predeal"]["transNorm"],
        scaleNorm = config["params"]["predeal"]["scaleNorm"],
        ref = config["params"]["predeal"]["ref"],
        filter = config["params"]["predeal"]["filter"],
        remainnum = int(config["params"]["predeal"]["remainnum"]),
        qcFilter = config["params"]["predeal"]["qcFilter"],
        rsd = int(config["params"]["predeal"]["rsd"]),
    shell:
        '''
        flow_predealdata_Analyze \
        -df '{input.rawdatafile}' \
        -cf '{input.classfile}' \
        -s '{output.resultfile}' \
        -pp {params.processsavepath} \
        -mrp '{params.missvarremovepercent}' \
        -mrg '{params.missvarremovebygroup}' \
        -mvf '{params.missvarfillmethod}' \
        -rn '{params.rowNorm}' \
        -tn '{params.transNorm}' \
        -sn '{params.scaleNorm}' \
        -r '{params.ref}' \
        -fi '{params.filter}' \
        -re '{params.remainnum}' \
        -qf '{params.qcFilter}' \
        -rs '{params.rsd}'
        '''

rule mul_result_all:
    input:
        # score、loading、3Dscore图
        expand("{projectpath}/{compare}/{mode}-{type}-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               mode=config["params"]["mulstatistics_mode"],
               compare=[i for i in compare if compare[i]["groupnum"]==2 and compare[i]["samplenum"]>2],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
               type=["score","3Dscore"]),
        expand("{projectpath}/{compare}/{mode}-{type}-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i != "OPLS-DA" and i != "PLS-DA"],
               compare=[i for i in compare if compare[i]["groupnum"]>2 and compare[i]["samplenum"]>2],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
               type=["score","3Dscore"]),
        expand("{projectpath}/{compare}/{mode}-{type}-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i != "OPLS-DA" and i != "PCA"],
               compare=[i for i in compare if compare[i]["groupnum"]>2 and compare[i]["maxsamplenum"]>1],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
               type=["score","loading","3Dscore"]),
        # loading图
        expand("{projectpath}/{compare}/{mode}-loading-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i != "PCA"],
               compare=[i for i in compare if compare[i]["groupnum"]==2  and compare[i]["samplenum"]>2],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"]),
        expand("{projectpath}/{compare}/{mode}-loading-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i == "PLS-DA"],
               compare=[i for i in compare if compare[i]["groupnum"]>2 and compare[i]["maxsamplenum"]>1],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"]),
        # splot图
        expand("{projectpath}/{compare}/{mode}-splot-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i != "PCA"],
               compare=[i for i in compare if compare[i]["groupnum"]==2  and compare[i]["samplenum"]>2],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"]),
        expand("{projectpath}/{compare}/{mode}-splot-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i == "PLS-DA"],
               compare=[i for i in compare if compare[i]["groupnum"]>2 and compare[i]["maxsamplenum"]>1],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"]),
        # 相应检测排序
        expand("{projectpath}/{compare}/permutation-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               compare=[i for i in compare if config["params"]["plot_permutation"]["permutation_mode"] in config["params"]["mulstatistics_mode"] and compare[i]["maxsamplenum"]>1],
               imagetype=config["params"]["imagetype"]),
        # 多元统计参数汇总
        expand("{projectpath}/{compare}/summarydata-{compare}.xls",
                 projectpath=config["path"]["projectpath"],
                 compare=[i for i in compare if ("PCA" in config["params"]["mulstatistics_mode"] and compare[i]["samplenum"]>2) or ("PLS-DA" in config["params"]["mulstatistics_mode"] and compare[i]["maxsamplenum"]>1) or ("OPLS-DA" in config["params"]["mulstatistics_mode"] and compare[i]["maxsamplenum"]>1)]),
        expand("{projectpath}/summarydata.xls",
               projectpath=config["path"]["projectpath"])
               
rule mul_result_score:
    input:
        # score、3Dscore图
        expand("{projectpath}/{compare}/{mode}-{type}-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               mode=config["params"]["mulstatistics_mode"],
               compare=[i for i in compare_all if compare_all[i]["groupnum"]==2 and compare_all[i]["samplenum"]>2],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
               type=["score","3Dscore"]),
        expand("{projectpath}/{compare}/{mode}-{type}-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i != "OPLS-DA" and i != "PLS-DA"],
               compare=[i for i in compare_all if compare_all[i]["groupnum"]>2],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
               type=["score","3Dscore"]),
        expand("{projectpath}/{compare}/{mode}-{type}-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i != "OPLS-DA" and i != "PCA"],
               compare=[i for i in compare_all if compare_all[i]["groupnum"]>2 and compare_all[i]["maxsamplenum"]>1],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
               type=["score","3Dscore"])
  
rule mul_result_score2:
    input:
        # score、3Dscore图
        expand("{projectpath}/{compare}/{mode}-{type}-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               mode=config["params"]["mulstatistics_mode"],
               compare=[i for i in compare_allsample if compare_allsample[i]["groupnum"]==2 and compare_allsample[i]["samplenum"]>2],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
               type=["score","3Dscore"]),
        expand("{projectpath}/{compare}/{mode}-{type}-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i != "OPLS-DA"],
               compare=[i for i in compare_allsample if compare_allsample[i]["groupnum"]>2],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
               type=["score","3Dscore"]),
        expand("{projectpath}/{compare}/{mode}-{type}-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               mode=[i for i in config["params"]["mulstatistics_mode"] if i != "OPLS-DA"],
               compare=[i for i in compare_allsample if compare_allsample[i]["groupnum"]>2 and compare_allsample[i]["maxsamplenum"]>1],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
               type=["score","3Dscore"]),
        # 相应检测排序
        expand("{projectpath}/{compare}/permutation-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               compare=[i for i in compare_allsample if config["params"]["plot_permutation"]["permutation_mode"] in config["params"]["mulstatistics_mode"] and compare_allsample[i]["maxsamplenum"]>1],
               imagetype=config["params"]["imagetype"])

if not os.path.exists(config["path"]["config_savepath"]): 
    os.makedirs(config["path"]["config_savepath"])
with open("{path}/config_multivariate.yaml".format(path=config["path"]["config_savepath"]),"w",encoding="utf-8") as f:
    yaml.dump(config,f,allow_unicode=True,sort_keys=False)