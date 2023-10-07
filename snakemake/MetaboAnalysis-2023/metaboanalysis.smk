import pandas as pd
import yaml
import os
from snakemake.utils import validate
from itertools import zip_longest
import glob
import re
import copy

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
    else:
        if not isinstance(config["compare"][compare]["group"], list):
            config["compare"][compare]["group"]=[config["compare"][compare]["group"]]
        
# 验证config
validate(config, "config_metabo_analysis.schema.yaml")

if "all_compare" not in config["params"]:
  config["params"]["all_compare"]={}

compare=config["params"]["all_compare"]
compare_all=dict({"All":config["compare"]["All"]})
compare_allsample=dict({"Allsample":config["compare"]["Allsample"]})
if config["params"]["omic"] == "M" or config["params"]["omic"] == "ML":
  if config["params"]["enrich"]["org"] in ["hsa","mmu","rno"]:
    enrichtype = ["KEGG","Reactome"]
  else:
     enrichtype = ["KEGG"]
  ppitype = []
  gseatype = ["KEGG"]
else:
  if config["params"]["enrich"]["org"] in ["hsa","mmu","rno"]:
    enrichtype = ["KEGG","GO","InterPro","Wikipathways","Reactome"]
  else:
    enrichtype = ["KEGG","GO","InterPro"]
  ppitype = ["PPI_network"]
  gseatype = ["KEGG","GO"]

if re.match("^UntargetBoth.*",config["params"]["projecttype"]["class"]):
  multype = ["~LCMS","~GCMS"]
  config["params"]["mulboth"] = True
else:
  multype = [""]
  config["params"]["mulboth"] = False

print(config["params"]["mulboth"])

module diff_analysis:
    snakefile: "../Diffanalysis/whole/diff_analysis.smk"
    config: config
    skip_validation: False

use rule * from diff_analysis

compare=copy.deepcopy(config["params"]["all_compare"])
compare_allsample=dict({"Allsample":config["compare"]["Allsample"]})
if compare_allsample["Allsample"]["groupnum"] > 2:
  compare.update(compare_allsample)

rule metabo_result_mul:
    input:
          # score、loading、3Dscore图
          expand("{projectpath}/{num}.多元统计分析/{mode}/{mode}-{type}-{compare}{multype}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(4+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 mode=["PCA","OPLS-DA"],
                 compare=[i for i in compare if compare[i]["groupnum"]==2 and compare[i]["samplenum"]>2],
                 imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
                 multype=multype,
                 type=["score","3Dscore"]),
          expand("{projectpath}/{num}.多元统计分析/{mode}/{mode}-{type}-{compare}{multype}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(4+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 mode=["PCA","PLS-DA"],
                 compare=[i for i in compare if compare[i]["groupnum"]>2 and compare[i]["samplenum"]>2],
                 imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
                 multype=multype,
                 type=["score","3Dscore"]),
          # loading图
          expand("{projectpath}/{num}.多元统计分析/{mode}/{mode}-loading-{compare}{multype}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(4+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 mode=["OPLS-DA"],
                 compare=[i for i in compare if compare[i]["groupnum"]==2  and compare[i]["samplenum"]>2],
                 imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
                 multype=multype),
          expand("{projectpath}/{num}.多元统计分析/{mode}/{mode}-loading-{compare}{multype}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(4+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 mode=["PLS-DA"],
                 compare=[i for i in compare if compare[i]["groupnum"]>2 and compare[i]["maxsamplenum"]>1],
                 imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
                 multype=multype),
          # splot图
          expand("{projectpath}/{num}.多元统计分析/{mode}/{mode}-splot-{compare}{multype}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(4+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 mode=["OPLS-DA"],
                 compare=[i for i in compare if compare[i]["groupnum"]==2  and compare[i]["samplenum"]>2],
                 imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
                 multype=multype),
          expand("{projectpath}/{num}.多元统计分析/{mode}/{mode}-splot-{compare}{multype}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(4+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 mode=["PLS-DA"],
                 compare=[i for i in compare if compare[i]["groupnum"]>2 and compare[i]["maxsamplenum"]>1],
                 imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
                 multype=multype),
          # 相应检测排序
          expand("{projectpath}/{num}.多元统计分析/permutation/permutation-{compare}{multype}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(4+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 compare=[i for i in compare if config["params"]["plot_permutation"]["permutation_mode"] in config["params"]["mulstatistics_mode"]  and compare[i]["maxsamplenum"]>1],
                 imagetype=config["params"]["imagetype"],
                 multype=multype),
          # 多元统计参数汇总
          expand("{projectpath}/{num}.多元统计分析/summarydata.xls",
                 projectpath=config["path"]["projectpath"] if len([i for i in config["params"]["all_compare"] if compare[i]["samplenum"]>2]) > 0 else [],
                 num=str(4+(1 if config["params"]["projecttype"]["qc"] else 0)))

compare=config["params"]["all_compare"]

if config["params"]["projecttype"]["qc"]:
  rule metabo_result_pre:
      input:
            expand("raw/质控/{mode}/label/{mode}-{type}-{compare}{multype}.{imagetype}",
                    mode=["PCA"],
                    compare=[i for i in compare_all if compare_all[i]["samplenum"]>2],
                    imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
                    type=["score","3Dscore"],
                    multype=multype),
            expand("raw/质控/Samplestreeplot/SamplesTreeplot.{imagetype}",
                   imagetype=config["params"]["imagetype"]),
            expand("raw/质控/Boxplot/Boxplot.{imagetype}",
                   imagetype=config["params"]["imagetype"]),
            expand("raw/质控/QCcorrplot/QCcorrplot.{imagetype}",
                   imagetype=config["params"]["imagetype"]),
            expand("raw/质控/数据矩阵.xlsx")
  
  rule metabo_result_qc:
      input:
            expand("{projectpath}/3.质量控制/{mode}/{mode}-{type}-{compare}{multype}.{imagetype}",
                    projectpath=config["path"]["projectpath"],
                    mode=["PCA"],
                    compare=[i for i in compare_all if compare_all[i]["samplenum"]>2],
                    imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
                    type=["score","3Dscore"],
                    multype=multype),
            expand("{projectpath}/3.质量控制/Samplestreeplot/SamplesTreeplot.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/3.质量控制/Boxplot/Boxplot.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/3.质量控制/QCcorrplot/QCcorrplot.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/4.数据矩阵/数据矩阵.xlsx",
                   projectpath=config["path"]["projectpath"])
else:
  rule metabo_result_pre:
      input:
            expand("raw/质控/{mode}/label/{mode}-{type}-{compare}{multype}.{imagetype}",
                    mode=["PCA"],
                    compare=[i for i in compare_all if compare_all[i]["samplenum"]>2],
                    imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
                    type=["score","3Dscore"],
                    multype=multype),
            expand("raw/质控/Samplestreeplot/SamplesTreeplot.{imagetype}",
                   imagetype=config["params"]["imagetype"]),
            expand("raw/质控/Boxplot/Boxplot.{imagetype}",
                   imagetype=config["params"]["imagetype"]),
            expand("raw/质控/数据矩阵.xlsx")
            
  rule metabo_result_qc:
      input:
            expand("{projectpath}/3.数据矩阵/数据矩阵.xlsx",
                   projectpath=config["path"]["projectpath"])

if config["params"]["projecttype"]["class"] == "UntargetLipid" or config["params"]["projecttype"]["class"] == "QtargetLipid" or "脂质" in config["params"]["projecttype"]["pro_type"]:
  print(config["params"]["projecttype"]["class"])
  rule metabo_result_lipid:
      input:
            expand("{projectpath}/{num}.差异代谢物/脂质绘图/分类气泡图/classfication_bubbleplot-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                   compare=[i for i in compare if compare[i]["groupnum"]==2],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/{num}.差异代谢物/脂质绘图/碳链长度统计/carbonheatmap-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                   compare=[i for i in compare if compare[i]["groupnum"]==2],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/{num}.差异代谢物/脂质绘图/不饱和度统计/unsaturationheatmap-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                   compare=[i for i in compare if compare[i]["groupnum"]==2],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/{num}.差异代谢物/脂质绘图/气泡分布图/{compare}",
                   projectpath=config["path"]["projectpath"],
                   num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                   compare=[i for i in compare if compare[i]["groupnum"]==2]),
else:
  rule metabo_result_lipid:
      input:
          
rule metabo_result_diff_xlsx:
    input:
          expand("{projectpath}/{num}.差异代谢物/差异表达矩阵.xlsx",
                  projectpath=config["path"]["projectpath"],
                  num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0))),
          expand("{projectpath}/{num}.差异代谢物/差异表达矩阵(未筛选).xlsx",
                  projectpath=config["path"]["projectpath"],
                  num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)))

rule metabo_result_diff_map:
    input:
          expand("{projectpath}/{num}.差异代谢物/Heatmap/heatmap_scale-{compare}.xls",
                 projectpath=config["path"]["projectpath"],
                 num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 compare=[i for i in compare if compare[i]["samplenum"]>2]),
          expand("{projectpath}/{num}.差异代谢物/Heatmap/cluster_none/heatmap-{compare}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 compare=compare,
                 imagetype=config["params"]["imagetype"]),
          expand("{projectpath}/{num}.差异代谢物/Heatmap/cluster_samples/heatmap-{compare}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 compare=compare,
                 imagetype=config["params"]["imagetype"]),
          expand("{projectpath}/{num}.差异代谢物/Heatmap_top50/cluster_none/heatmap_top50-{compare}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 compare=compare,
                 imagetype=config["params"]["imagetype"]),
          expand("{projectpath}/{num}.差异代谢物/Heatmap_top50/cluster_samples/heatmap_top50-{compare}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 compare=compare,
                 imagetype=config["params"]["imagetype"]),
          expand("{projectpath}/{num}.差异代谢物/Volcano/volcano-{compare}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 compare=[i for i in compare if compare[i]["groupnum"]==2 and compare[i]["maxsamplenum"]>1],
                 imagetype=config["params"]["imagetype"]),
          expand("{projectpath}/{num}.差异代谢物/Boxplot_top50/{compare}/",
                 projectpath=config["path"]["projectpath"],
                 num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 compare=[i for i in compare if compare[i]["maxsamplenum"]>=3]),
          expand("{projectpath}/{num}.差异代谢物/Lolipopmap/lolipopmap-{compare}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 compare=[i for i in compare if compare[i]["groupnum"]==2 and compare[i]["maxsamplenum"]>1],
                 imagetype=config["params"]["imagetype"]),
          expand("{projectpath}/{num}.差异代谢物/相关性分析/Correlation_top20/Correlation_top20-{compare}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 compare=[i for i in compare if compare[i]["samplenum"]>5],
                 imagetype=config["params"]["imagetype"]),
          # expand("{projectpath}/{num}.差异代谢物/相关性分析/Corrnetwork_top20/Corrnetwork_top20-{compare}.{imagetype}",
          #        projectpath=config["path"]["projectpath"],
          #        num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
          #        compare=[i for i in compare if compare[i]["samplenum"]>5],
          #        imagetype=config["params"]["imagetype"]),
          expand("{projectpath}/{num}.差异代谢物/Venn/{venn}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 venn=[] if len(config["params"]["all_compare"]) == 1 else ["Venn"],
                 imagetype=config["params"]["imagetype"])
          #expand("{projectpath}/{num}.差异代谢物/Venn/{venn}.{imagetype}",
          #       projectpath=config["path"]["projectpath"],
          #       num=str(5+(1 if config["params"]["projecttype"]["qc"] else 0)),
          #       venn=[] if len(config["params"]["all_compare"]) < 6 else ["Upset"],
          #       imagetype=config["params"]["imagetype"])


rule metabo_result_diff_enrich:
    input:
          expand("{projectpath}/{num}.通路富集/{type}/{compare}/",
                 projectpath=config["path"]["projectpath"],
                 num=str(6+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 compare=compare,
                 type=enrichtype)
                 
rule metabo_result_diff_gsea:
    input:
          expand("{projectpath}/{num}.通路富集/GSEA/{compare}/",
                 projectpath=config["path"]["projectpath"],
                 num=str(6+(1 if config["params"]["projecttype"]["qc"] else 0)),
                 compare=[i for i in compare if compare[i]["groupnum"]==2])

if config["params"]["projecttype"]["class"] == "Target":
  rule metabo_result_diff_function:
    input:
        rules.metabo_result_diff_enrich.input
else:
  rule metabo_result_diff_function:
    input:
        rules.metabo_result_diff_enrich.input,
        rules.metabo_result_diff_gsea.input

if config["params"]["projecttype"]["compare"]:
  rule metabo_result_diff:
      input:
            expand("{projectpath}/1.数据统计/Foldchange/",
                   projectpath=config["path"]["projectpath"]),
            rules.metabo_result_mul.input,
            rules.metabo_result_lipid.input,
            rules.metabo_result_diff_xlsx.input,
            rules.metabo_result_diff_map.input,
            rules.metabo_result_diff_function.input
else:
  rule metabo_result_diff:
      input:
if config["params"]["projecttype"]["compare"]:
  rule metabo_result_oecloud:
      input:
            expand("{projectpath}/OECloud_tools/background.oecloud",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/OECloud_tools/data_for_oecloud.xlsx",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/OECloud_tools/diff-data.oecloud",
                   projectpath=config["path"]["projectpath"])
else:
  rule metabo_result_oecloud:
      input:
            expand("{projectpath}/OECloud_tools/background.oecloud",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/OECloud_tools/data_for_oecloud.xlsx",
                   projectpath=config["path"]["projectpath"])

if config["params"]["projecttype"]["class"] == "Target":
  rule target_qc_file:
    input:
          tic_file = "Tic.docx",
          data_file = "data.xlsx",
          data_note = "data说明文件.docx",
          confirm_file = "分析确认单.xlsx",
          info_file = expand("{pro_type}-标品信息表.xlsx",
          pro_type = config["params"]["projecttype"]["pro_type"])
    output:
          qc_dir = directory(expand("{projectpath}/2.谱图分析/",
                        projectpath=config["path"]["projectpath"])),
          confirm_dir = directory(expand("{projectpath}/7.项目登记单/",
                        projectpath=config["path"]["projectpath"])),
          tic_file = expand("{projectpath}/2.谱图分析/Tic.docx",
                        projectpath=config["path"]["projectpath"]),
          data_file = expand("{projectpath}/OECloud_tools/data.xlsx",
                        projectpath=config["path"]["projectpath"]),
          data_note = expand("{projectpath}/OECloud_tools/data说明文件.docx",
                        projectpath=config["path"]["projectpath"]),
          info_file = expand("{projectpath}/OECloud_tools/{pro_type}-标品信息表.xlsx",
                        projectpath=config["path"]["projectpath"],
                        pro_type = config["params"]["projecttype"]["pro_type"]),
          qc_plot = expand("{projectpath}/2.谱图分析/CV_plot.{imagetype}",
                        projectpath=config["path"]["projectpath"],
                        imagetype=config["params"]["imagetype"])
    threads: 1
    shell:
        '''          
        path=$(pwd)
        if [ ! -d $path/'{output.qc_dir}' ]; then
          mkdir -p $path/'{output.qc_dir}'
        fi
        cp $path/'{input.tic_file}' '{output.tic_file}'
        cp $path/Calibration\ Curve_*.docx '{output.qc_dir}'
        cp $path/'{input.data_file}' '{output.data_file}'
        cp $path/'{input.data_note}' '{output.data_note}'
        cp $path/'{input.info_file}' '{output.info_file}'

        if [ ! -d $path/'{output.confirm_dir}' ]; then
          mkdir -p $path/'{output.confirm_dir}'
        fi
        cp $path/'{input.confirm_file}' '{output.confirm_dir}'
        
        CV_plot \
        --inputfile '{input.data_file}' \
        --outpath '{output.qc_dir}'
        '''
        
  rule metabo_result:
    input:
          expand("{projectpath}/1.数据统计/Pie_Graph/",
                 projectpath=config["path"]["projectpath"]),
          expand("{projectpath}/2.谱图分析/",
                 projectpath=config["path"]["projectpath"]),
          rules.metabo_result_qc.input,
          rules.metabo_result_diff.input,
          rules.metabo_result_oecloud.input,
          expand("{projectpath}/OECloud_tools/data.xlsx",
                 projectpath=config["path"]["projectpath"]),
          expand("{projectpath}/OECloud_tools/data说明文件.docx",
                 projectpath=config["path"]["projectpath"]),
          expand("{projectpath}/7.项目登记单/",
                 projectpath=config["path"]["projectpath"]),
          expand("{projectpath}/OECloud_tools/{pro_type}-标品信息表.xlsx",
                 projectpath=config["path"]["projectpath"],
                 pro_type = config["params"]["projecttype"]["pro_type"])
else:
  rule metabo_result:
    input:
          expand("{projectpath}/1.数据统计/Pie_Graph/",
                 projectpath=config["path"]["projectpath"]),
          expand("{projectpath}/2.色谱图",
                 projectpath=config["path"]["projectpath"]),
          rules.metabo_result_qc.input,
          rules.metabo_result_diff.input,
          rules.metabo_result_oecloud.input

if not os.path.exists(config["path"]["config_savepath"]): 
    os.makedirs(config["path"]["config_savepath"])
with open("{path}/config_metabo_analysis.yaml".format(path=config["path"]["config_savepath"]),"w",encoding="utf-8") as f:
    yaml.dump(config,f,allow_unicode=True,sort_keys=False)
