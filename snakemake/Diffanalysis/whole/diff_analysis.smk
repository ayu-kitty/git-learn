import pandas as pd
import yaml
import os
from snakemake.utils import validate
from itertools import zip_longest

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
validate(config, "config_diff_analysis.schema.yaml")

module diff_ana:
    snakefile: "../diffana/diff_ana.smk"
    config: config
    skip_validation: False

use rule * from diff_ana

module diff_getdata:
    snakefile: "../diffgetdata/diff_getdata.smk"
    config: config
    skip_validation: False

use rule * from diff_getdata

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

if config["params"]["omic"] == "M" or config["params"]["omic"] == "ML":
  compare=config["params"]["all_compare"]
  # compare=copy.deepcopy(config["params"]["all_compare"])
  # compare_allsample=dict({"Allsample":config["compare"]["Allsample"]})
  # if compare_allsample["Allsample"]["groupnum"] > 2:
  #   compare.update(compare_allsample)
  rule result_mul:
      input:
            # score、loading、3Dscore图
            expand("{projectpath}/2.多元统计分析/{mode}/{mode}-{type}-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   mode=["PCA","OPLS-DA"],
                   compare=[i for i in compare if compare[i]["groupnum"]==2 and compare[i]["samplenum"]>2],
                   imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
                   type=["score","3Dscore"]),
            expand("{projectpath}/2.多元统计分析/{mode}/{mode}-{type}-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   mode=["PCA","PLS-DA"],
                   compare=[i for i in compare if compare[i]["groupnum"]>2 and compare[i]["samplenum"]>2],
                   imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
                   type=["score","3Dscore"]),
            # loading图
            expand("{projectpath}/2.多元统计分析/{mode}/{mode}-loading-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   mode=["PCA","OPLS-DA"],
                   compare=[i for i in compare if compare[i]["groupnum"]==2  and compare[i]["samplenum"]>2],
                   imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"]),
            expand("{projectpath}/2.多元统计分析/{mode}/{mode}-loading-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   mode=["PCA","PLS-DA"],
                   compare=[i for i in compare if compare[i]["groupnum"]>2 and compare[i]["maxsamplenum"]>1],
                   imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"]),
            # splot图
            expand("{projectpath}/2.多元统计分析/{mode}/{mode}-splot-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   mode=["OPLS-DA"],
                   compare=[i for i in compare if compare[i]["groupnum"]==2  and compare[i]["samplenum"]>2],
                   imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"]),
            expand("{projectpath}/2.多元统计分析/{mode}/{mode}-splot-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   mode=["PLS-DA"],
                   compare=[i for i in compare if compare[i]["groupnum"]>2 and compare[i]["maxsamplenum"]>1],
                   imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"]),
            # 相应检测排序
            expand("{projectpath}/2.多元统计分析/permutation/permutation-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   compare=[i for i in compare if config["params"]["plot_permutation"]["permutation_mode"] in config["params"]["mulstatistics_mode"]  and compare[i]["maxsamplenum"]>1],
                   imagetype=config["params"]["imagetype"]),
            # 多元统计参数汇总
            expand("{projectpath}/2.多元统计分析/summarydata.xls",
                   projectpath=config["path"]["projectpath"] if len([i for i in config["params"]["all_compare"] if compare[i]["samplenum"]>2]) > 0 else [])
else:
  rule result_mul:
      input:
          expand("{projectpath}/result/1.Reliable_results/PCA/none_label/PCA-{type}-{compare}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 compare=[i for i in compare if compare[i]["samplenum"]>2],
                 imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
                 type=["score","3Dscore"]),
          expand("{projectpath}/result/1.Reliable_results/PCA/label/PCA-{type}-{compare}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 compare=[i for i in compare if compare[i]["samplenum"]>2],
                 imagetype=config["params"]["imagetype"],
                 type=["score"]),
          expand("{projectpath}/result/1.Reliable_results/PCA/none_label/PCA-{type}-{compare}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 compare=[i for i in compare_all if compare_all[i]["samplenum"]>2],
                 imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
                 type=["score","3Dscore"]),
          expand("{projectpath}/result/1.Reliable_results/PCA/label/PCA-{type}-{compare}.{imagetype}",
                 projectpath=config["path"]["projectpath"],
                 compare=[i for i in compare_all if compare_all[i]["samplenum"]>2],
                 imagetype=config["params"]["imagetype"],
                 type=["score"])

compare=config["params"]["all_compare"]

if config["params"]["omic"] == "ML":
  rule result_lipid:
      input:
            expand("{projectpath}/{num}.差异代谢物/脂质绘图/分类气泡图/classfication_bubbleplot-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   num=["4"],
                   compare=[i for i in compare if compare[i]["groupnum"]==2],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/{num}.差异代谢物/脂质绘图/碳链长度统计/carbonheatmap-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   num=["4"],
                   compare=[i for i in compare if compare[i]["groupnum"]==2],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/{num}.差异代谢物/脂质绘图/不饱和度统计/unsaturationheatmap-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   num=["4"],
                   compare=[i for i in compare if compare[i]["groupnum"]==2],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/{num}.差异代谢物/脂质绘图/气泡分布图/{compare}",
                   projectpath=config["path"]["projectpath"],
                   num=["4"],
                   compare=[i for i in compare if compare[i]["groupnum"]==2]),
else:
  rule result_lipid:
      input:

if config["params"]["omic"] == "M" or config["params"]["omic"] == "ML":
  rule diff_getdata_result:
      input:
            expand("{projectpath}/1.数据矩阵/数据矩阵.xlsx",
                   projectpath=config["path"]["projectpath"]),
            rules.result_mul.input,
            rules.result_lipid.input,
            expand("{projectpath}/3.差异代谢物/Foldchange/",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/3.差异代谢物/差异表达矩阵.xlsx",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/3.差异代谢物/差异表达矩阵(未筛选).xlsx",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/3.差异代谢物/Heatmap/heatmap_scale-{compare}.xls",
                   projectpath=config["path"]["projectpath"],
                   compare=[i for i in compare if compare[i]["samplenum"]>2]),
            expand("{projectpath}/3.差异代谢物/Heatmap/cluster_none/heatmap-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   compare=compare,
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/3.差异代谢物/Heatmap/cluster_samples/heatmap-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   compare=compare,
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/3.差异代谢物/Heatmap_top50/cluster_none/heatmap_top50-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   compare=compare,
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/3.差异代谢物/Heatmap_top50/cluster_samples/heatmap_top50-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   compare=compare,
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/3.差异代谢物/Volcano/volcano-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   compare=[i for i in compare if compare[i]["groupnum"]==2 and compare[i]["maxsamplenum"]>1],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/3.差异代谢物/Boxplot_top50/{compare}/",
                   projectpath=config["path"]["projectpath"],
                   compare=[i for i in compare if compare[i]["maxsamplenum"]>=3]),
            expand("{projectpath}/3.差异代谢物/Lolipopmap/lolipopmap-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   compare=[i for i in compare if compare[i]["groupnum"]==2 and compare[i]["maxsamplenum"]>1],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/3.差异代谢物/相关性分析/Correlation_top20/Correlation_top20-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   compare=[i for i in compare if compare[i]["samplenum"]>5],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/3.差异代谢物/Venn/{venn}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   venn=[] if len(config["params"]["all_compare"]) == 1 else ["Venn"],
                   imagetype=config["params"]["imagetype"]),
            #expand("{projectpath}/3.差异代谢物/Venn/{venn}.{imagetype}",
            #       projectpath=config["path"]["projectpath"],
            #       venn=[] if len(config["params"]["all_compare"]) < 6 else ["Upset"],
            #       imagetype=config["params"]["imagetype"]),
            # 通路富集
            expand("{projectpath}/4.通路富集/{type}/{compare}/",
                   projectpath=config["path"]["projectpath"],
                   compare=compare,
                   type=enrichtype),
            expand("{projectpath}/4.通路富集/GSEA/{compare}/",
                   projectpath=config["path"]["projectpath"],
                   compare=[i for i in compare if compare[i]["groupnum"]==2]),
            expand("{projectpath}/OECloud_tools/background.oecloud",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/OECloud_tools/data_for_oecloud.xlsx",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/OECloud_tools/diff-data.oecloud",
                   projectpath=config["path"]["projectpath"])
else:
  rule diff_getdata_result:
      input:
            expand("{projectpath}/result/1.Reliable_results/Samplestreeplot/SamplesTreeplot.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/result/1.Reliable_results/Samplecorrplot/Samplecorrplot.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/result/1.Reliable_results/Boxplot/Boxplot.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/result/1.Reliable_results/Densityplot/Densityplot.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/result/1.Reliable_results/RSD/RSD-boxplot.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   imagetype=config["params"]["imagetype"]),
            rules.result_mul.input,
            expand("{projectpath}/result/2.Different_expression/Foldchange",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/result/2.Different_expression/表达矩阵.xlsx",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/result/2.Different_expression/差异表达矩阵.xlsx",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/result/2.Different_expression/差异表达矩阵(未筛选).xlsx",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/result/2.Different_expression/Heatmap/heatmap_scale-{compare}.xls",
                   projectpath=config["path"]["projectpath"],
                   compare=[i for i in compare if compare[i]["samplenum"]>2]),
            expand("{projectpath}/result/2.Different_expression/Heatmap/cluster_none/heatmap-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   compare=compare,
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/result/2.Different_expression/Heatmap/cluster_samples/heatmap-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   compare=compare,
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/result/2.Different_expression/Volcano/volcano-{compare}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   compare=[i for i in compare if compare[i]["groupnum"]==2 and compare[i]["minsamplenum"]>1],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/result/2.Different_expression/Venn/{venn}.{imagetype}",
                   projectpath=config["path"]["projectpath"],
                   venn=[] if len(config["params"]["all_compare"]) == 1 else ["Venn"],
                   imagetype=config["params"]["imagetype"]),
            expand("{projectpath}/result/3.Enrichment/{type}/{compare}/",
                   projectpath=config["path"]["projectpath"],
                   compare=compare,
                   type=enrichtype),
            expand("{projectpath}/result/4.{type}/{compare}/",
                   projectpath=config["path"]["projectpath"],
                   compare=compare,
                   type=ppitype),
            expand("{projectpath}/result/5.GSEA/{type}/{compare}/",
                   projectpath=config["path"]["projectpath"],
                   compare=[i for i in compare if compare[i]["groupnum"]==2],
                   type=gseatype),
            expand("{projectpath}/result/OECloud_tools/background.oecloud",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/result/OECloud_tools/data_for_oecloud.xlsx",
                   projectpath=config["path"]["projectpath"]),
            expand("{projectpath}/result/OECloud_tools/diff-data.oecloud",
                   projectpath=config["path"]["projectpath"])

if not os.path.exists(config["path"]["config_savepath"]): 
    os.makedirs(config["path"]["config_savepath"])
with open("{path}/config_diff.yaml".format(path=config["path"]["config_savepath"]),"w",encoding="utf-8") as f:
    yaml.dump(config,f,allow_unicode=True,sort_keys=False)
