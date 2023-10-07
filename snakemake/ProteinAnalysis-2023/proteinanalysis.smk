import pandas as pd
import yaml
import os
from snakemake.utils import validate
from itertools import zip_longest
import glob
import re

for file in os.listdir(config["path"]["rawdatapath"]):
  if re.match("(DLM|ZLM|LM|DQD|DOE|DZLM|DZQD|确认单).*.xlsx",file):
    registrationfile = file
    confirm_file = config["path"]["rawdatapath"]+file

print(confirm_file)

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
validate(config, "config_protein_analysis.schema.yaml")

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


module diff_analysis:
    snakefile: "../Diffanalysis/whole/diff_analysis.smk"
    config: config
    skip_validation: False

use rule * from diff_analysis

rule rawdata_fetch:
    input:
        rawpath = config["path"]["rawdatapath"],
        confirm_file = confirm_file
    output:
        project_path = expand("{{savepath}}/1.Project_information/{registrationfile}",registrationfile = registrationfile)
    shell:
        '''
        cp '{input.confirm_file}' '{output.project_path}'
        '''

rule rawdata_information:
    input:
        config["path"]["diff_xls"],
        expand("{projectpath}/result/2.Qualitative/{my_dir}",
               projectpath=config["path"]["projectpath"],
               my_dir = ["Annotation"]),
    output:
        "{savepath}/overview_information.xlsx"
    params:
        rawdatapath = config["path"]["rawdatapath"],
        projectpath = config["path"]["projectpath"]
    shell:
        '''
        /opt/conda/bin/Rscript -e "lmbio::over_table( \
        rawdatapath = '{params.rawdatapath}', \
        projectpath = '{params.projectpath}', \
        savepath = '{output}')"
        '''
        
rule protein_motif:
    input:
        config["path"]["diff_xls"]
    output:
        directory("{savepath}/Motif/{compare}")
    params:
        rawdatapath = config["path"]["rawdatapath"],
        projectpath = config["path"]["projectpath"]
    shell:
        '''
        runmotif \
        -if {input}/diff-data-{wildcards.compare}.xls \
        -g {wildcards.compare} \
        -sp '{wildcards.savepath}/Motif/'
        '''


rule rawdata_plot:
    input:
        rawpath = config["path"]["rawdatapath"],
        bg = config["params"]["enrich"]["database"]
    output:
        directory("{savepath}/Supplementary"),
        anno_path = directory("{savepath}/2.Qualitative/Annotation"),
        stat_path = directory("{savepath}/2.Qualitative/Statistics"),
        Identified_path = directory("{savepath}/2.Qualitative/Identified"),
        out_plot = directory(expand("{{savepath}}/3.Reliable_results/{plot_dir}",plot_dir=["Boxplot","Densityplot","Sampletreeplot"]))
    params:
        fontfamily = config["params"]["plot"]["family"],
        classfile = config["path"]["classtypefile"]
    shell:
        '''
        unionplot -ip '{input.rawpath}/' \
        -sh '{input.bg}/' \
        -sp '{wildcards.savepath}/2.Qualitative/Statistics/' \
        -sg '{wildcards.savepath}/2.Qualitative/Identified/' \
        -ap '{wildcards.savepath}/2.Qualitative/Annotation/' \
        -cp '{wildcards.savepath}/3.Reliable_results/' \
        --classfile {params.classfile} \
        --fontfamily '{params.fontfamily}'
        cp '{input.rawpath}/绘图数据.xlsx' '{wildcards.savepath}/3.Reliable_results'
        cp -r /data/hstore3/public/propip/Supplementary '{wildcards.savepath}/'
        if [ -e '{input.rawpath}/Protein quantitation.xlsx' ] && [ ! -d '{input.rawpath}/Label set1' ]; \
        then \
            cp '{input.rawpath}/Protein quantitation.xlsx' '{wildcards.savepath}/2.Qualitative/'; \
        fi
        if [ -d '{input.rawpath}/Label set1' ]; \
        then \
            cp -r {input.rawpath}/Label* '{wildcards.savepath}/2.Qualitative/'; \
        fi
        if [ -e '{input.rawpath}/Peptides.xlsx' ]; \
        then \
            cp '{input.rawpath}/Peptides.xlsx' '{wildcards.savepath}/2.Qualitative/'; \
        fi
        ptm=`find {input.rawpath}/ -name *Sites.xlsx | wc -l`
        if [ $ptm -eq 1 ]; \
        then \
            cp {input.rawpath}/*Sites.xlsx '{wildcards.savepath}/2.Qualitative/'; \
        fi
        if [ -d '{input.rawpath}/Label*' ]; \
        then \
            cp -r '{input.rawpath}/Label*' '{wildcards.savepath}/2.Qualitative/'; \
        fi
        if [ -e {input.rawpath}/DDA* ]; \
        then \
            cp {input.rawpath}/DDA* '{wildcards.savepath}/2.Qualitative/'; \
        fi
        if [ -e {input.rawpath}/QC.docx ]; \
        then \
            cp {input.rawpath}/QC.docx '{wildcards.savepath}/Supplementary/'; \
        fi
        if [ -e {input.rawpath}/*.fa* ]; \
        then \
            cp {input.rawpath}/*.fa* '{wildcards.savepath}/2.Qualitative/'; \
        fi
        '''
        
rule protein_result_pca:
    input:
        expand("{projectpath}/result/3.Reliable_results/PCA/none_label/PCA-{type}-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               compare=[i for i in compare if compare[i]["samplenum"]>2],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
               type=["score","3Dscore"]),
        expand("{projectpath}/result/3.Reliable_results/PCA/label/PCA-{type}-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               compare=[i for i in compare if compare[i]["samplenum"]>2],
               imagetype=config["params"]["imagetype"],
               type=["score"]),
        expand("{projectpath}/result/3.Reliable_results/PCA/none_label/PCA-{type}-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               compare=[i for i in compare_all if compare_all[i]["samplenum"]>2],
               imagetype=config["params"]["imagetype"]+config["params"]["imagetype_other"],
               type=["score","3Dscore"]),
        expand("{projectpath}/result/3.Reliable_results/PCA/label/PCA-{type}-{compare}.{imagetype}",
               projectpath=config["path"]["projectpath"],
               compare=[i for i in compare_all if compare_all[i]["samplenum"]>2],
               imagetype=config["params"]["imagetype"],
               type=["score"])
               
rule protein_result_oecloud:
  input:
    expand("{projectpath}/result/OECloud_tools/background.oecloud",
           projectpath=config["path"]["projectpath"]),
    expand("{projectpath}/result/OECloud_tools/data_for_oecloud.xlsx",
           projectpath=config["path"]["projectpath"]),
    expand("{projectpath}/result/OECloud_tools/diff-data.oecloud",
           projectpath=config["path"]["projectpath"])

rule protein_result_qc:
  input:
    expand("{projectpath}/result/2.Qualitative/{my_dir}",
           projectpath=config["path"]["projectpath"],
           my_dir = ["Annotation","Identified","Statistics"]),
    expand("{projectpath}/result/3.Reliable_results/{plot_dir}",
           projectpath=config["path"]["projectpath"],
           plot_dir=["Boxplot","Densityplot","Sampletreeplot"])

rule protein_result_diff_cor:
  input:
    expand("{projectpath}/result/4.Different_expression/相关性分析/Correlation_top20/Correlation_top20-{compare}.{imagetype}",
           projectpath=config["path"]["projectpath"],
           compare=[i for i in compare if compare[i]["samplenum"]>5],
           imagetype=config["params"]["imagetype"])

if config["params"]["omic"] == "PF":
  rule protein_result_diff_cor_1:
    input:
      rules.protein_result_diff_cor.input
else:
  rule protein_result_diff_cor_1:
    input:
    
if config["params"]["MStype"] in ["PT","PL"]:
  rule protein_result_diff_motif:
    input:
      expand("{projectpath}/result/4.Different_expression/Motif/{compare}",
             projectpath=config["path"]["projectpath"],
             compare=compare)
else:
  rule protein_result_diff_motif:
    input:

rule protein_result_diff_gsea:
  input:
    expand("{projectpath}/result/7.GSEA/{type}/{compare}/",
           projectpath=config["path"]["projectpath"],
           compare=[i for i in compare if compare[i]["groupnum"]==2],
           type=gseatype)
           
rule protein_result_diff_ppi:
  input:
    expand("{projectpath}/result/6.{type}/{compare}/",
           projectpath=config["path"]["projectpath"],
           compare=compare,
           type=ppitype)

rule protein_result_diff_enrich:
  input:
    expand("{projectpath}/result/5.Enrichment/{type}/{compare}/",
           projectpath=config["path"]["projectpath"],
           compare=compare,
           type=enrichtype)


rule protein_result_diff_function:
  input:
    rules.protein_result_diff_enrich.input,
    rules.protein_result_diff_ppi.input,
    rules.protein_result_diff_gsea.input

rule protein_result_diff_xlsx:
  input:          
    expand("{projectpath}/result/4.Different_expression/Foldchange",
           projectpath=config["path"]["projectpath"]),
    expand("{projectpath}/result/4.Different_expression/表达矩阵.xlsx",
           projectpath=config["path"]["projectpath"]),
    expand("{projectpath}/result/4.Different_expression/差异表达矩阵.xlsx",
           projectpath=config["path"]["projectpath"]),
    expand("{projectpath}/result/4.Different_expression/差异表达矩阵(未筛选).xlsx",
           projectpath=config["path"]["projectpath"])

rule protein_result_diff_map:
  input:   
    expand("{projectpath}/result/4.Different_expression/Heatmap/cluster_none/heatmap-{compare}.{imagetype}",
           projectpath=config["path"]["projectpath"],
           compare=compare,
           imagetype=config["params"]["imagetype"]),
    expand("{projectpath}/result/4.Different_expression/Heatmap/cluster_samples/heatmap-{compare}.{imagetype}",
           projectpath=config["path"]["projectpath"],
           compare=compare,
           imagetype=config["params"]["imagetype"]),
    expand("{projectpath}/result/4.Different_expression/Volcano/volcano-{compare}.{imagetype}",
           projectpath=config["path"]["projectpath"],
           compare=[i for i in compare if compare[i]["groupnum"]==2 and compare[i]["minsamplenum"]>1],
           imagetype=config["params"]["imagetype"]),
    expand("{projectpath}/result/4.Different_expression/Venn/{venn}.{imagetype}",
           projectpath=config["path"]["projectpath"],
           venn=[] if len(config["params"]["all_compare"]) == 1 else ["Venn"],
           imagetype=config["params"]["imagetype"])

rule protein_result_diff:
  input:
    rules.protein_result_diff_xlsx.input,
    rules.protein_result_diff_map.input,
    rules.protein_result_diff_cor_1.input,
    rules.protein_result_diff_function.input,
    rules.protein_result_diff_motif.input

rule protein_result:
  input:
    expand("{projectpath}/result/1.Project_information/{registrationfile}",
           projectpath=config["path"]["projectpath"],
           registrationfile=[registrationfile,"overview_information.xlsx"]),
    rules.protein_result_oecloud.input,
    rules.protein_result_qc.input,
    rules.protein_result_pca.input,
    rules.protein_result_diff.input,
    expand("{projectpath}/result/Supplementary",
           projectpath=config["path"]["projectpath"])

if not os.path.exists(config["path"]["config_savepath"]): 
    os.makedirs(config["path"]["config_savepath"])
with open("{path}/config_protein_analysis.yaml".format(path=config["path"]["config_savepath"]),"w",encoding="utf-8") as f:
    yaml.dump(config,f,allow_unicode=True,sort_keys=False)
