
# enrich_db
module enrich_db:
    snakefile: "enrich_db.smk"
    config: config

use rule * from enrich_db

enclass = {"G":"GO","K":"KEGG","I":"InterPro","R":"Reactome","W":"Wikipathways","ppi":"PPI_network"}
if config["params"]["omic"] == "M" or config["params"]["omic"] == "ML":
  gseatype = ["KEGG"]
else:
  gseatype = ["KEGG","GO"]

if "gseatype" in config["params"]:
  gseatype2 = [enclass[key] for key in config["params"]["gseatype"].split(",")]
  gseatype = list(set(gseatype2) & set(gseatype))
if "gseadatapath" not in config["params"]["enrich"]:
  config["params"]["enrich"]["gseadatapath"]="rds"
if "gseamin" not in config["params"]["enrich"] or config["params"]["enrich"]["gseamin"]==0:
  if config["params"]["omic"] == "M" or config["params"]["omic"] == "ML":
    config["params"]["enrich"]["gseamin"]=5
  else:
    config["params"]["enrich"]["gseamin"]=15
if "graph" not in config["params"]["enrich"]:
  config["params"]["enrich"]["graph"]=20
if "termlist" not in config["params"]["enrich"]:
  config["params"]["enrich"]["termlist"]=""
compare=config["params"]["all_compare"]

rule gsea_result:
    input:
        expand("{projectpath}/GSEA/{type}/{compare}/",
               projectpath=config["path"]["projectpath"],
               compare=[i for i in compare if compare[i]["groupnum"]==2],
               type=gseatype)

if config["params"]["enrich"]["gseadatapath"] == "rds":
  rule gsea_getdata:
      input:
          gseafile = expand("{diff_filter_rds_savepath}/diff_filter-{{compare}}.rds",
                            diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"]),
          classfile = config["path"]["classfile"]
      output:
          "{savepath}/{compare}-gsea-data.xls",
          "{savepath}/{compare}-gsea-data.cls"
      params:
          rdspath = config["path"]["diff_filter_rds_savepath"],
          omic = config["params"]["omic"]
      threads: 1
      shell:
          '''
          makegsea \
          --inputfile '{input.gseafile}' \
          --outputpath '{wildcards.savepath}' \
          --comparename '{wildcards.compare}'  \
          --groupfile '{input.classfile}'\
          --omic "{params.omic}"
        '''
else:
    rule gsea_getdata:
      input:
          gseafile = config["params"]["enrich"]["gseadatapath"],
          classfile = config["path"]["classfile"]
      output:
          "{savepath}/{compare}-gsea-data.xls",
          "{savepath}/{compare}-gsea-data.cls"
      params:
          rdspath = config["path"]["diff_filter_rds_savepath"],
          omic = config["params"]["omic"]
      threads: 1
      shell:
          '''
          makegsea \
          --inputfile '{input.gseafile}' \
          --outputpath '{wildcards.savepath}' \
          --comparename '{wildcards.compare}'  \
          --groupfile '{input.classfile}'\
          --omic "{params.omic}"
        '''

rule gsea_go:
    input:
        gseafile = config["path"]["enrich_xls_savepath"]+"/{compare}-gsea-data.xls",
        clsfile = config["path"]["enrich_xls_savepath"]+"/{compare}-gsea-data.cls",
        db = config["params"]["enrich"]["database"]
    output:
        directory("{savepath,.*GSEA}/GO/{compare}")
    params:
        fontfamily = config["params"]["plot"]["family"],
        min = config["params"]["enrich"]["gseamin"],
        graph = config["params"]["enrich"]["graph"],
        termlist = config["params"]["enrich"]["termlist"]
    shell:
        '''
        oe gsea \
        gsea '{input.gseafile}' \
        -g '{input.db}/go.gmt' \
        -c '{input.clsfile}' \
        -m signal_to_noise \
        -o 'go' \
        -min {params.min} -graph {params.graph} {params.termlist} >/dev/null 2>&1
        /opt/conda/bin/Rscript -e "lmbio::heatmap_in_gsea('{wildcards.compare}-gsea-data-go/go')"
        mv '{wildcards.compare}-gsea-data-go/go' '{output}'
        rm -rf '{wildcards.compare}-gsea-data-go'
        '''


rule gsea_kegg:
    input:
        gseafile = config["path"]["enrich_xls_savepath"]+"/{compare}-gsea-data.xls",
        clsfile = config["path"]["enrich_xls_savepath"]+"/{compare}-gsea-data.cls",
        db = config["params"]["enrich"]["database"]
    output:
        directory("{savepath,.*GSEA}/KEGG/{compare}")
    params:
        fontfamily = config["params"]["plot"]["family"],
        min = config["params"]["enrich"]["gseamin"],
        graph = config["params"]["enrich"]["graph"],
        termlist = config["params"]["enrich"]["termlist"]
    shell:
        '''
        oe gsea \
        gsea '{input.gseafile}' \
        -g '{input.db}/kegg.gmt' \
        -c '{input.clsfile}' \
        -m signal_to_noise \
        -o 'kegg' \
        -min {params.min} -graph {params.graph} {params.termlist} >/dev/null 2>&1
        /opt/conda/bin/Rscript -e "lmbio::heatmap_in_gsea('{wildcards.compare}-gsea-data-kegg/kegg')"
        mv '{wildcards.compare}-gsea-data-kegg/kegg' '{output}'
        rm -rf '{wildcards.compare}-gsea-data-kegg'
        '''

rule gsea_kegg2:
    input:
        gseafile = config["path"]["enrich_xls_savepath"]+"/{compare}-gsea-data.xls",
        clsfile = config["path"]["enrich_xls_savepath"]+"/{compare}-gsea-data.cls",
        db = config["params"]["enrich"]["database"]
    output:
        directory("{savepath,.*GSEA}/{compare,((?!/).)*}")
    params:
        fontfamily = config["params"]["plot"]["family"],
        min = config["params"]["enrich"]["gseamin"],
        graph = config["params"]["enrich"]["graph"],
        termlist = config["params"]["enrich"]["termlist"]
    shell:
        '''
        oe gsea \
        gsea '{input.gseafile}' \
        -g '{input.db}/kegg.gmt' \
        -c '{input.clsfile}' \
        -m signal_to_noise \
        -o 'kegg' \
        -min {params.min} -graph {params.graph} {params.termlist} >/dev/null 2>&1
        /opt/conda/bin/Rscript -e "lmbio::heatmap_in_gsea('{wildcards.compare}-gsea-data-kegg/kegg')"
        mv '{wildcards.compare}-gsea-data-kegg/kegg' '{output}'
        rm -rf '{wildcards.compare}-gsea-data-kegg'
        '''
