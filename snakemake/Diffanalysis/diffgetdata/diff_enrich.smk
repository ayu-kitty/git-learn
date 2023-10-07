
enclass = {"G":"GO","K":"KEGG","I":"InterPro","R":"Reactome","W":"Wikipathways","ppi":"PPI_network"}
if config["params"]["omic"] == "M" or config["params"]["omic"] == "ML":
  if config["params"]["enrich"]["org"] in ["hsa","mmu","rno"]:
    enrichtype = ["KEGG","KEGG_map","Reactome"]
  else:
    enrichtype = ["KEGG","KEGG_map"]
else:
  if config["params"]["enrich"]["org"] in ["hsa","mmu","rno"]:
    enrichtype = ["KEGG","KEGG_map","GO","InterPro","Wikipathways","Reactome","PPI_network"]
  else:
    enrichtype = ["KEGG","KEGG_map","GO","InterPro","PPI_network"]

if "enrichtype" in config["params"]:
  enrichtype2 = [enclass[key] for key in config["params"]["enrichtype"].split(",")]
  enrichtype = list(set(enrichtype2) & set(enrichtype))
  print(enrichtype)
if "ppitopn" not in config["params"]["enrich"]:
  config["params"]["enrich"]["ppitopn"]=25
  
if "projectppipath" not in config["path"]:
  config["path"]["projectppipath"]="PPI_network"
  
# enrich_db
module enrich_db:
    snakefile: "enrich_db.smk"
    config: config

use rule * from enrich_db

rule enrich_result:
    input:
        expand("{projectpath}/{type}/{compare}/",
               projectpath=config["path"]["projectpath"],
               compare=config["params"]["all_compare"],
               type=[i for i in enrichtype if i != "PPI_network" ]),
        expand("{projectppipath}/{type}/{compare}/",
               projectppipath=config["path"]["projectppipath"],
               compare=config["params"]["all_compare"],
               type=[i for i in enrichtype if i == "PPI_network" ])

rule enrich_diff_getdata_filter:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{{compare}}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"])
    output:
        "{savepath}/{compare}-diff-data.xls"
    params:
        rdspath = config["path"]["diff_filter_rds_savepath"],
        omic = config["params"]["omic"]
    threads: 1
    shell:
        '''
        makediff \
        --inputfile '{input}' \
        --outputpath '{wildcards.savepath}' \
        --comparename '{wildcards.compare}'  \
        --omic "{params.omic}"
        '''

rule enrich_diff_getdata_filter2:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{{compare}}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"])
    output:
        "{savepath}/{compare}-diff-data_chebi.xls"
    params:
        rdspath = config["path"]["diff_filter_rds_savepath"],
        omic = config["params"]["omic"]
    threads: 1
    shell:
        '''
        makediff \
        --inputfile '{input}' \
        --outputpath '{wildcards.savepath}' \
        --comparename '{wildcards.compare}'  \
        --omic "{params.omic}" \
        --chebi
        '''

rule enrich_kegg:
    input:
        comparefile=config["path"]["enrich_xls_savepath"]+"/{compare}-diff-data.xls",
        dbpath=config["params"]["enrich"]["database"]
    output:
        kegg=directory("{savepath,((?!GSEA).)*}/KEGG/{compare}/")
    threads: 1
    params:
        path=config["path"]["enrich_xls_savepath"],
        org=config["params"]["enrich"]["org"],
        db=config["params"]["enrich"]["database"],
        kd=config["params"]["enrich"]["keggdb"],
        fontfamily = config["params"]["plot"]["family"],
        projectpath = config["path"]["projectpath"]
    shell:
        '''
        path=$(pwd)
        cd "{params.path}"
        enrichpip -t K \
        -if "{wildcards.compare}-diff-data.xls" \
        -db "$path/{params.db}" \
        -kd "{params.kd}" \
        -o "{params.org}" \
        -s "$path/{wildcards.savepath}/" \
        -fa "{params.fontfamily}"
        cd "$path"
        
        while [ -f "文件压缩中.log" ]
        do
          # echo 1
          sleep 10
        done
        
        echo "压缩中" > 文件压缩中.log
        
        # if [ -d '{wildcards.savepath}/KEGG_map/{wildcards.compare}' ]; then
        #   if [ ! $(find '{wildcards.savepath}/KEGG_map' -name KEGG_map.tar) ]; then
        #     tar -cf '{wildcards.savepath}/KEGG_map/KEGG_map.tar' --directory='{wildcards.savepath}/KEGG_map' {wildcards.compare}/.
        #   else
        #     tar -rf '{wildcards.savepath}/KEGG_map/KEGG_map.tar' --directory='{wildcards.savepath}/KEGG_map' {wildcards.compare}/.
        #   fi
        #     rm -rf '{wildcards.savepath}/KEGG_map/{wildcards.compare}'
        # fi
        
        if [ -f '{output.kegg}/enrichment-kegg-{wildcards.compare}-Total.xls' ]; then
          if [ ! $(find '{params.projectpath}' -name diff-data.oecloud) ]; then
            echo "diff-data.oecloud is null"
          else
            tar -rf "$(find '{params.projectpath}' -name diff-data.oecloud)" --directory='{output.kegg}' 'enrichment-kegg-{wildcards.compare}-Total.xls'
          fi
        fi
        
        rm 文件压缩中.log
        '''

rule enrich_go:
    input:
        config["path"]["enrich_xls_savepath"]+"/{compare}-diff-data.xls",
        config["params"]["enrich"]["database"]
    output:
        directory("{savepath,((?!GSEA).)*}/GO/{compare}/"),
    threads: 1
    priority: 1
    params:
        path=config["path"]["enrich_xls_savepath"],
        org=config["params"]["enrich"]["org"],
        db=config["params"]["enrich"]["database"],
        kd=config["params"]["enrich"]["keggdb"],
        fontfamily = config["params"]["plot"]["family"]
    shell:
        '''
        path=$(pwd)
        cd "{params.path}"
        enrichpip -t G \
        -if "{wildcards.compare}-diff-data.xls" \
        -db "$path/{params.db}" \
        -kd "{params.kd}" \
        -o "{params.org}" \
        -s "$path/{wildcards.savepath}/" \
        -fa "{params.fontfamily}"
        '''

rule enrich_InterPro:
    input:
        config["path"]["enrich_xls_savepath"]+"/{compare}-diff-data.xls",
        config["params"]["enrich"]["database"]
    output:
        directory("{savepath}/InterPro/{compare}/"),
    threads: 1
    params:
        path=config["path"]["enrich_xls_savepath"],
        org=config["params"]["enrich"]["org"],
        db=config["params"]["enrich"]["database"],
        kd=config["params"]["enrich"]["keggdb"],
        fontfamily = config["params"]["plot"]["family"]
    shell:
        '''
        path=$(pwd)
        cd "{params.path}"
        enrichpip -t I \
        -if "{wildcards.compare}-diff-data.xls" \
        -db "$path/{params.db}" \
        -kd "{params.kd}" \
        -o "{params.org}" \
        -s "$path/{wildcards.savepath}/" \
        -fa "{params.fontfamily}"
        '''

if config["params"]["omic"] == "M" or config["params"]["omic"] == "ML":
  rule enrich_Reactome:
      input:
          config["path"]["enrich_xls_savepath"]+"/{compare}-diff-data_chebi.xls",
          config["params"]["enrich"]["database"]
      output:
          directory("{savepath}/Reactome/{compare}/"),
      threads: 1
      params:
          path=config["path"]["enrich_xls_savepath"],
          org=config["params"]["enrich"]["org"],
          db=config["params"]["enrich"]["database"],
          kd=config["params"]["enrich"]["keggdb"],
          fontfamily = config["params"]["plot"]["family"]
      shell:
          '''
          path=$(pwd)
          cd "{params.path}"
          enrichpip -t R \
          -if "{wildcards.compare}-diff-data_chebi.xls" \
          -db "$path/{params.db}" \
          -kd "{params.kd}" \
          -o "{params.org}" \
          -s "$path/{wildcards.savepath}/" \
          -fa "{params.fontfamily}"
          '''
else:
  rule enrich_Reactome:
      input:
          config["path"]["enrich_xls_savepath"]+"/{compare}-diff-data.xls",
          config["params"]["enrich"]["database"]
      output:
          directory("{savepath}/Reactome/{compare}/"),
      threads: 1
      params:
          path=config["path"]["enrich_xls_savepath"],
          org=config["params"]["enrich"]["org"],
          db=config["params"]["enrich"]["database"],
          kd=config["params"]["enrich"]["keggdb"],
          fontfamily = config["params"]["plot"]["family"]
      shell:
          '''
          path=$(pwd)
          cd "{params.path}"
          enrichpip -t R \
          -if "{wildcards.compare}-diff-data.xls" \
          -db "$path/{params.db}" \
          -kd "{params.kd}" \
          -o "{params.org}" \
          -s "$path/{wildcards.savepath}/" \
          -fa "{params.fontfamily}"
        '''

rule enrich_Wikipathways:
    input:
        config["path"]["enrich_xls_savepath"]+"/{compare}-diff-data.xls",
        config["params"]["enrich"]["database"]
    output:
        directory("{savepath}/Wikipathways/{compare}/"),
    threads: 1
    params:
        path=config["path"]["enrich_xls_savepath"],
        org=config["params"]["enrich"]["org"],
        db=config["params"]["enrich"]["database"],
        kd=config["params"]["enrich"]["keggdb"],
        fontfamily = config["params"]["plot"]["family"]
    shell:
        '''
        path=$(pwd)
        cd "{params.path}"
        enrichpip -t W \
        -if "{wildcards.compare}-diff-data.xls" \
        -db "$path/{params.db}" \
        -kd "{params.kd}" \
        -o "{params.org}" \
        -s "$path/{wildcards.savepath}/" \
        -fa "{params.fontfamily}"
        '''

rule enrich_ppi:
    input:
        config["path"]["enrich_xls_savepath"]+"/{compare}-diff-data.xls",
        config["params"]["enrich"]["database"]
    output:
        directory("{savepath,.*PPI_network}/{compare}/")
    params:
        path=config["path"]["enrich_xls_savepath"],
        db=config["params"]["enrich"]["database"],
        topn=config["params"]["enrich"]["ppitopn"],
        fontfamily = config["params"]["plot"]["family"]
    shell:
        '''
        path=$(pwd)
        cd "{params.path}"
        proppi \
        -if "{wildcards.compare}-diff-data.xls" \
        -bp "$path/{params.db}" \
        -fa "{params.fontfamily}" \
        -s "$path/{wildcards.savepath}/" \
        --number {params.topn}
        '''

