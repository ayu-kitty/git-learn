import os

projectname = config["projectname"] 
mode = config["mode"]
analysis_id = config["analysis_id"]

rule all:
    input:
        projectname+"/质谱流式单细胞蛋白组分析报告.html"

rule report:
    input:
        "./项目登记单.xlsx",
        projectname+"/1.QC/bar_info.xlsx",
        directory(projectname+"/2.SPADE"),
        projectname+"/3.Cluster/heatmap_info.xlsx",
        projectname+"/4.Reducedim/reduced_dim_info.xlsx",
        projectname+"/5.Celltype/cellanno.xlsx"
    output:
        "{path}/质谱流式单细胞蛋白组分析报告.html"
    params:
        mode=mode
    shell:
        "cytof_2.9_report -s {wildcards.path} -m {params.mode}"

# rule get_project_info:
#     input:
#     output:
#         "{path}/项目登记单.xlsx"
#     params:
#         analysis_id=analysis_id
#     shell:
#         "cytof_1.1_GetAnalystInfo -ad {params.analysis_id} -sp {wildcards.path}"
    
rule fcs_to_rds:
    input:
    output:
        "rawdata/flowSet_object.rds"
    shell:
        "cytof_1.2_readfcs"

rule fcs_to_spade:
    input:
        "rawdata/flowSet_object.rds"
    output:
        directory(projectname+"/2.SPADE")
    shell:
        "cytof_1.3_spade -sp {output}"


rule metadata_panel:
    input:
    output:
        "rawdata/metadata.xlsx",
        "rawdata/panel.xlsx"
    shell:
        "cytof_1.4_metadata_panel"

rule create_sce:
    input:
        "rawdata/metadata.xlsx",
        "rawdata/panel.xlsx",
        "rawdata/flowSet_object.rds"
    output:
        "rawdata/flowSet_object_sce.rds"
    shell:
        "cytof_2.1_create_sce"

rule qc_plot:
    input:
        "rawdata/flowSet_object_sce.rds"
    output:
        "{path}/bar_info.xlsx"
    shell:
        "cytof_2.2_sce_qcplot -sp {wildcards.path}"

rule sce_cluster:
    input:
        "rawdata/flowSet_object_sce.rds"
    output:
        "rawdata/flowSet_object_sce_cluster.rds"
    shell:
        "cytof_2.3_sce_cluster"

rule sce_cluster_plot:
    input:
        "rawdata/flowSet_object_sce_cluster.rds"
    output:
        "{path}/heatmap_info.xlsx"
    shell:
        "cytof_2.4_sce_cluster_plot -sp {wildcards.path}"

rule sce_reducedim:
    input:
        "rawdata/flowSet_object_sce_cluster.rds"
    output:
        "rawdata/flowSet_object_sce_cluster_reducedDim.rds"
    shell:
        "cytof_2.5_sce_reducedim"

rule dim_plot:
    input:
        "rawdata/flowSet_object_sce_cluster_reducedDim.rds"
    output:
        "{path}/reduced_dim_info.xlsx"
    shell:
        "cytof_2.6_sce_reducedim_plot -sp {wildcards.path}"

rule sce_anno:
    input:
        cell="rawdata/celltype.xlsx",
        heatmap=projectname+"/3.Cluster/heatmap_info.xlsx",
        rds="rawdata/flowSet_object_sce_cluster_reducedDim.rds"
    output:
        "{path}/cellanno.xlsx"
    params:
        diffpath=projectname+"/6.Diffanalysis"
    shell:
        '''
        cytof_2.7_sce_anno -sp {wildcards.path} -clust {input.heatmap} -cell {input.cell}
        cytof_2.8_sce_diffprops -clust meta20 -sp {params.diffpath}
        cytof_2.8_sce_diffprops -clust merging1 -sp {params.diffpath}
        '''
