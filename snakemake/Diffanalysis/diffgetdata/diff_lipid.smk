rule diff_plot_lipid_classfication_bubbleplot:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{{compare}}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"])
    output:
        "{savepath}/classfication_bubbleplot-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)}"
    params:
        rdspath = config["path"]["diff_filter_rds_savepath"],
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"]
    threads: 1
    shell:
        '''
        map_lipid_classfication_bubbleplot \
        --filename '{input}' \
        --savepath '{wildcards.savepath}' \
        --mapname 'classfication_bubbleplot-{wildcards.compare}' \
        --imagetype '{wildcards.imagetype}' \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''

rule diff_plot_lipid_carbonheatmap:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{{compare}}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"])
    output:
        expand("{{savepath}}/carbonheatmap-{{compare}}.{imagetype}",
               imagetype=config["params"]["imagetype"])
    params:
        rdspath = config["path"]["diff_filter_rds_savepath"],
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_lipid_carbonheatmap \
        --filename '{input}' \
        --savepath '{wildcards.savepath}' \
        --mapname 'carbonheatmap-{wildcards.compare}' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''

rule diff_plot_lipid_unsaturationheatmap:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{{compare}}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"])
    output:
        expand("{{savepath}}/unsaturationheatmap-{{compare}}.{imagetype}",
               imagetype=config["params"]["imagetype"])
    params:
        rdspath = config["path"]["diff_filter_rds_savepath"],
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_lipid_unsaturationheatmap \
        --filename '{input}' \
        --savepath '{wildcards.savepath}' \
        --mapname 'unsaturationheatmap-{wildcards.compare}' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''
        
rule diff_plot_lipid_point:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{{compare}}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"])
    output:
        directory("{savepath}/气泡分布图/{compare}")
    params:
        rdspath = config["path"]["diff_filter_rds_savepath"],
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_lipid_point \
        --filename '{input}' \
        --savepath '{output}' \
        --mapname 'bubbleplot-{wildcards.compare}' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''