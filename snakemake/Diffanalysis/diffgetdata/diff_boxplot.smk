rule diff_getdata_violin_top50:
    input:
        rdspath=config["path"]["diff_filter_rds_savepath"]+"/diff_filter-{compare}.rds"
    output:
        "{savepath}/Violin_top50-{compare}.xls"
    threads: 1
    shell:
        '''
        flow_diff_getheatmapdata \
        --rdspath '{input.rdspath}' \
        --savepath '{output}' \
        --top 50
        '''

rule diff_plot_violin_top50:
    input:
        "{savepath}/Violin_top50/Violin_top50-{compare}.xls"
    output:
        directory("{savepath}/Violin_top50/{compare}")
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        classfile = config["path"]["classtypefile"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_common_boxplot \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/Violin_top50/{wildcards.compare}' \
        --mapname 'Violin-{wildcards.compare}' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height} \
        --classfile {params.classfile}
        '''

rule diff_getdata_boxplot_top50:
    input:
        rdspath=config["path"]["diff_filter_rds_savepath"]+"/diff_filter-{compare}.rds"
    output:
        "{savepath}/Boxplot_top50-{compare}.xls"
    threads: 1
    shell:
        '''
        flow_diff_getheatmapdata \
        --rdspath '{input.rdspath}' \
        --savepath '{output}' \
        --top 50
        '''

rule diff_plot_boxplot_top50:
    input:
        "{savepath}/Boxplot_top50/Boxplot_top50-{compare}.xls"
    output:
        directory("{savepath}/Boxplot_top50/{compare}")
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        classfile = config["path"]["classtypefile"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_common_boxplot2 \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/Boxplot_top50/{wildcards.compare}' \
        --mapname 'Boxplot-{wildcards.compare}' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height} \
        --classfile {params.classfile}
        '''