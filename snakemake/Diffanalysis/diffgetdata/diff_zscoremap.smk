rule diff_getdata_zscore_top20:
    input:
        rdspath=config["path"]["diff_filter_rds_savepath"]+"/diff_filter-{compare}.rds"
    output:
        "{savepath}/Zscore_top20-{compare}.xls"
    threads: 1
    shell:
        '''
        flow_diff_getheatmapdata \
        --rdspath '{input.rdspath}' \
        --savepath '{output}' \
        --top 20
        '''

rule diff_plot_zscore_top20:
    input:
        "{savepath}/Zscore_top20-{compare}.xls"
    output:
        "{savepath}/Zscore_top20-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        classfile = config["path"]["classtypefile"],
        scale = lambda wildcards: "'none' --log" if config["compare"][wildcards.compare]["samplenum"] == 2 else "row"
    threads: 1
    shell:
        '''
        map_common_zscore \
        --filename '{input}' \
        --savepath '{wildcards.savepath}' \
        --mapname 'Zscore_top20-{wildcards.compare}' \
        --imagetype '{wildcards.imagetype}' \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height} \
        --classfile {params.classfile}
        '''
