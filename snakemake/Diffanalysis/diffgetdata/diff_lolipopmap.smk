rule diff_plot_lolipopmap:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{{compare}}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"])
    output:
        "{savepath}/lolipopmap-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)}"
    params:
        rdspath = config["path"]["diff_filter_rds_savepath"],
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"]
    threads: 1
    shell:
        '''
        map_common_lolipopmap \
        --filename '{input}' \
        --savepath '{wildcards.savepath}' \
        --mapname 'lolipopmap-{wildcards.compare}' \
        --imagetype '{wildcards.imagetype}' \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''
