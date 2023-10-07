
rule diff_getdata_heatmap:
    input:
        rdspath=config["path"]["diff_filter_rds_savepath"]+"/diff_filter-{compare}.rds"
    output:
        "{savepath}/heatmap-{compare}.xls"
    threads: 1
    priority: 1
    shell:
        '''
        flow_diff_getheatmapdata \
        --rdspath '{input.rdspath}' \
        --savepath '{output}'
        '''
        
rule diff_getdata_heatmap_scale:
    input:
        rdspath=config["path"]["diff_filter_rds_savepath"]+"/diff_filter-{compare}.rds"
    output:
        "{savepath}/heatmap_scale-{compare}.xls"
    threads: 1
    priority: 1
    shell:
        '''
        flow_diff_getheatmapdata \
        --rdspath '{input.rdspath}' \
        --savepath '{output}' \
        --scale
        '''
        
rule diff_plot_heatmap:
    input:
        "{savepath}/heatmap-{compare}.xls"
    output:
        "{savepath,((?!cluster).)*}/heatmap-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        classfile = config["path"]["classtypefile"],
        scale = lambda wildcards: "'none' --log" if config["compare"][wildcards.compare]["samplenum"] == 2 else "row"
    threads: 1
    priority: 1
    shell:
        '''
        map_common_heatmap \
        --filename '{input}' \
        --savepath '{wildcards.savepath}' \
        --mapname 'heatmap-{wildcards.compare}' \
        --imagetype '{wildcards.imagetype}' \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height} \
        --show_rownames \
        --show_colnames \
        --cluster_rows \
        --cluster_cols \
        --scale {params.scale} \
        --classfile {params.classfile} \
        --colgroup 1
        '''
        
rule diff_plot_heatmap_cluster_none:
    input:
        "{savepath}/heatmap-{compare}.xls"
    output:
        "{savepath}/cluster_none/heatmap-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        classfile = config["path"]["classtypefile"],
        scale = lambda wildcards: "'none' --log" if config["compare"][wildcards.compare]["samplenum"] == 2 else "row"
    threads: 1
    priority: 1
    shell:
        '''
        map_common_heatmap \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/cluster_none/' \
        --mapname 'heatmap-{wildcards.compare}' \
        --imagetype '{wildcards.imagetype}' \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height} \
        --show_rownames \
        --show_colnames \
        --cluster_rows \
        --scale {params.scale} \
        --classfile {params.classfile} \
        --colgroup 1
        '''

rule diff_plot_heatmap_cluster_samples:
    input:
        "{savepath}/heatmap-{compare}.xls"
    output:
        "{savepath}/cluster_samples/heatmap-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        classfile = config["path"]["classtypefile"],
        scale = lambda wildcards: "'none' --log" if config["compare"][wildcards.compare]["samplenum"] == 2 else "row"
    threads: 1
    priority: 1
    shell:
        '''
        map_common_heatmap \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/cluster_samples/' \
        --mapname 'heatmap-{wildcards.compare}' \
        --imagetype '{wildcards.imagetype}' \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height} \
        --show_rownames \
        --show_colnames \
        --cluster_rows \
        --cluster_cols \
        --scale {params.scale} \
        --classfile {params.classfile} \
        --colgroup 1
        '''
     
rule diff_getdata_heatmap_top50:
    input:
        rdspath=config["path"]["diff_filter_rds_savepath"]+"/diff_filter-{compare}.rds"
    output:
        "{savepath}/heatmap_top50-{compare}.xls"
    threads: 1
    priority: 1
    shell:
        '''
        flow_diff_getheatmapdata \
        --rdspath '{input.rdspath}' \
        --savepath '{output}' \
        --top 50
        '''
        
rule diff_getdata_heatmap_top50_scale:
    input:
        rdspath=config["path"]["diff_filter_rds_savepath"]+"/diff_filter-{compare}.rds"
    output:
        "{savepath}/heatmap_top50_scale-{compare}.xls"
    threads: 1
    priority: 1
    shell:
        '''
        flow_diff_getheatmapdata \
        --rdspath '{input.rdspath}' \
        --savepath '{output}' \
        --top 50 \
        --scale
        '''

rule diff_plot_heatmap_top50:
    input:
        "{savepath}/heatmap_top50-{compare}.xls"
    output:
        "{savepath,((?!cluster).)*}/heatmap_top50-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        classfile = config["path"]["classtypefile"],
        scale = lambda wildcards: "'none' --log" if config["compare"][wildcards.compare]["samplenum"] == 2 else "row"
    threads: 1
    priority: 1
    shell:
        '''
        map_common_heatmap \
        --filename '{input}' \
        --savepath '{wildcards.savepath}' \
        --mapname 'heatmap_top50-{wildcards.compare}' \
        --imagetype '{wildcards.imagetype}' \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height} \
        --show_rownames \
        --show_colnames \
        --cluster_rows \
        --cluster_cols \
        --scale {params.scale} \
        --classfile {params.classfile} \
        --colgroup 1
        '''

rule diff_plot_heatmap_top50_cluster_none:
    input:
        "{savepath}/heatmap_top50-{compare}.xls"
    output:
        "{savepath}/cluster_none/heatmap_top50-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        classfile = config["path"]["classtypefile"],
        scale = lambda wildcards: "'none' --log" if config["compare"][wildcards.compare]["samplenum"] == 2 else "row"
    threads: 1
    priority: 1
    shell:
        '''
        map_common_heatmap \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/cluster_none/' \
        --mapname 'heatmap_top50-{wildcards.compare}' \
        --imagetype '{wildcards.imagetype}' \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height} \
        --show_rownames \
        --show_colnames \
        --cluster_rows \
        --scale {params.scale} \
        --classfile {params.classfile} \
        --colgroup 1
        '''

rule diff_plot_heatmap_top50:
    input:
        "{savepath}/heatmap_top50-{compare}.xls"
    output:
        "{savepath}/cluster_samples/heatmap_top50-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        classfile = config["path"]["classtypefile"],
        scale = lambda wildcards: "'none' --log" if config["compare"][wildcards.compare]["samplenum"] == 2 else "row"
    threads: 1
    priority: 1
    shell:
        '''
        map_common_heatmap \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/cluster_samples/' \
        --mapname 'heatmap_top50-{wildcards.compare}' \
        --imagetype '{wildcards.imagetype}' \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height} \
        --show_rownames \
        --show_colnames \
        --cluster_rows \
        --cluster_cols \
        --scale {params.scale} \
        --classfile {params.classfile} \
        --colgroup 1
        '''
