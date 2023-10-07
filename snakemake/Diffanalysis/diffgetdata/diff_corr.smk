
rule diff_getdata_cor:
    input:
        rdspath=config["path"]["diff_filter_rds_savepath"]+"/diff_filter-{compare}.rds"
    output:
        "{savepath}/Correlation_expression-{compare}.xls",
    threads: 1
    shell:
        '''
        flow_diff_getheatmapdata \
        --rdspath '{input.rdspath}' \
        --savepath '{output}' \
        --addgroup 'F'
        '''

rule diff_plot_cor:
    input:
        "{savepath}/Correlation_expression-{compare}.xls"
    output:
        expand("{{savepath}}/Correlation/Correlation-{{compare}}.{imagetype}",
               imagetype=config["params"]["imagetype"])
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_common_corrplot \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/Correlation' \
        --mapname 'Correlation-{wildcards.compare}' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''

rule diff_plot_corrnetwork:
    input:
        "{savepath}/Correlation_expression-{compare}.xls"
    output:
        expand("{{savepath}}/Corrnetwork/Corrnetwork-{{compare}}.{imagetype}",
               imagetype=config["params"]["imagetype"])
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_common_corrnetwork \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/Corrnetwork' \
        --mapname 'Corrnetwork-{wildcards.compare}' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''

  
rule diff_getdata_cor_top50:
    input:
        rdspath=config["path"]["diff_filter_rds_savepath"]+"/diff_filter-{compare}.rds"
    output:
        "{savepath}/Correlation_expression_top50-{compare}.xls"
    threads: 1
    shell:
        '''
        flow_diff_getheatmapdata \
        --rdspath '{input.rdspath}' \
        --savepath '{output}' \
        --addgroup 'F' \
        --top 50
        '''

rule diff_plot_cor_top50:
    input:
        "{savepath}/Correlation_expression_top50-{compare}.xls"
    output:
        expand("{{savepath}}/Correlation_top50/Correlation_top50-{{compare}}.{imagetype}",
               imagetype=config["params"]["imagetype"])
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_common_corrplot \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/Correlation_top50' \
        --mapname 'Correlation_top50-{wildcards.compare}' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''
        
rule diff_plot_corrnetwork_top50:
    input:
        "{savepath}/Correlation_expression_top50-{compare}.xls"
    output:
        expand("{{savepath}}/Corrnetwork_top50/Corrnetwork_top50-{{compare}}.{imagetype}",
               imagetype=config["params"]["imagetype"])
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_common_corrnetwork \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/Corrnetwork_top50' \
        --mapname 'Corrnetwork_top50-{wildcards.compare}' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''
        
rule diff_getdata_cor_top20:
    input:
        rdspath=config["path"]["diff_filter_rds_savepath"]+"/diff_filter-{compare}.rds"
    output:
        "{savepath}/Correlation_expression_top20-{compare}.xls"
    threads: 1
    shell:
        '''
        flow_diff_getheatmapdata \
        --rdspath '{input.rdspath}' \
        --savepath '{output}' \
        --addgroup 'F' \
        --top 20
        '''

rule diff_plot_cor_top20:
    input:
        "{savepath}/Correlation_expression_top20-{compare}.xls"
    output:
        expand("{{savepath}}/Correlation_top20/Correlation_top20-{{compare}}.{imagetype}",
               imagetype=config["params"]["imagetype"])
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_common_corrplot \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/Correlation_top20' \
        --mapname 'Correlation_top20-{wildcards.compare}' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''
        
rule diff_plot_corrnetwork_top20:
    input:
        "{savepath}/Correlation_expression_top20-{compare}.xls"
    output:
        expand("{{savepath}}/Corrnetwork_top20/Corrnetwork_top20-{{compare}}.{imagetype}",
               imagetype=config["params"]["imagetype"])
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_common_corrnetwork \
        --filename '{input}' \
        --savepath '{wildcards.savepath}/Corrnetwork_top20' \
        --mapname 'Corrnetwork_top20-{wildcards.compare}' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''
