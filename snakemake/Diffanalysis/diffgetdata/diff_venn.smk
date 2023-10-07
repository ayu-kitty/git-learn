
rule diff_plot_venn_top5:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{compare}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"],
               compare=list(config["params"]["all_compare"].keys())[:5] if len(config["params"]["all_compare"]) > 5 else config["params"]["all_compare"])
    output:
        expand("{{savepath}}/Venn.{imagetype}",
               imagetype=config["params"]["imagetype"])
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        filename = lambda wildcards, output, input: "'"+"' '".join(input)+"'",
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_common_venn \
        --filename {params.filename} \
        --savepath '{wildcards.savepath}' \
        --mapname 'Venn' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''

rule diff_plot_upset:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{compare}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"],
               compare=config["params"]["all_compare"])
    output:
        expand("{{savepath}}/Upset.{imagetype}",
               imagetype=config["params"]["imagetype"])
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        filename = lambda wildcards, output, input: "'"+"' '".join(input)+"'",
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_common_venn \
        --filename {params.filename} \
        --savepath '{wildcards.savepath}' \
        --mapname 'Upset' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height} \
        --upset
        '''

rule diff_plot_flower:
    input:
        expand("{diff_filter_rds_savepath}/diff_filter-{compare}.rds",
               diff_filter_rds_savepath=config["path"]["diff_filter_rds_savepath"],
               compare=config["params"]["all_compare"])
    output:
        expand("{{savepath}}/Flower.{imagetype}",
               imagetype=config["params"]["imagetype"])
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        filename = lambda wildcards, output, input: "'"+"' '".join(input)+"'",
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_common_venn \
        --filename {params.filename} \
        --savepath '{wildcards.savepath}' \
        --mapname 'Flower' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''
