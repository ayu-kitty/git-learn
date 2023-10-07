rule Ksea:
    input:
        expand(config["savepath"]+"{group}/",group=config['compare'])

rule Ksea_phos:
    input:
        datafile=config["path"]+config["ptmdata"]
    output:
        directory(config["savepath"]+"{group}/")
    threads: 1
    priority: 0
    params:
        org=config["org"],
        inputpath=config["path"],
        savepath=config["savepath"]
    shell:
        "KSEA_Activity -ic '{wildcards.group}' -s '{params.savepath}' -ip '{params.inputpath}' -or '{params.org}' -p '{input.datafile}'"


