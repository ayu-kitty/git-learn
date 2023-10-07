
if "volcanonum" not in config["params"]["plot"]:
  config["params"]["plot"]["volcanonum"]=5

rule diff_getdata_volcano:
    input:
        rdspath=config["path"]["diff_filter_rds_savepath"]+"/diff_filter-{compare}.rds"
    output:
        "{savepath}/volcano-{compare}.xls"
    threads: 1
    priority: 1
    shell:
        '''
        flow_diff_getvolcanodata \
        --rdspath '{input.rdspath}' \
        --savepath '{output}'
        '''

if config["params"]["omic"] == "PF" :
  if config["params"]["diff_filter"]["adjpfilter"] > 0:
    rule diff_plot_volcano:
        input:
            "{savepath}/volcano-{compare}.xls"
        output:
            "{savepath}/volcano-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
        params:
            family = config["params"]["plot"]["family"],
            width = config["params"]["plot"]["width"],
            height = config["params"]["plot"]["height"],
            fcfilter = config["params"]["diff_filter"]["fcfilter"],
            vipfilter = config["params"]["diff_filter"]["vipfilter"],
            pfilter = config["params"]["diff_filter"]["adjpfilter"],
            volcanonum = config["params"]["plot"]["volcanonum"]
        threads: 1
        priority: 1
        shell:
            '''
            map_common_volcano \
            --filename '{input}' \
            --savepath '{wildcards.savepath}' \
            --mapname 'volcano-{wildcards.compare}' \
            --imagetype '{wildcards.imagetype}' \
            --fcfilter '{params.fcfilter}' \
            --vipfilter {params.vipfilter} \
            --pfilter {params.pfilter} \
            --family '{params.family}' \
            --width {params.width} \
            --height {params.height} \
            --color "blue-grey-grey-grey-red" \
            --showname T \
            --namenum {params.volcanonum} \
            --labelname "projectdata/靶蛋白.xlsx" \
            --adjust
            '''
  else:
    rule diff_plot_volcano:
        input:
            "{savepath}/volcano-{compare}.xls"
        output:
            "{savepath}/volcano-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
        params:
            family = config["params"]["plot"]["family"],
            width = config["params"]["plot"]["width"],
            height = config["params"]["plot"]["height"],
            fcfilter = config["params"]["diff_filter"]["fcfilter"],
            vipfilter = config["params"]["diff_filter"]["vipfilter"],
            pfilter = config["params"]["diff_filter"]["pfilter"],
            volcanonum = config["params"]["plot"]["volcanonum"]
        threads: 1
        priority: 1
        shell:
            '''
            map_common_volcano \
            --filename '{input}' \
            --savepath '{wildcards.savepath}' \
            --mapname 'volcano-{wildcards.compare}' \
            --imagetype '{wildcards.imagetype}' \
            --fcfilter '{params.fcfilter}' \
            --vipfilter {params.vipfilter} \
            --pfilter {params.pfilter} \
            --family '{params.family}' \
            --width {params.width} \
            --height {params.height} \
            --color "blue-grey-grey-grey-red" \
            --showname T \
            --namenum {params.volcanonum} \
            --labelname "projectdata/靶蛋白.xlsx"
            '''
else:
  if config["params"]["diff_filter"]["adjpfilter"] > 0:
    rule diff_plot_volcano:
        input:
            "{savepath}/volcano-{compare}.xls"
        output:
            "{savepath}/volcano-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
        params:
            family = config["params"]["plot"]["family"],
            width = config["params"]["plot"]["width"],
            height = config["params"]["plot"]["height"],
            fcfilter = config["params"]["diff_filter"]["fcfilter"],
            vipfilter = config["params"]["diff_filter"]["vipfilter"],
            pfilter = config["params"]["diff_filter"]["adjpfilter"]
        threads: 1
        priority: 1
        shell:
            '''
            map_common_volcano \
            --filename '{input}' \
            --savepath '{wildcards.savepath}' \
            --mapname 'volcano-{wildcards.compare}' \
            --imagetype '{wildcards.imagetype}' \
            --fcfilter '{params.fcfilter}' \
            --vipfilter {params.vipfilter} \
            --pfilter {params.pfilter} \
            --family '{params.family}' \
            --width {params.width} \
            --height {params.height} \
            --adjust
            '''
  else:
    rule diff_plot_volcano:
        input:
            "{savepath}/volcano-{compare}.xls"
        output:
            "{savepath}/volcano-{compare}.{imagetype,(pdf)|(jpg)|(png)|(tiff)|(html)}"
        params:
            family = config["params"]["plot"]["family"],
            width = config["params"]["plot"]["width"],
            height = config["params"]["plot"]["height"],
            fcfilter = config["params"]["diff_filter"]["fcfilter"],
            vipfilter = config["params"]["diff_filter"]["vipfilter"],
            pfilter = config["params"]["diff_filter"]["pfilter"],
            adjpfilter = config["params"]["diff_filter"]["adjpfilter"]
        threads: 1
        priority: 1
        shell:
            '''
            map_common_volcano \
            --filename '{input}' \
            --savepath '{wildcards.savepath}' \
            --mapname 'volcano-{wildcards.compare}' \
            --imagetype '{wildcards.imagetype}' \
            --fcfilter '{params.fcfilter}' \
            --vipfilter {params.vipfilter} \
            --pfilter {params.pfilter} \
            --family '{params.family}' \
            --width {params.width} \
            --height {params.height}
            '''
