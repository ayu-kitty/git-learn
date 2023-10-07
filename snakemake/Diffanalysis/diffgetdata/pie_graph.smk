rule get_plot_pie:
    input:
       "raw/数据矩阵.xlsx"
    output:
       directory("{savepath,.*Pie_Graph}")
    params:
        family = config["params"]["plot"]["family"],
        width = config["params"]["plot"]["width"],
        height = config["params"]["plot"]["height"],
        classfile = config["path"]["classtypefile"],
        imagetype = config["params"]["imagetype"]
    threads: 1
    shell:
        '''
        map_common_pie \
        --filename '{input}' \
        --sheet 数据矩阵 \
        --savepath '{output}' \
        --mapname 'pie' \
        --imagetype {params.imagetype} \
        --family '{params.family}' \
        --width {params.width} \
        --height {params.height}
        '''
