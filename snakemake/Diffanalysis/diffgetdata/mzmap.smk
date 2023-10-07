rule get_mzmap:
    input:
        infofile="内部分析单.xlsx"
    output:
        directory("{savepath,.*质谱图}")
    threads: 1
    priority: 8
    shell:
        '''
        /opt/conda/bin/Rscript -e "lmbio::mzmltoimage_obj( \
        filename = '{input.infofile}', \
        savepath = '{output}')"
        '''

rule get_mzmap:
    input:
        infofile="内部分析单.xlsx"
    output:
        directory("{savepath,.*色谱图}")
    threads: 1
    priority: 8
    shell:
        '''
        /opt/conda/bin/Rscript -e "lmbio::mzmltoimage_obj( \
        filename = '{input.infofile}', \
        savepath = '{output}')"
        '''
