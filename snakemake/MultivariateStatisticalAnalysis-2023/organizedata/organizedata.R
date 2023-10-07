#!/opt/conda/bin/Rscript

suppressMessages(library("lmbio"))

# 参数
rawfile=snakemake@input$rawfile
datafile=snakemake@output$datafile
infofile=snakemake@output$infofile
classfile=snakemake@output$classfile
classtypefile=snakemake@output$classtypefile
# rawfile="/luming/test/20220817-mulstatisticsflow/test1/数据矩阵.xlsx"
# datafile="oecloud/rawdata/datafile.txt"
# infofile="oecloud/rawdata/infofile.txt"
# classfile="oecloud/rawdata/classfile.yaml"

# logfile <- file("log.txt", open = "a")
# sink(logfile,append = T)
# sink(logfile, type = "message",append = T)

organizedata(rawfile = rawfile,
             datafile = datafile,
             infofile = infofile,
             classfile = classfile,
             classtypefile = classtypefile)
