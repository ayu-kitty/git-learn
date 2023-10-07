#!/opt/conda/bin/Rscript

suppressMessages(library("lmbio"))

# 参数
# classfile="oecloud/rawdata/classfile.yaml"
# datafile="oecloud/rawdata/datafile.txt"
# mulstatistics_rds_savepath="oecloud/mulstatisticsanalyst/mulstatistic_pca-WT-vs-KO"
# group=c("WT","KO")
# name="WT/KO"
# log10L=F
# orthoI=0
# permI=0
# predI=-1
# scaleC="standard"
# mode="PCA"
# amount=3
datafile=snakemake@input$datafile
classfile=snakemake@input$classfile
mulstatistics_rds_savepath=snakemake@output$rdssavepath
group=snakemake@params$group
name=snakemake@params$name
log10L=snakemake@params$log10L
orthoI=snakemake@params$orthoI
permI=snakemake@params$permI
predI=snakemake@params$predI
scaleC=snakemake@params$scaleC
mode=snakemake@params$mode
amount=snakemake@params$amount
both=snakemake@params$both

# logfile <- file("log.txt", open = "a")
# sink(logfile,append = T)
# sink(logfile, type = "message",append = T)

result <- mulstatistics_file(datafile = datafile,
                             classfile = classfile,
                             group = group,
                             name = gsub(pattern = "-vs-",replacement = "/",x = name),
                             mulstatistics_rds_savepath = mulstatistics_rds_savepath,
                             log10L = log10L,
                             orthoI = orthoI,
                             permI = permI,
                             predI = predI,
                             scaleC = scaleC,
                             mode = mode,
                             both = both,
                             amount = amount) 

