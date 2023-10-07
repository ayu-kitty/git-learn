#!/opt/conda/bin/Rscript

suppressMessages(library("lmbio"))

# 参数
# rdspath="oecloud/mulstatisticsanalyst/mulstatistic_opls-KO-Air-vs-WT-Air.rds"
# datasavepath="./"
# savename="KO-Air-vs-WT-Air"
# infofile="oecloud/rawdata/infofile.txt"

rdspath=snakemake@input$rdspath
savepath=snakemake@output[[1]]
infofile=snakemake@config$path$infofile

# logfile <- file("log.txt", open = "a")
# sink(logfile,append = T)
# sink(logfile, type = "message",append = T)

result <- getloading_file(rdspath = rdspath,
                          infofile = infofile,
                          savepath = savepath)
