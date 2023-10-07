#!/opt/conda/bin/Rscript

suppressMessages(library("lmbio"))

# 参数
# rdspath=c("oecloud/mulstatisticsanalyst/mulstatistic_opls-KO-Air-vs-WT-Air.rds",
#           "oecloud/mulstatisticsanalyst/mulstatistic_pls-KO-Air-vs-WT-Air.rds",
#           "oecloud/mulstatisticsanalyst/mulstatistic_pca-KO-Air-vs-WT-Air.rds")
# datasavepath="./"
# savename="KO-Air-vs-WT-Air"

rdspath=Reduce(c,snakemake@input)
savepath=snakemake@output[[1]]

# logfile <- file("log.txt", open = "a")
# sink(logfile,append = T)
# sink(logfile, type = "message",append = T)

result <- getsummary_file(rdspath = rdspath,
                          savepath = savepath)
