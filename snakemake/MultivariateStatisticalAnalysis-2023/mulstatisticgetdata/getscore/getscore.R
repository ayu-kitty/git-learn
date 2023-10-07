#!/opt/conda/bin/Rscript

suppressMessages(library("lmbio"))

rdspath=snakemake@input$rdspath
savepath=snakemake@output[[1]]
savename=snakemake@wildcards$compare

# logfile <- file("log.txt", open = "a")
# sink(logfile,append = T)
# sink(logfile, type = "message",append = T)

result <- getscore_file(rdspath = rdspath,
                        savepath = savepath)
