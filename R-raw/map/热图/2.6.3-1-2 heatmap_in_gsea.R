#!/opt/conda/bin/Rscript

#' @export
heatmap_in_gsea <- function(path = "./"){
  allheatmapfile <- list.files(path = path,
                               pattern = 'heatmap-data.xls$',full.names = T)
  if(length(allheatmapfile) == 0){
    savetxt(data = "富集结果过少，GSEA无结果",filename = paste0(path,"/说明.txt"))
    return()
  }
  
  allheatmapfile2 <- gsub(pattern = 'heatmap-data.xls$',replacement = 'heatmap',x = allheatmapfile)
  for ( i in 1:length(allheatmapfile)){
    
    map_common_heatmap(filename = allheatmapfile[i],
                       mapname = allheatmapfile2[i],
                       imagetype = c("png","pdf"), 
                       show_rownames = T,
                       show_colnames =T)
    
    # system(command = paste0('source /etc/profile;source ~/.bashrc;',packagepath(path = 'command'),'/map_common_heatmap -f "',
    #                         allheatmapfile[i],'" -mn "',allheatmapfile2[i],'" -i png pdf -sr -sc'))
  }
}
