#!/opt/conda/bin/Rscript

#' SPADE分析
#'
#' @param flowSetObject 读取fcs的结果文件：**.rds
#' @param data_file_path fcs文件所在的目录，默认rawdata
#' @param Desired_number_of_clusters SPADE目标聚类群数，默认70
#'
#' @export
plot_spade <- function(flowsetObject = "./rawdata/flowSet_object.rds",
                       data_file_path = "rawdata",
                       Desired_number_of_clusters = 70,
                       savemiddlepath = "./rawdata/spade",
                       savepath = "2.SPADE",
                       downsampling_target_number = NULL, 
                       downsampling_target_pctile = NULL, 
                       downsampling_target_percent = 0.1,
                       ...){
  print("开始SPADE聚类分析")
  suppressMessages(library("spade"))
  suppressMessages(library("stringr"))
  
  # get sample data marker used to cluster
  fs_object <- readRDS(flowsetObject)
  marker <- colnames(fs_object)
  marker_input <- marker[-which(tolower(marker) %in% tolower(c("Event_length", "Center", "Offset","Residual", "Time", "Width", "191", "193")))]
  
  # Run basic SPADE analyses, clustering on two parameters. 
  for(fcs in rownames(fs_object@phenoData)){
    output_dir <- paste0(savemiddlepath, "/", unlist(str_split(fcs,".fcs"))[1])
    createdir(output_dir,force = T)
    spadetry <- try({
      SPADE.driver(file.path(data_file_path,fcs), 
                   out_dir = output_dir, 
                   cluster_cols = marker_input, 
                   k = Desired_number_of_clusters,
                   downsampling_target_number = downsampling_target_number, 
                   downsampling_target_pctile = downsampling_target_pctile, 
                   downsampling_target_percent = downsampling_target_percent,
                   ...)
    },silent = F)
    
    if("try-error" %in% class(spadetry)){
      createdir(output_dir,force = T)
      SPADE.driver(file.path(data_file_path,fcs), 
                   out_dir = output_dir, 
                   cluster_cols = marker_input, 
                   k = Desired_number_of_clusters,
                   downsampling_target_number = NULL, 
                   downsampling_target_pctile = NULL, 
                   downsampling_target_percent = 1,
                   ...)
    }

    mst_graph <- igraph:::read.graph(paste(output_dir,"mst.gml",sep=.Platform$file.sep),format="gml")
    # Generate PDFs of annotated graphs (into output_dir/pdf)
    plor_dir <- paste0(savepath,"/",unlist(str_split(fcs,".fcs"))[1])
    createdir(plor_dir)
    SPADE.plot.trees(graph = mst_graph, 
                     files = output_dir, 
                     out_dir = plor_dir, 
                     layout=igraph:::layout.kamada.kawai(mst_graph))
  }
  
  system(paste0("for i in `ls ",savepath,"`;do mogrify -path ",savepath,"/$i/ -format jpg  ",savepath,"/$i/*.pdf;done"))
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-fo","--flowsetObject",default = "./rawdata/flowSet_object.rds", help = "读取fcs的结果文件：**.rds")
  parser$add_argument("-dp","--data_file_path",default = "./rawdata", help = "fcs文件所在的目录，默认rawdata")
  parser$add_argument("-num","--Desired_number_of_clusters", default = 70, type = "integer", help = "SPADE目标聚类群数，默认70")
  parser$add_argument("-sp","--savepath",default = "2.SPADE", help = "SPADE保存路径")
  
  args <- parser$parse_args()
  
  spadeargs <- do.call(what = plot_spade, args = args)
}
## 注意：每个目录和文件必须严格存在！rawdata是固定的下机数据存放目录 flowSet_object.rds是固定的文件名