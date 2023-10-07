#!/opt/conda/bin/Rscript

#' SingleCellExperiment对象操作之：分群
#' 
#' @param sceObject sce对象的变量名
#' @param features 抗体类型，features = "type" 表示谱系标签/ "state" 表示状态标签/ NULL 无指定,默认NULL
#' @param xdim 分群的X维度，默认10
#' @param ydim 分群的Y维度，默认10
#' @param maxK 分群的最大群数，默认不超过20个群
#' @param savesce 逻辑，是否保存分群后的sce对象
#' @param savepath 保存聚类结果文件的路径
#' @param savename 保存为sce对象的名字
#' @param ...
#' 
#' @export
sce_cluster <- function(sceObject,
                        features = NULL,
                        xdim = 10,
                        ydim = 10,
                        maxK = 20,
                        savepath = "rawdata",
                        savename = "flowSet_object_sce_cluster.rds",
                        seed = 1,
                        ...){
  suppressMessages(library("CATALYST"))
  
  sceObject <- readdata(filename = sceObject)
  sceObject <- cluster(sceObject, features = features,
                       xdim = xdim, ydim = ydim, maxK = maxK,seed = seed)
  
  sceObject@colData$cluster <- factor(cluster_ids(sceObject, paste0("meta", maxK)))
  # print(paste0(date(), " ", "prepare to save RDS..."))
  saverds(data = sceObject,filename = paste0(savepath,"/",savename))
  # print("save RDS success")
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-sce","--sceObject",default = "rawdata/flowSet_object_sce.rds", help = "sce对象的rds名：flowSet_object_sce.rds")
  parser$add_argument("-fe","--features",default = NULL, help = " 抗体类型，features = 'type' 表示谱系标签/ 'state' 表示状态标签/ NULL 无指定,默认NULL")
  parser$add_argument("-x","--xdim", default = 10, type = "integer", help = "分群的X维度，默认10")
  parser$add_argument("-y","--ydim", default = 10, type = "integer", help = "分群的Y维度，默认10")
  parser$add_argument("-max","--maxK", default = 20, type = "integer", help = "分群的最大群数，默认不超过20个群")
  parser$add_argument("-se","--seed", default = 1, type = "integer", help = "随机种子数")
  parser$add_argument("-sp","--savepath", default = "rawdata", help = "保存聚类结果文件的路径")
  parser$add_argument("-sn","--savename", default = "flowSet_object_sce_cluster.rds", help = "保存为sce对象的名字")
  
  args <- parser$parse_args()
  
  clusterargs <- do.call(what = sce_cluster, args = args)
}
