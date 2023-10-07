#!/opt/conda/bin/Rscript

#' SingleCellExperiment对象操作之：降维
#' 
#' @param sceObject sce对象的变量名
#' @param cells_number 运行t-SNE/UMAP时选取的细胞数，默认1e4个
#' @param features 分群的抗体，features = "type" 表示谱系标签/ "state" 表示功能标签/ NULL 无指定,默认NULL
#' @param savesce 逻辑，是否保存降维后的sce对象
#' @param savename 保存为sce对象的名字
#' @param ...
#' 
#' @export
sce_reducedim <- function(sceObject,
                          cells_number = 1e4,
                          features = NULL,
                          savepath = "rawdata",
                          savename = "flowSet_object_sce_cluster_reducedDim.rds",
                          seed = 1,
                          ...){
  suppressMessages(library("CATALYST"))
  
  sceObject <- readdata(filename = sceObject)
  
  # run t-SNE/UMAP on at most 500/1000 cells per sample
  set.seed(seed = seed)
  sceObject <- runDR(sceObject, "TSNE", cells = cells_number, features = features)
  sceObject <- runDR(sceObject, "UMAP", cells = cells_number, features = features)

  # print(paste0(date(), " ", "prepare to save RDS..."))
  saverds(data = sceObject,filename = paste0(savepath,"/",savename))
  # print("save RDS success")

}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-sce","--sceObject",default = "rawdata/flowSet_object_sce_cluster.rds", help = "sce对象的聚类结果rds名：flowSet_object_sce_cluster.rds")
  parser$add_argument("-fe","--features",default = NULL, help = " 抗体类型，features = 'type' 表示谱系标签/ 'state' 表示状态标签/ NULL 无指定,默认NULL")
  parser$add_argument("-cell_num","--cells_number",default = 10000, type = "integer", help = "运行t-SNE/UMAP时选取的细胞数，默认1e4个")
  parser$add_argument("-sp","--savepath", default = "rawdata", help = "保存聚类结果文件的路径")
  parser$add_argument("-sn","--savename", default = "flowSet_object_sce_cluster_reducedDim.rds", help = "保存为sce对象的名字")
  parser$add_argument("-se","--seed", default = 1, type = "integer", help = "随机种子数")
  
  args <- parser$parse_args()
  
  dimargs <- do.call(what = sce_reducedim, args = args)
}

