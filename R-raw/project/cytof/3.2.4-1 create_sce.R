#!/opt/conda/bin/Rscript

#' 构建SingleCellExperiment对象
#' 
#' @param flowsetObject flowSet对象的变量名
#' @param easy_sce 逻辑，是否需要简单构建sce对象，默认不做
#' @param metadata metadata的文件名
#' @param panel  panel的文件名
#' @param savesce 逻辑，是否保存SingleCellExperiment对象的rds文件
#' @param savepath 保存flowset_Object_se.rds文件的路径
#' @param scename 保存的sce文件名
#' 
#' @export
create_sce <- function(flowsetObject = "rawdata/flowSet_object.rds",
                       easy_sce = F,
                       metadata = "rawdata/metadata.xlsx",
                       panel = "rawdata/panel.xlsx",
                       savepath = "rawdata",
                       scename = "flowSet_object_sce.rds",
                       by_time = FALSE,
                       ...){
  
  suppressMessages(library("stringr"))
  suppressMessages(library("CATALYST"))
  
  flowsetObject <- readdata(filename = flowsetObject)
  
  # 简单构建SingleCellExperiment对象
  if(easy_sce){
    sce <- prepData(flowsetObject, by_time = by_time)
    print(paste0(date(), " ", "easy SingleCellExperiment object create success"))
  }
  
  if(!easy_sce){
    # flowsetObject、md、panel共同构建SingleCellExperiment对象
    panel <- readdata(panel)
    metadata <- readdata(metadata)

    # spot check that all panel columns are in the flowSet object
    all(panel$fcs_colname %in% colnames(flowsetObject))
    # specify levels for groups & sample IDs to assure desired ordering
    metadata$condition <- factor(metadata$condition, levels = names(table(metadata$condition)))
    metadata$sample_id <- factor(metadata$sample_id, 
                                 levels = metadata$sample_id[order(metadata$condition)])
    sceObject <- prepData(flowsetObject, panel, metadata, by_time = by_time) # features = NULL
    print(paste0(date(), " ", "whole SingleCellExperiment object create success"))
  }
  sceObject@colData$sample <- factor(sample_ids(sceObject))
  sceObject@colData$group <- factor(sceObject@colData$condition)
  # 保存RDS

  # print(paste0(date(), " " , "prepare to save RDS..."))
  saverds(data = sceObject,filename = paste0(savepath,"/",scename))
  # print("save RDS success")
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-fo","--flowsetObject",default = "rawdata/flowSet_object.rds", help = "读取fcs的结果文件：rawdata/flowSet_object.rds")
  parser$add_argument("-sp","--savepath",default = "rawdata", help = "保存SingleCellExperiment对象的rds文件的路径，默认当前路径（项目路径）")
  parser$add_argument("-meta","--metadata", default = "rawdata/metadata.xlsx", help = "metadata的文件名")
  parser$add_argument("-pl","--panel", default = "rawdata/panel.xlsx", help = "panel的文件名")
  parser$add_argument("-sn","--scename", default = "flowSet_object_sce.rds", help = " 保存的sce文件名")
  
  # 特殊参数
  parser$add_argument("-easy","--easy_sce", default = F, action = "store_true", help = "是否需要简单构建sce对象，默认不做")
  
  args <- parser$parse_args()
  
  sceargs <- do.call(what = create_sce, args = args)
}
## 注意：每个目录和文件必须严格存在！rawdata是固定的下机数据存放目录 flowSet_object.rds是固定的文件名
