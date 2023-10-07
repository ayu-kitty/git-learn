#!/opt/conda/bin/Rscript

#' 亚群比例差异分析
#' 
#' @param sceObject sce降维后的结果
#' @param clustering_to_use 用于亚群分析的细胞群，例如：meta12
#' @param batch 是否加入病人效应
#' @param cols_random_add 如果加入病人效应，则应该填入：patient_id
#' @param method_DA 差异分析方法
#' @param FDR_cutoff 精确度
#' @param saveplotpath 数据保存的路径
#' @param width 图片高度
#' @param height 图片宽度
#' @param dpi 图片精确度
#' @param family 图片字体
#' 
#' @export
sce_diffprops <- function(sceObject,
                          clustering_to_use = "meta20",
                          batch = F,
                          cols_random_add = NULL,
                          method_DA = NULL,
                          FDR_cutoff = NULL,
                          saveplotpath = "./",
                          color = SelectColors("corona"),
                          width = 8,height = 8,
                          ...){

  suppressMessages(library("CATALYST"))
  suppressMessages(library("cytofWorkflow"))
  suppressMessages(library("ggplot2"))
  suppressMessages(library("ComplexHeatmap"))
  suppressMessages(library("reshape2"))
  
  sceObject <- readdata(filename = sceObject)
  
  # 统计亚群的细胞数
  ns <- table(cluster_id = cluster_ids(sceObject, clustering_to_use), 
              sample_id = sample_ids(sceObject))
  # 计算亚群细胞数的比例
  fq <- proportions(ns, 2) * 100
  # 宽数据变长数据
  df <- as.data.frame(fq)
  # 长数据变宽数据
  dat <- dcast(df,cluster_id~sample_id)
  # sce对象的实验信息
  ei <- metadata(sceObject)$experiment_info
  # 对分组数据进行t.test检验
  dat$p <- apply(dat,1,function(x){t.test(as.numeric(x[-1])~ei$condition)$p.value})
  
  # method_DA
  if(is.null(method_DA)){
    method_DA <- "diffcyt-DA-GLMM"
  }
  # 精确度
  if(is.null(FDR_cutoff)){
    FDR_cutoff <- 0.05
  }
  contrast <- createContrast(c(0, 1))
  # 亚群的差异丰度 differential abundance (DA) of clusters
  da_formula <- createFormula(ei, 
                              cols_fixed = "condition", 
                              cols_random = "sample_id")
  # method_DA
  da_res <- diffcyt(sceObject, 
                    formula = da_formula, contrast = contrast,
                    analysis_type = "DA", method_DA = method_DA,
                    clustering_to_use = clustering_to_use, verbose = FALSE)
  # 去批次(加入病人效应)
  cols_random <- "sample_id"
  if(is.null(cols_random_add)){
    print("no capture batch infomation")
  }else{
    print("perform de-batching")
    cols_random <- c(cols_random, cols_random_add)
  }
  if(batch){
    da_formula <- createFormula(ei, 
                                cols_fixed = "condition", 
                                cols_random = cols_random)
    da_res <- diffcyt(sceObject, 
                      formula = da_formula, contrast = contrast,
                      analysis_type = "DA", method_DA = method_DA,
                      clustering_to_use = clustering_to_use, verbose = FALSE)
  }
  
  plotfile(savepath = saveplotpath,
           mapname = ifelse(clustering_to_use == "meta20","sce_subpopulation_ExprHeatmap","sce_celltype_ExprHeatmap"),
           width =  width,height = height,...)
  
  picture5 <- plotDiffHeatmap(sceObject, rowData(da_res$res), all = TRUE, fdr = FDR_cutoff, col_anno = c("group"))
  picture5@top_annotation@anno_list[["group"]]@fun@var_env[["color_mapping"]]@full_col[1:length(table(sceObject$cluster))] <- color[1:length(table(sceObject$cluster))]
  picture5@top_annotation@anno_list[["group"]]@color_mapping@full_col[1:length(table(sceObject$cluster))] <- color[1:length(table(sceObject$cluster))]
  draw(picture5)

  plotsave()
  
  # 保存亚群细胞占比信息
  df = topTable(da_res, show_props = TRUE, format_vals = TRUE, digits = 2,show_all_cols=T)
  df = as.data.frame(df)
  
  savexlsx1(data = df,
            filename = paste0(saveplotpath,"/",ifelse(clustering_to_use == "meta20","sce_subpopulation_proportions.xlsx","sce_celltype_proportions.xlsx")),
            sheet = "data")

  print("plot finished")
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-sce","--sceObject",default = "rawdata/flowSet_object_sce_reducedDim_celltype.rds", help = "读取sce对象的rds文件名")
  parser$add_argument("-clust","--clustering_to_use",default = 'meta20', help = " 用于亚群分析的细胞群，例如：meta12，默认meta20")
  parser$add_argument("-b","--batch",default = F, action = "store_true", help = "是否加入病人（批次）效应，默认不加")
  parser$add_argument("-ba","--cols_random_add",default = NULL, help = "如果加入病人效应，则应该填入：patient_id")
  parser$add_argument("-da","--method_DA",default = NULL, help = "差异分析方法，默认diffcyt-DA-GLMM")
  parser$add_argument("-fdr","--FDR_cutoff",default = NULL, help = "精确度，默认0.05")
  
  parser$add_argument("-sp","--saveplotpath",default = "6.Diffanalysis/", help = "细胞比例差异分析结果保存路径，默认在5.diffanalysis目录下")
  parser$add_argument("-wi","--width", default = 8, type = "integer", help = "图片宽度")
  parser$add_argument("-he","--height", default = 8, type = "integer", help = "图片长度")
  parser$add_argument("-i","--dpi", default = 300, type = "integer", help = "图片分辨率")
  parser$add_argument("-fa","--family", default = "sans", help = "字体")
  
  args <- parser$parse_args()
  
  try({
    diffargs <- do.call(what = sce_diffprops, args = args)
  },silent = F)
  
}

