#!/opt/conda/bin/Rscript

#' SingleCellExperiment对象操作之：降维绘图
#' 
#' @param sceObject sce对象的变量名
#' @param saveplotpath 保存图片的目录
#' @param color_by 绘图的颜色以该参数区分，可选：sample_id、condition、patient_id
#' @param facet_by 设置降维的分面，可选：sample_id、condition、patient_id
#' @param hm1 同features，抗体类型type、state或NULL
#' @param hm2 感兴趣的抗体，例如CD4
#' @param PCA_k 绘制代表由 FlowSOM 推断的 tSNE 和 PCA时，指定集群，默认meta20
#' @param som_k 由 FlowSOM 选取的xdim*ydim，默认som100
#' @param som_m 由 FlowSOM 推断的亚群数，默认meta20
#' @param width 图片宽度
#' @param height 图片高度
#' @param dpi 图片清晰度
#' @param ...
#' 
#' @export
sce_reducedim_plot <- function(sceObject,
                               saveplotpath = "4.Reducedim/",
                               color_by = "cluster",
                               facet_by = "sample",
                               hm1 = NULL,
                               hm2 = "abundances",
                               PCA_k = "meta20",
                               som_k = "som100",
                               som_m = "meta20",
                               samplecolor = SelectColors("customecol2"),
                               ...){

  suppressMessages(library("CATALYST"))
  suppressMessages(library("cytofWorkflow"))
  suppressMessages(library("ggplot2"))
  suppressMessages(library("HDCytoData"))
  
  sceObject <- readdata(filename = sceObject)
  
  # 绘制降维后的图片
  #plotDR(sce, "TSNE", color_by = color_by) # 可以根据metadata的列名绘图、抗体名字、"meta20"，默认"meta20"
  # t-SNE、UMAP
  picture1 <- plotDR(sceObject, "TSNE", color_by = color_by) + 
    theme(legend.position = "right") + 
    scale_color_manual(values = samplecolor) + labs(color = color_by) + 
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  picture2 <- plotDR(sceObject, "UMAP", color_by = color_by) + 
    theme(legend.position = "right") + 
    scale_color_manual(values = samplecolor) + labs(color = color_by) + 
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  ggplotsave(plot = picture1, savepath = saveplotpath, mapname = "TSNE",
             width = 8, height = 8,...)
  ggplotsave(plot = picture2, savepath = saveplotpath, mapname = "UMAP",
             width = 8, height = 8,...)
  
  # facet by sample
  picture3 <- plotDR(sceObject, "TSNE", color_by = color_by, facet_by = facet_by) +
    scale_color_manual(values = samplecolor) + labs(color = color_by) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  picture4 <- plotDR(sceObject, "UMAP", color_by = color_by, facet_by = facet_by) +
    scale_color_manual(values = samplecolor) + labs(color = color_by) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  ggplotsave(plot = picture3, savepath = saveplotpath, mapname = "sce_reducedDim_TSNE_facet",
             width = 8, height = 8,...)
  ggplotsave(plot = picture4, savepath = saveplotpath, mapname = "sce_reducedDim_UMAP_facet",
             width = 8, height = 8,...)
  
  # TSNE for every antigen
  for(i in 1:length(names(sceObject))){
    temp_plot <- plotDR(sceObject, "TSNE", color_by = names(sceObject)[i])
    ggplotsave(plot = temp_plot, mapname = paste0(names(sceObject)[i], "_", "TSNE"), 
               savepath = paste0(saveplotpath, "/", "tsne4antigen"),
               width = 8, height = 8,...)
  }
  
  # UMAP for every antigen
  for(i in 1:length(names(sceObject))){
    temp_plot <- plotDR(sceObject, "UMAP", color_by = names(sceObject)[i])
    ggplotsave(plot = temp_plot, mapname = paste0(names(sceObject)[i], "_", "UMAP"), 
               savepath = paste0(saveplotpath, "/", "umap4antigen"),
               width = 8, height = 8,...)
  }
  
  savexlsx1(data = picture1$data,
            filename = paste0(saveplotpath, "/reduced_dim_info.xlsx"),
            sheet = "TSNE")
  savexlsx1(data = picture2$data,
            filename = paste0(saveplotpath, "/reduced_dim_info.xlsx"),
            sheet = "UMAP")
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-sce","--sceObject",default = "rawdata/flowSet_object_sce_cluster_reducedDim.rds", help = "读取sce对象的rds文件名")
  parser$add_argument("-color","--color_by",default = "cluster", help = "绘图的颜色以该参数区分，可选：sample_id、condition、patient_id")
  parser$add_argument("-facet","--facet_by",default = "sample", help = "设置降维的分面，可选：sample_id、condition、patient_id")
  parser$add_argument("-h1","--hm1",default = NULL, help = "同features，抗体类型type、state或NULL")
  parser$add_argument("-h2","--hm2",default = "abundances", help = "感兴趣的抗体，例如CD4，默认abundances")
  parser$add_argument("-pca","--PCA_k",default = "meta20", help = "绘制代表由 FlowSOM 推断的 tSNE 和 PCA时，指定集群，默认meta20")
  parser$add_argument("-k","--som_k",default = "som100", help = "som_k 由 FlowSOM 选取的xdim*ydim，默认som100")
  parser$add_argument("-m","--som_m",default = "meta20", help = "som_m 由 FlowSOM 推断的亚群数，默认meta20")
  
  parser$add_argument("-sp","--saveplotpath",default = "4.Reducedim/", help = "保存图片的目录")
  parser$add_argument("-i","--dpi", default = 300, type = "integer", help = "图片清晰度")
  parser$add_argument("-fa","--family", default = "sans", help = "图片字体")
  
  args <- parser$parse_args()
  
  dimplotargs <- do.call(what = sce_reducedim_plot, args = args)
}
