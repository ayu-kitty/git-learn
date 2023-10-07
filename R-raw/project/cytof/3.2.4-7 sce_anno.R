#!/opt/conda/bin/Rscript

#' 质谱流式细胞注释
#' 
#' @param celltype 细胞类型数据
#' @param clusterdata 聚类数据
#' @param highlimit 高表达限制值
#' @param lowlimit 低表达限制值
#' @param dimlimit 中表达限制值
#'
#' @export 
cellannotation <- function(celltype,
                           clusterdata,
                           highlimit = 0.3,
                           lowlimit = 0.2,
                           dimlimit = c(0.5,0.2)){
  print("细胞类型去重中")
  celltype <- celltype[!duplicated(celltype[,1]),]
  celltype <- celltype[!duplicated(celltype[,-1]),]
  clusterscore <- clusterdata[,0,drop=F]
  
  print("细胞类型打分中")
  for(i in 1:dim(celltype)[1]){
    clusterscore[,celltype[i,1]] <- 0
    for (j in 1:dim(clusterdata)[1]) {
      score <- 0
      for (k in colnames(celltype)[-1]) {
        ant <- celltype[i,k]
        if(is.na(ant)){next}
        ant <- unlist(strsplit(ant,split = ","))
        num <- length(ant)
        ant <- ant[ant %in% colnames(clusterdata)]
        newnum <- length(ant)
        score <- score-(num-newnum)*0.01
        if(length(ant) == 0){next}
        
        for (antname in ant) {
          if(k == "high"){
            # antscore <- clusterdata[j,antname]-highlimit
            usescore <- quantile(clusterdata[,antname],probs = 0.75)*0.8
            if(usescore < highlimit){usescore <- highlimit}
            antscore <- clusterdata[j,antname]-usescore
          }else if(k == "low"){
            # antscore <- lowlimit-clusterdata[j,antname]
            usescore <- quantile(clusterdata[,antname],probs = 0.25)*1.25
            if(usescore > lowlimit){usescore <- lowlimit}
            antscore <- usescore-clusterdata[j,antname]
          }else if(k == "dim"){
            # antscore <- abs(clusterdata[j,antname]-dimlimit[1])-dimlimit[2]
            usescore <- quantile(clusterdata[,antname],probs = 0.5)
            if(abs(usescore-dimlimit[1]) > dimlimit[2]){usescore <- dimlimit[1]}
            antscore <- abs(clusterdata[j,antname]-quantile(clusterdata[,antname],probs = 0.5))-dimlimit[2]
          }else{
            antscore <- 0
          }
          
          if(antscore < 0){
            antscore <- -10
          }
          score <- score+antscore
        }
      }
      clusterscore[j,celltype[i,1]] <- score
    }
  }
  
  print("细胞类型打分整理中")
  clustertype <- clusterdata[,0,drop=F]
  clustertype[,"cluster_id"] <- paste0("cluster",rownames(clustertype)) 
  clustertype[,"score"] <- apply(clusterscore, 1, max)
  clustertype[,"celltype"] <- apply(clusterscore, 1, function(data){names(which.max(data))})
  clustertype[clustertype[,"score"] <= 0,"celltype"] <- "unkown" 
  
  return(clustertype)
}

#' SingleCellExperiment对象操作之：人工注释细胞类型
#' 
#' @param sceObject sceObject对象rds名
#' @param saveplotpath 细胞注释结果保存路径，默认在4.celltype目录下
#' @param highlimit 高表达限制值
#' @param lowlimit 低表达限制值
#' @param dimlimit 中表达限制值
#' @param maxK 分群数
#' @param celltypedir 读取整理好的celltype表格，默认在rawdata目录下
#' @param clusterdatadir 抗体在每个簇的表达热图表格，默认在2.cluster目录下
#' @param mapname1 注释结果可视化保存的图片名，默认"manual_annotation_celltype"
#' @param mapname2 注释结果分面图，默认"manual_annotation_celltype_facet"
#' @param mapname3 注释结果堆叠柱状图，默认"sce_celltype_ratio"
#' @param savexlsxname 注释结果保存的表格名，默认"cellanno.xlsx"
#' @param scename 细胞注释后的sceObject结果保存名，默认"flowSet_object_sce_reducedDim_celltype.rds"
#' @param family 字体
#' @param width 图片宽度
#' @param height 图片长度
#' @param dpi 图片分辨率
#' 
#' @export
sce_anno <- function(sceObject = "rawdata/flowSet_object_sce_cluster_reducedDim.rds",
                     saveplotpath = "5.Celltype/",
                     highlimit = 0.3,
                     lowlimit = 0.2,
                     dimlimit = c(0.5, 0.2),
                     maxK = 20,
                     celltypedir = "rawdata/celltype.xlsx",
                     clusterdatadir = "2.Cluster/heatmap_info.xlsx",
                     scename = "flowSet_object_sce_reducedDim_celltype.rds",
                     samplecolor = SelectColors("customecol2"),
                     ...){
  
  suppressMessages(library("CATALYST"))
  suppressMessages(library("cytofWorkflow"))
  suppressMessages(library("ggplot2"))
  suppressMessages(library("HDCytoData"))
  
  sceObject <- readdata(filename = sceObject)
  
  celltype <- readdata(celltypedir)
  clusterdata <- readdata(clusterdatadir)
  clustertype <- cellannotation(celltype = celltype,
                                clusterdata = clusterdata,
                                highlimit = highlimit,
                                lowlimit = lowlimit,
                                dimlimit = dimlimit)
  
  merging_table1 <- data.frame(original_cluster=1:maxK, new_cluster=clustertype$celltype)
  cell <- unique(merging_table1$new_cluster)
  cell <- sort(cell, decreasing = F)
  if("unkown" %in% cell){
    cell_1 <- cell[cell!="unkown"]
    cell <- c(cell_1, "unkown")
  }
  merging_table1$new_cluster <- factor(merging_table1$new_cluster, cell)
  
  sceObject <- mergeClusters(sceObject, k = paste0("meta", maxK),
                             table = merging_table1, id = "merging1",overwrite = TRUE)
  
  plot1 <- plotDR(sceObject, "TSNE", k_pal=samplecolor, color_by = "merging1") + 
    labs(color="celltype") + guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  plot1_1 <- plotDR(sceObject, "UMAP", k_pal=samplecolor, color_by = "merging1") + 
    labs(color="celltype") + guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  # 分面图
  plot2 <- plotDR(sceObject, "TSNE", k_pal=samplecolor, color_by = "merging1", facet_by = ifelse(length(unique(sceObject$group)) == 1,"sample","group")) + 
    labs(color="celltype") + guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  plot2_1 <- plotDR(sceObject, "UMAP", k_pal=samplecolor, color_by = "merging1", facet_by = ifelse(length(unique(sceObject$group)) == 1,"sample","group")) + 
    labs(color="celltype") + guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  # 堆叠柱状图
  plot3 <- plotAbundances(sceObject, k = "merging1", k_pal=samplecolor, by = "sample_id") + labs(color="celltype") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) + 
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  df3 <- data.frame(plot3$data)
  colnames(df3) <- c("sample", "celltype", "Freq(%)", "group")
  plot3$data$cluster_id -> df3$celltype
  plot3$data$sample_id -> df3$sample
  
  ggplotsave(plot = plot1, savepath = saveplotpath, mapname = "manual_annotation_celltype",...)
  ggplotsave(plot = plot2, savepath = saveplotpath, mapname = "manual_annotation_celltype_facet",...)
  ggplotsave(plot = plot1_1, savepath = saveplotpath, mapname = "manual_annotation_celltype_umap",...)
  ggplotsave(plot = plot2_1, savepath = saveplotpath, mapname = "manual_annotation_celltype_facet_umap",...)
  ggplotsave(plot = plot3, savepath = saveplotpath, mapname = "sce_celltype_ratio",...)
  
  savexlsx1(data = clustertype,
            filename = paste0(saveplotpath, "/cellanno.xlsx"),
            sheet = "cellanno")
  savexlsx1(data = df3,
            filename = paste0(saveplotpath, "/cellanno_ratio.xlsx"),
            sheet = "cellanno_ratio")
  
  # save sce_cellannotion
  # print(paste0(date(), " ", "prepare to save RDS..."))
  saverds(data = sceObject,filename = paste0(dirname(celltypedir),"/",scename))
  # print("save sceObject successfully")
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-sce","--sceObject",default = "rawdata/flowSet_object_sce_cluster_reducedDim.rds", help = "读取sce对象的rds文件名")
  parser$add_argument("-max","--maxK",default = 20, type = "integer", help = " 分群数（聚类数）")
  parser$add_argument("-cell","--celltypedir",default = "rawdata/celltype.xlsx", help = "读取整理好的celltype表格，默认在rawdata目录下，命名celltype.xlsx")
  parser$add_argument("-clust","--clusterdatadir",default = "2.Cluster/heatmap_info.xlsx", help = "抗体在每个簇的表达热图表格heatmap_info.xlsx，默认在2.cluster目录下")
  parser$add_argument("-high","--highlimit",default = 0.3, type = "double", help = "高表达限制值")
  parser$add_argument("-low","--lowlimit",default = 0.2, type = "double", help = "低表达限制值")
  parser$add_argument("-dim","--dimlimit",default = c(0.5,0,2), type = "double", nargs = 2,help = "中表达限制值")
  
  parser$add_argument("-sn","--scename",default = "flowSet_object_sce_reducedDim_celltype.rds", help = "细胞注释后的sceObject结果保存名，默认'flowSet_object_sce_reducedDim_celltype.rds")
  parser$add_argument("-sp","--saveplotpath",default = "5.Celltype/", help = "细胞注释结果保存路径，默认在4.celltype目录下")
  parser$add_argument("-wi","--width", default = 8, type = "integer", help = "图片宽度")
  parser$add_argument("-he","--height", default = 8, type = "integer", help = "图片长度")
  parser$add_argument("-i","--dpi", default = 300, type = "integer", help = "图片分辨率")
  parser$add_argument("-fa","--family", default = "sans", help = "字体")
  
  args <- parser$parse_args()
  
  cellargs <- do.call(what = sce_anno, args = args)
}
