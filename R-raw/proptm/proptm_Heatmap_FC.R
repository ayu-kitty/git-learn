#!/opt/conda/bin/Rscript
#' 去除本底差异数据FC热图
#'
#' @param savepath 保存路径
#' @param incompare 比较组名称
#' @param imagetype 保存图片类型
#' @param height 高度
#' @param width 宽度
#' @param dpi 分辨率
#' @param fontfamily 字体类型
#' @param units 
#' @param pro 蛋白表
#' @param mapname 图片名
#' @param ... 
#' @export
heatmap_count<-function(savepath = "Heatmap_FC/", imagetype = c("pdf", "png"),
                     height = 12, width = 8, dpi = 300, fontfamily = "sans",units = "in",input = "./rawdata/",
                     pro = "去除本底差异数据.xlsx",mapname = 'Heatmap-FC',...){
  pacman::p_load(dplyr,pheatmap,stringr)
  incompare <- getsheetname(paste0(input,pro))
  for(i in 1:length(incompare)){
    protein <- readdata(paste0(input,pro),sheet = incompare[i])
    #绘图
    fcname <- colnames(protein)[grep("FC",colnames(protein))]
    pname <- colnames(protein)[grep("P-value",colnames(protein))]
    heatmapdata <- protein[fcname]
    rownames(heatmapdata)<-protein$prosite
    colnames(heatmapdata) <- gsub("FC_","",fcname)
    if(nrow(heatmapdata)!=0){
      p<-pheatmap(log2(heatmapdata),
                  #main = paste0("FC_",incompare[i]),
                  fontsize_col =25,
                  fontsize_row =ifelse(nrow(heatmapdata)>25,580/nrow(heatmapdata),ifelse(nrow(heatmapdata)<10,150/nrow(heatmapdata),260/nrow(heatmapdata))),
                  angle_col = 90,cluster_cols = F,
                  scale = "none",
                  show_rownames = T,
                  cellwidth = 150,
                  color = SelectColors(palette = "heatmapcol",n = 255),
                  cellheight = 600/dim(heatmapdata)[1])
      ggplotsave(plot = p,
                 mapname = paste0(mapname,"-",incompare[i]),
                 height = height,
                 width =width,
                 savepath = savepath,
                 imagetype = imagetype,
                 dpi = dpi,
                 family = fontfamily,
                 units = units)
      heatlist <- list(protein[c("prosite","Accession",fcname,pname)],cbind(protein["prosite"],log2(heatmapdata)))
      names(heatlist) <- c("聚类热图数据","绘图数据")
      savexlsx(heatlist,paste0(savepath,mapname,"-",incompare[i],".xlsx"))
    }
  }
}
  if(this.path::is.main() & !is.null(whereami::thisfile_r())){
    options(warn=-1)
    suppressMessages(library("lmbio"))
    parser <- ArgumentParser()
    # 基本参数
    parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
    parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
    parser$add_argument("-wi","--width",default = 8, type= "double",help = "图片宽度")
    parser$add_argument("-he","--height",default = 12, type= "double",help = "图片高度")
    parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
    
    # 此图参数
    parser$add_argument("-pr","--pro", default="去除本底差异数据.xlsx", help="输入数据，默认为去除本底差异数据.xlsx")
    parser$add_argument("-ip","--input", default="./rawdata/", help="读取数据路径，默认为./rawdata/")
    parser$add_argument("-m","--mapname",default = 'Heatmap-FC', help = "图片名称，默认为Heatmap-FC")
    parser$add_argument("-s","--savepath",default = "Heatmap_FC/", help = "热图保存图片路径，默认为Heatmap_FC/")
    args <- parser$parse_args()
    Heatmap_count <- do.call(heatmap_count,args = args)
  }