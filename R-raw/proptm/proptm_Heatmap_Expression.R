#!/opt/conda/bin/Rscript

#' 表达量聚类热图
#'
#' @param savepath 
#' @param inputpath 
#' @param inputfile 
#' @param imagetype 
#' @param height 
#' @param width 
#' @param dpi 
#' @param fontfamily 
#' @param units 
#' @param mapname 
#' @param ... 
#'
#' @export
heatmap_Expre <- function(savepath = "Heatmap_Expression/", inputpath="./",inputfile="去除本底差异数据.xlsx",
                          rawpath = "../../../rawdata/",imagetype = c("pdf", "png"),
                        height = 12, width = 10, dpi = 300, fontfamily = "sans",
                        units = "in",mapname = 'heatmap',...){
  pacman::p_load(dplyr,pheatmap,stringr)
  samp <- readdata(paste0(rawpath,"Sample_Group.xlsx"),sheet=1)
  grp <- readdata(paste0(rawpath,"Sample_Group.xlsx"),sheet=2)
  classfile = paste0(rawpath,"classtype-ptm.xlsx")
  cg <- getsheetname(paste0(inputpath,inputfile))
  for(i in 1:length(cg)){
    alldata <- readdata(paste0(inputpath,inputfile),sheet = cg[i])
    #gsub("-vs-","_",cg[i]) %>% strsplit(.,"_") %>% unlist(.) -> compa
    grp_na<- gsub('-vs-','/',cg[i])
    if(grp_na%in%grp$比较组){
      strsplit(grp_na,"/") %>% unlist(.) -> compa
    }else{
      strsplit(grp$比较组[which(grp_na==grp$new_group)],"/") %>% unlist(.) -> compa
    }
    sams <- samp[samp$Group %in% compa,]
    annotation_col <- data.frame(Group = sams$Group)
    rownames(annotation_col) <- sams$Sample
    ncol <- which(names(alldata) == "class")
    ptmn<- c("_Phosphoprotein","_Acetylprotein","_Ubiquprotein","_Deamidprotein")
    sitefc <- intersect(ptmn,gsub("FC","",names(alldata)[grep("FC",names(alldata))]))
    names <- sams$Sample
    alldata <- alldata[, c(1,(ncol+1):(ncol+2*nrow(sams)))]
    casename <- samp[samp$Group %in% compa[1],"Sample"]
    controlname <- samp[samp$Group %in% compa[2],"Sample"]
    sitedata <- select(alldata,matches(sitefc))
    names(sitedata) <- gsub(sitefc,"",names(sitedata))
    prodata <- select(alldata,matches("_protein"))
    names(prodata) <- gsub("_protein","",names(prodata))
    heatdata <- data.frame(prosite=alldata$prosite)
    for (j in 1:length(names)) {
      sitesample <- sitedata[,names[j]]
      prosample <- prodata[,names[j]]
      heatdata[,names[j]] <- sitesample-prosample
    }
    #绘图
    heatmapdata <- heatdata[-1]
    rownames(heatmapdata)<-heatdata$prosite
    heatmapdata <- heatmapdata[,c(controlname,casename)]
    if(nrow(heatmapdata)!=0){
      cellwidth = ifelse(600 / dim(heatmapdata)[2] > 40, 40, 600 / dim(heatmapdata)[2])
      cellheight = ifelse(600 / dim(heatmapdata)[1] > 15, 15, 600 / dim(heatmapdata)[1])
      treeheight_row = 40 + dim(heatmapdata)[1] * 0.05
      fontsize_col = ifelse(0.5 * cellwidth > 15, 15, 0.5 * cellwidth)
      fontsize_row = 0.75 * cellheight
      p <- pheatmap(heatmapdata,
                    treeheight_row = treeheight_row, # 聚类树的列长度
                    treeheight_col = 40 + dim(heatmapdata)[2] * 0.05, # 聚类树行长度
                    scale = ifelse(nrow(heatmapdata)>2,"row","none"), # 矩阵有没有进行标准化
                    cluster_cols = F, # 按列聚类
                    cluster_rows = ifelse(nrow(heatmapdata)>1,T,F), # 按行聚类
                    fontsize_row = fontsize_row, # 行字体大小
                    fontsize_col = fontsize_col, # ；列字体大小
                    cellwidth = cellwidth, 
                    cellheight = cellheight, 
                    annotation_col = annotation_col,
                    annotation_colors = {
                      list(Group = stylefun_group(classfile = classfile,
                                                  object = annotation_col,
                                                  styletype = "fill"))},
                    color = SelectColors(palette = "heatmapcol",n = 255), # 设置渐变色
                    border_color = F, # 每个小块间是否要用颜色分开
                    show_rownames = T,
                    show_colnames = T,
                    filename = NA,
                    angle_col = 90)
      ggplotsave(plot = p,
                 mapname = paste0(mapname,"-",cg[i]),
                 height = (dim(heatmapdata)[1] * cellheight + 200 + max(nchar(names(heatmapdata))) * 0.5 * fontsize_col)/72,
                 width =(dim(heatmapdata)[2] * cellwidth + 200 + max(nchar(row.names(heatmapdata))) * 0.5 * fontsize_row + treeheight_row) / 72,
                 savepath = paste0(savepath),
                 imagetype = imagetype,
                 dpi = dpi,
                 family = fontfamily,
                 units = units)
      mat_scale <- data.frame(t(apply(heatmapdata, 1, scale)))
      colnames(mat_scale) <- colnames(heatmapdata)
      if(nrow(mat_scale)>1){
        mat_scale <- mat_scale[p$tree_row$order,]}
      mat_scale <- cbind(rownames(mat_scale),mat_scale)
      colnames(mat_scale)[1] <- "prosite"
      heatlist <- list(heatdata,mat_scale)
      names(heatlist) <- c("聚类热图数据","绘图数据")
      savexlsx(heatlist,paste0(savepath,mapname,"-",cg[i],".xlsx"))
    }
  }
  if (file.exists(paste0(inputpath,"Rplots.pdf"))) {
    file.remove(paste0(inputpath,"Rplots.pdf"))
  }
}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 10, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 12, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-if","--inputfile", default = "去除本底差异数据.xlsx",help = "显著差异数据文件名")
  parser$add_argument("-m","--mapname",default = 'heatmap', help = "图片名称，默认为heatmap")
  parser$add_argument("-ip","--inputpath", default = "./",help = "输入文件路径，默认为当前路径")
  parser$add_argument("-sp","--savepath",default = "Heatmap_Expression/", help = "热图保存图片路径，默认为./heatmap-expression/")
  parser$add_argument("-rp","--rawpath", default = "../../../rawdata/",help = "rawdata路径，该目录下需有Sample_Group和classtype-ptm.xlsx")
  args <- parser$parse_args()
  Heatmap_count <- do.call(heatmap_Expre,args = args)
}