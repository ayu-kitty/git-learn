#!/opt/conda/bin/Rscript
#' 激酶活性多比较组时绘制热图
#'
#' @param savepath 保存路径
#' @param m.cutoff 底物数量
#' @param stats p值
#' @param imagetype 图片类型
#' @param height 图片高度
#' @param width 图片宽度
#' @param dpi 分辨率
#' @param fontfamily 字体
#' @param ... 
#' @export
KSEA_pheatmap<-function(savepath = "./KSEA_Heatmap/", inputpath = './KSEA_Activity/',
                        m.cutoff=2,stats='p.value' ,imagetype = c("pdf", "png") , 
                        dpi = 300, fontfamily = "sans", ...){
  
  pacman::p_load(dplyr,stringr,Hmisc,readxl,openxlsx,KSEAapp,ggplot2,optparse,pheatmap,tidyr,lmbio)
  
  con<- dir(path = inputpath, pattern = "*")
  createdir(savepath)
  if(length(con)>1){
    KS<- list()
    for (g in 1:length(con)) {
      KS[[g]]<-read.xlsx(paste0(inputpath,con[g],"/KSEA Kinase Scores.xlsx"))
    }
    sample.labels<- con
    
    filter.m = function(dataset, m.cutoff) {
      filtered = dataset[(dataset$m >= m.cutoff), ]
      return(filtered)
    }
    KS.m = lapply(KS, function(...) filter.m(...,m.cutoff))
    for (i in 1:length(KS.m)) {
      names = colnames(KS.m[[i]])[c(2:7)]
      colnames(KS.m[[i]])[c(2:7)] = paste(names, i,sep = ".")
    }
    master = Reduce(function(...) merge(..., by = "Kinase.Gene",all = F), KS.m)
    if(!nrow(master)==0){
      row.names(master) = master$Kinase.Gene
      columns = as.character(colnames(master))
      merged.scores = as.matrix(master[, grep("z.score", columns)])
      colnames(merged.scores) = sample.labels
      merged.stats = as.matrix(master[, grep(stats, columns)])
      data_mark=merged.stats
      for(i in 1:nrow(merged.stats)){
        for(j in 1:ncol(merged.stats)){
          if (merged.stats[i,j]<0.0001){data_mark[i,j]="****"}else{
            if (merged.stats[i,j]<0.001){data_mark[i,j]="***"}else{
              if (merged.stats[i,j]<0.01){data_mark[i,j]="**"}else{
                if (merged.stats[i,j]<0.05){data_mark[i,j]="*"}else{
                  data_mark[i,j]=""}}}}}}
      cellwidth =  ifelse(600 / dim(merged.scores)[2] > 40, 40, 600 / dim(merged.scores)[2])
      cellheight = ifelse(600 / dim(merged.scores)[1] > 15, 15, 600 / dim(merged.scores)[1])
      fontsize_col = ifelse(0.5 * cellwidth > 15, 15, 0.5 * cellwidth)
      fontsize_row = 0.75 * cellheight
      p<-pheatmap(merged.scores,display_numbers=data_mark,
                   cluster_rows = F,
                   cluster_cols = F,
                   clustering_method = "complete",
                   fontfamily= fontfamily,
                   fontsize_number=14,
                   number_format ="%.3f",
                   number_color = "black",
                   show_colnames = T,
                   show_rownames = T,
                   angle_col=90,
                   color =  SelectColors(palette = "heatmapcol",n = 255),
                   cellwidth =  cellwidth,
                   cellheight = cellheight,
                   scale = "none",
                   fontsize_col = fontsize_col,
                   fontsize_row = fontsize_row)
      height = (dim(merged.scores)[1] * cellheight +  250 + max(nchar(colnames(merged.scores))) * (0.5 * fontsize_col)) / 72
      width = (dim(merged.scores)[2] * cellwidth + 250+ max(nchar(row.names(merged.scores))) * (0.5 * fontsize_row)) / 72
      ggplotsave(plot = p,
               mapname = "KSEA.Merged.Heatmap",
               height = height,
               width = width,
               savepath = savepath,
               imagetype = imagetype,
               dpi = dpi,
               family = fontfamily,
               units = "in")
      ms<-cbind(as.data.frame(rownames(merged.scores)),as.data.frame(merged.scores))
      colnames(ms)=c("kinase",paste0(sample.labels,"-z.scores"))
      write.xlsx(ms,paste0(savepath,"KSEA热图数据.xlsx"))
    }else{
      savetxt(data = "无交集激酶，不生成交集激酶热图。",
              filename = paste0(savepath,"/说明.txt"),append = T)
    }
    
    ###并集蛋白
    KS.m_nrow = lapply(KS.m, function(...)nrow(...))%>%unlist()
    if(0%in%KS.m_nrow){
      KS.m2 = KS.m[lapply(KS.m,nrow)>0]
      master2 = Reduce(function(...) merge(..., by = "Kinase.Gene",all = T), KS.m2)
      sample.labels2 = sample.labels[-(which(KS.m_nrow==0))]
      savetxt(data = "该项目有比较组底物数量均<2，因此筛选后激酶数量为0，绘图时删去该比较组。",
              filename = paste0(savepath,"/说明.txt"),append = T)
    }else{
      master2 = Reduce(function(...) merge(..., by = "Kinase.Gene",all = T), KS.m)
      sample.labels2 = sample.labels
    }
    row.names(master2) = master2$Kinase.Gene
    columns = as.character(colnames(master2))
    merged.scores2 = as.matrix(master2[, grep("z.score", columns)])
    colnames(merged.scores2) = sample.labels2
    merged.stats2 = as.matrix(master2[, grep(stats, columns)])
    data_mark2=merged.stats2
    for(i in 1:nrow(merged.stats2)){
      for(j in 1:ncol(merged.stats2)){
        if (is.na(merged.stats2[i,j])){data_mark2[i,j]=""}else{
          if (merged.stats2[i,j]<0.0001){data_mark2[i,j]="****"}else{
            if (merged.stats2[i,j]<0.001){data_mark2[i,j]="***"}else{
              if (merged.stats2[i,j]<0.01){data_mark2[i,j]="**"}else{
                if (merged.stats2[i,j]<0.05){data_mark2[i,j]="*"}else{
                  data_mark2[i,j]=""}}}}}}}
    cellwidth2 =  ifelse(600 / dim(merged.scores2)[2] > 40, 40, 600 / dim(merged.scores2)[2])
    cellheight2 = ifelse(600 / dim(merged.scores2)[1] > 15, 15, 600 / dim(merged.scores2)[1])
    fontsize_col2 = ifelse(0.5 * cellwidth > 15, 15, 0.5 * cellwidth2)
    fontsize_row2 = 0.75 * cellheight2
    p2<-pheatmap(merged.scores2,display_numbers=data_mark2,
                 cluster_rows = F,
                 cluster_cols = F,
                 clustering_method = "complete",
                 fontfamily= 'sans',
                 fontsize_number=ifelse(nrow(merged.scores2)>25,680/nrow(merged.scores2),ifelse(nrow(merged.scores2)<10,250/nrow(merged.scores2),360/nrow(merged.scores2))),
                 number_format ="%.3f",
                 number_color = "black",
                 show_colnames = T,
                 show_rownames = T,
                 angle_col=90,
                 color = SelectColors(palette = "heatmapcol",n = 255),
                 cellwidth =  cellwidth2,
                 cellheight = cellheight2,
                 scale = "none",
                 fontsize_col = fontsize_col2,
                 fontsize_row = fontsize_row2)
    height2 = (dim(merged.scores2)[1] * cellheight2 +  250 + max(nchar(colnames(merged.scores2))) * (0.5 * fontsize_col2)) / 72
    width2 = (dim(merged.scores2)[2] * cellwidth2 + 250+ max(nchar(row.names(merged.scores2))) * (0.5 * fontsize_row2)) / 72
    ggplotsave(plot = p2,
               mapname = "KSEA.Merged_union.Heatmap",
               height = height2,
               width = width2,
               savepath = savepath,
               imagetype = imagetype,
               dpi = dpi,
               family = fontfamily,
               units = "in")
    ms2<-cbind(as.data.frame(rownames(merged.scores2)),as.data.frame(merged.scores2))
    colnames(ms2)=c("kinase",paste0(sample.labels2,"-z.scores"))
    write.xlsx(ms2,paste0(savepath,"KSEA热图数据-union.xlsx"))
    
    write.xlsx(KS, file = paste0(savepath,"phos_kinase-source.xlsx"),sheetName=sample.labels)
  }else{
    savetxt(data = "该项目为单个比较组，不进行热图绘制",
            filename = paste0(savepath,"/说明.txt"),append = T)
    return()
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
  # 此图参数
  parser$add_argument("-s","--savepath",default = "./KSEA_Heatmap/", help = "绘图存放分析文件夹路径，默认./KSEA_Heatmap/")
  parser$add_argument("-ip","--inputpath",default = "./KSEA_Activity/", help = "绘图数据导入路径，默认./KSEA_Activity/")
  parser$add_argument("-m","--m.cutoff",default=2,help="激酶对应底物的最小数量,默认2")
  parser$add_argument("-v","--stats",default='p.value',help="富集展示方法,默认'p.value'")
  args <- parser$parse_args()
  KSEA_pheatmap <- do.call(KSEA_pheatmap,args = args)
}
