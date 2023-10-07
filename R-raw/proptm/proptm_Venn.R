#!/opt/conda/bin/Rscript

#' 韦恩图
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
#' @param ptm 位点表
#' @param col 颜色
#' @param mapname 图片名
#' @param ... 
#' @export
ptmc<-function(colna){
  if("Phospho (STY) Probabilities"%in% colna){
    ptmz<-"磷酸化"
  }else if("Acetyl (K) Probabilities"%in% colna){
    ptmz<-"乙酰化"
  }else if("GlyGly (K) Probabilities"%in% colna){
    ptmz<-"泛素化"
  }else if("Deamidation 18O (N) Probabilities"%in% colna){
    ptmz<-"糖基化"
  }
  ptmz
}
#' @export
ptma<-function(colna){
  if("Phospho (STY) Probabilities"%in% colna){
    ptmn<-"Phospho"
  }else if("Acetyl (K) Probabilities"%in% colna){
    ptmn<-"Acetyl"
  }else if("GlyGly (K) Probabilities"%in% colna){
    ptmn<-"Ubiqu"
  }else if("Deamidation 18O (N) Probabilities"%in% colna){
    ptmn<-"Deamid"
  }
  ptmn
}
#' @export
venn_count<-function(savepath = "Venn/", incompare = "sheet1", imagetype = c("pdf", "png"),
                     height = 10, width = 10, dpi = 300, fontfamily = "sans",units = "in",input = "./rawdata/",
                     pro = "表达矩阵-pro.xlsx",ptm = "表达矩阵-ptm.xlsx",col = c("#4ea8de","#e4c1f9"),
                     mapname = 'credit_venn',...){
  pacman::p_load(dplyr,ggplot2,stringr,VennDiagram,lmbio)
  ##可信
  protein<-readdata(paste0(input,pro),sheet = incompare)
  site<-readdata(paste0(input,ptm),sheet = incompare)
  ptmn<-ptma(colna=colnames(site))
  ptmz<-ptmc(colnames(site))
  type<- c('修饰类型'=ptmn,'修饰'=ptmz)
  file.remove(paste0(input,'ptm_type.txt'))
  savetxt(data = type,
           filename = paste0(input,"/ptm_type.txt"),append = T)
  ##venn图
  list <- list(protein$Accession,unique(site$Accession))
  names(list) <- c("Proteins",paste0(ptmn,"proteins"))
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") 
  p <- venn.diagram(list,filename = NULL,
                    output=TRUE,fill=col,margin=0.1,alpha = 0.5,
                    col="transparent",fontfamily = "ArialMT",
                    cat.cex=1.5,cat.dist = c(0.055, 0.085),cat.fontfamily ="ArialMT")
  unlink(list.files(pattern = "^VennDiagram.*\\.log"))
  ggplotsave(plot = p,
             mapname = mapname,
             height = height,
             width =width,
             savepath = savepath,
             imagetype = imagetype,
             dpi = dpi,
             family = fontfamily,
             units = units)
  interpro<-intersect(protein$Accession,unique(site$Accession))
  sta<-data.frame(t(c(nrow(protein),length(unique(site$Accession)),length(interpro))))
  names(sta)<-c("蛋白数量",paste0(ptmz,"蛋白数量"),"关联蛋白数量")
  write.table(sta,paste0(savepath,"关联蛋白统计.txt"),sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
  prolist<-list(protein,site,protein[protein$Accession%in%interpro,],site[site$Accession%in%interpro,])
  names(prolist) <- c("Proteins",paste0(ptmn,"proteins"),"关联蛋白-pro",paste0("关联蛋白-",ptmn))
  savexlsx(prolist,paste0(savepath,"关联蛋白.xlsx"))    
} 
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 10, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 10, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-pr","--pro", default="表达矩阵-pro.xlsx", help="输入蛋白表，默认为表达矩阵-pro.xlsx")
  parser$add_argument("-pt","--ptm", default="表达矩阵-ptm.xlsx", help="输入位点表，默认为表达矩阵-ptm.xlsx")
  parser$add_argument("-ip","--input", default="./rawdata/", help="读取数据路径，默认为./rawdata/")
  parser$add_argument("-ic","--incompare",default = "sheet1", help = "比较组，也是输入文件的sheet名")
  parser$add_argument("-c","--col",default = c("#4ea8de","#e4c1f9"), help = "图片颜色",nargs="+")
  parser$add_argument("-m","--mapname",default = 'credit_venn', help = "图片名称，默认为venn")
  parser$add_argument("-s","--savepath",default = "Venn/", help = "venn保存图片路径，默认为Venn/")
  args <- parser$parse_args()
  Venn_count <- do.call(venn_count,args = args)
}