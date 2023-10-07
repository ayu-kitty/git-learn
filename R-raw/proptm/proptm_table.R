#!/opt/conda/bin/Rscript
#' 输出数据表格
#' @param savepath 
#' @param input 
#' @param pro 
#' @param ptm 
#' @param ... 
#' @export
satastic_count<-function(savepath = "./",input = "./rawdata/",
                     pro = "差异表达矩阵(未筛选)-pro.xlsx",ptm = "差异表达矩阵(未筛选)-ptm.xlsx",...){
  pacman::p_load(dplyr,stringr,openxlsx,lmbio)
  #修饰名称
  site <- readdata(paste0(input,"表达矩阵-ptm.xlsx"))
  ptmn<-ptma(colnames(site))
  #常规蛋白和位点合并
  crepro <- readdata(paste0(input,"表达矩阵-pro.xlsx"))[c("Accession")]
  creptm <- readdata(paste0(input,"表达矩阵-ptm.xlsx"))[c("ID","Accession")]
  data <- merge(crepro,creptm,by="Accession",all = T)
  data1 <- data
  incompare <- getSheetNames(paste0(input,pro))
  comproname <- c("Accession","log2FoldChange","FoldChange","Regulation","p-value")
  comptmname <- c("ID","log2FoldChange","FoldChange","Regulation","p-value")
  for (i in 1:length(incompare)) {
    compro <- readdata(paste0(input,pro),sheet = incompare[i])[comproname]
    colnames(compro) <- c("Accession",paste0(incompare[i],"_",comproname[-1],"_protein"))
    data1 <- merge(data1,compro,by = "Accession",all.x = T)
    comptm <- readdata(paste0(input,ptm),sheet = incompare[i])[comptmname]
    colnames(comptm) <- c("ID",paste0(incompare[i],"_",comptmname[-1],"_",ptmn))
    data1 <- merge(data1,comptm,by = "ID",all.x = T)
  }
  ann <- readdata("background/Annotation_Data.xlsx")[,-1]
  colnames(ann)[1] <- "Accession"
  out <- merge(data1,ann,by='Accession',all.x = T)
  #保存颜色
  wb <- createWorkbook()
  addWorksheet(wb,sheet = "result")
  for(i in 1:length(incompare)){
    #保存表格添加颜色
    addStyle(wb,sheet = 1,style = createStyle(fgFill = "#f47d8c"),rows=(which(out[paste0(incompare[i],"_Regulation_protein")] == "Up")+1),
             cols = which(names(out)%in%paste0(incompare[i],"_",comproname[c(-1,-5)],"_protein")),gridExpand = T)
    addStyle(wb,sheet = 1,style = createStyle(fgFill = "#727fb5"),rows=(which(out[paste0(incompare[i],"_Regulation_protein")] == "Down")+1),
             cols = which(names(out)%in%paste0(incompare[i],"_",comproname[c(-1,-5)],"_protein")),gridExpand = T)
    addStyle(wb,sheet = 1,style = createStyle(fgFill = "#f47d8c"),rows=(which(out[paste0(incompare[i],"_Regulation","_",ptmn)] == "Up")+1),
             cols = which(names(out)%in%paste0(incompare[i],"_",comptmname[c(-1,-5)],"_",ptmn)),gridExpand = T)
    addStyle(wb,sheet = 1,style = createStyle(fgFill = "#727fb5"),rows=(which(out[paste0(incompare[i],"_Regulation","_",ptmn)] == "Down")+1),
             cols = which(names(out)%in%paste0(incompare[i],"_",comptmname[c(-1,-5)],"_",ptmn)),gridExpand = T)
  }
  addStyle(wb,sheet = 1,style = createStyle(textDecoration = "bold"),rows = 1, cols = 1:ncol(out))
  #addStyle(wb,sheet = 1,style = createStyle(fontName = "Arial"),rows = 1:nrow(out),cols = 1:ncol(out),gridExpand = T)
  writeData(wb, "result", out)
  saveWorkbook(wb,paste0(savepath,"result.xlsx"),overwrite = T)
  #savexlsx(out,paste0(savepath,"result.xlsx"))
}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  
  # 此图参数
  parser$add_argument("-pr","--pro", default="差异表达矩阵(未筛选)-pro.xlsx", help="输入数据，默认为差异表达矩阵(未筛选)-pro.xlsx")
  parser$add_argument("-pt","--ptm", default="差异表达矩阵(未筛选)-ptm.xlsx", help="输入数据，默认为差异表达矩阵(未筛选)-ptm.xlsx")
  parser$add_argument("-ip","--input", default="./rawdata/", help="读取数据路径，默认为./rawdata/")
  parser$add_argument("-s","--savepath",default = "./", help = "数据保存路径，默认为./")
  args <- parser$parse_args()
  Satastic_count <- do.call(satastic_count,args = args)
}
