#!/opt/conda/bin/Rscript
#' @export
diffanno<-function(reportpath="./",...){
  pacman::p_load(openxlsx,dplyr)
  if(file.exists(paste0(reportpath,"/result/2.Qualitative/Annotation/annotation.xlsx")) &
     file.exists(paste0(reportpath,"/result/4.Different_expression/表达矩阵.xlsx")) &
     file.exists(paste0(reportpath,"/result/4.Different_expression/差异表达矩阵.xlsx"))){
    anno<-read.xlsx(paste0(reportpath,"/result/2.Qualitative/Annotation/annotation.xlsx"))
    nwb<-createWorkbook(paste0(reportpath,"/result/4.Different_expression/表达矩阵1.xlsx"))
    com<-getsheetname(paste0(reportpath,"/result/4.Different_expression/表达矩阵.xlsx"))
    difdata <- readxlsx(paste0(reportpath,"/result/4.Different_expression/表达矩阵.xlsx"))
    newdif<-left_join(difdata,anno,by="Accession")
    addsheet(data = newdif,wb = nwb)
    saveWorkbook(nwb, paste0(reportpath,"/result/4.Different_expression/表达矩阵1.xlsx"), overwrite = T)
    file.rename(paste0(reportpath,"/result/4.Different_expression/表达矩阵1.xlsx"),paste0(reportpath,"/result/4.Different_expression/表达矩阵.xlsx"))
    
    nwb<-createWorkbook(paste0(reportpath,"/result/4.Different_expression/差异表达矩阵1.xlsx"))
    com<-getsheetname(paste0(reportpath,"/result/4.Different_expression/差异表达矩阵.xlsx"))
    for(i in 1:length(com)){
      difdata <- readxlsx(paste0(reportpath,"/result/4.Different_expression/差异表达矩阵.xlsx"), sheet = com[i])
      newdif<-left_join(difdata,anno,by="Accession")
      if("FlankingRegion" %in% colnames(newdif) | "Sequence window" %in% colnames(newdif)){
        if("FlankingRegion" %in% colnames(newdif)){
          sfas<-gsub(";.*","",newdif$FlankingRegion)
        }else sfas<-gsub(";.*","",newdif$`Sequence window`)
        acid<-strsplit(sfas,"")%>%as.data.frame()%>%t()
        rownames(acid)<-c(1:nrow(acid))
        if(file.exists(paste0(reportpath,"/result/4.Different_expression/Motif/",com[i],"/momo.tsv"))){
          motif<-read.table(paste0(reportpath,"/result/4.Different_expression/Motif/",com[i],"/momo.tsv"),header = F)
		  if(motif[1,1]=="mod"){
			names(motif)<-motif[1,1:12]
			motif<-motif[-1,]
		  }else names(motif)<-c("mod","motif","regexp","score","fg_match","fg_size","bg_match","bg_size","fg/bg","unadjusted_p-value","tests","adjusted_p-value")
          newdif[,"Motif"]<-modiff(acid,motif)
        }
        
      }
      addsheet(data = newdif,
               wb = nwb,
               sheet = com[i],
               sheetmoudle = addsheet_fc)
      saveWorkbook(nwb, paste0(reportpath,"/result/4.Different_expression/差异表达矩阵1.xlsx"), overwrite = T)
    }
    file.rename(paste0(reportpath,"/result/4.Different_expression/差异表达矩阵1.xlsx"),paste0(reportpath,"/result/4.Different_expression/差异表达矩阵.xlsx"))
    
  }
}
#' @export
modiff<-function(acid,motif){
  M<-as.data.frame(matrix(ncol=nrow(motif),nrow=nrow(acid)))
  for(m in 1:nrow(motif)){
    gregexpr(".",motif$regexp[m],fixed = T)[[1]]->l
    aa<-acid
    aa[,l]<-"."
    p<-sapply(1:nrow(aa),function(x){
      paste(aa[x,],collapse = "")
    })
    M[p%in%motif$regexp[m],m]<-motif$regexp[m]
  }
  sapply(1:nrow(M),function(x){
    paste(M[x,which(M[x,]!="NA")],collapse = ";")
  })
}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  
  parser$add_argument("-rp","--reportpath", default = "./",help = "report文件夹所在路径，默认为当前路径")
  args <- parser$parse_args()
  diffanno <- do.call(diffanno,args = args)
}
