#!/opt/conda/bin/Rscript
#' 激酶活性分析
#'
#' @param savepath 保存路径
#' @param incompare 比较组名称
#' @param or 物种
#' @param phosfile 磷酸化差异表格
#' @param m.cutoff 底物数量
#' @param p.cutoff P值
#' @param imagetype 图片类型
#' @param height 图片高度
#' @param width 图片宽度
#' @param dpi 图片分辨率
#' @param fontfamily 字体
#' @param ... 
#' @export
KSEA_Activity<-function(savepath = "KSEA_Activity/",incompare = 'A-vs-B', species='hsa' ,inputpath = './' ,fc= 2,p.value=0.05,
                    phosfile = "差异表达矩阵(未筛选)-ptm.xlsx" , m.cutoff=2 , p.cutoff=0.05 ,
                    imagetype = c("pdf", "png") , height = 12, width = 14, dpi = 300, fontfamily = "sans", ...){
  pacman::p_load(dplyr,stringr,Hmisc,readxl,openxlsx,KSEAapp,ggplot2,optparse,pheatmap,tidyr)
  
  if(species=="hsa"){
    org<-"human"
  }else if(species=="mmu"){
    org<-"mouse"
  }else org<-"rat"
  print(paste0("~物种：",org))
  print(paste0("~激酶对应底物的最小数量：",m.cutoff))
  print(paste0("~p值：",p.cutoff))

  backdata<- read_excel("/data/hstore3/public/propip/Kinase_Substrate_Dataset-2023.xlsx")
  KSData<- backdata[which(backdata$SUB_ORGANISM==org&backdata$KIN_ORGANISM==org),]
  
  path=paste0(savepath,incompare)
  createdir(path)
  diff_table<- read_excel(paste0(inputpath,gsub('.*/','',phosfile)),sheet = incompare)
  
  if(file.exists(paste0(inputpath,"ratio.xlsx"))){
    data_fc = readxlsx(paste0(inputpath,'ratio.xlsx'))
    fc = data_fc$FC_ptm
    p.value = data_fc$`P-value_ptm`
  }else print("~！！！！rawdata路径下不存在ratio.xlsx文件，请输入FC和p值参数！！！")
  
  if ('FoldChange'%in%colnames(diff_table)) {
    diff_site<- diff_table[intersect(which(diff_table$`p-value`<p.value),
                                       union(which(diff_table$FoldChange>=fc),which(diff_table$FoldChange<=(1/fc)))),]
    colname<- c("ID","Gene Name","log2FoldChange","p-value")
    diff_phos<- diff_site[,colnames(diff_site)%in%colname]
    phos<- separate(diff_phos,ID,into = c("Protein","site"),sep = ":",remove = T)
  }else{
    diff_site<- diff_table[intersect(which(diff_table$`P-value`<0.05),
                                       union(which(diff_table$FC>=fc),which(diff_table$FC<=(1/fc)))),]
    colname<- c("Accession","Gene Name",'Positions within proteins' ,'Amino acid' ,"log2(FC)","P-value")
    diff_phos<- diff_site[,colnames(diff_site)%in%colname]
    diff_phos$site<- paste0(diff_phos$`Amino acid`,diff_phos$`Positions within proteins`)
    phos<- diff_phos[,c('Accession','site','Gene Name',"log2(FC)","P-value")]
    colnames(phos)<- c('Protein','site','Gene Name',"log2FoldChange","p-value")
  }
  ##计算
  new = phos
  colnames(new)<- c("protein","SUB_MOD_RSD","SUB_GENE","log2FC",'p')
  new = new[complete.cases(new$log2FC), ]
  
  KSData.filtered = KSData
  KSData.dataset = merge(KSData.filtered, new)
  KSData.dataset = KSData.dataset[order(KSData.dataset$GENE),]
  KSData.dataset.abbrev = KSData.dataset[, c("GENE","KINASE","SUB_GENE", "SUB_MOD_RSD","log2FC","p" )]
  colnames(KSData.dataset.abbrev) = c("Kinase.Gene","Kinase.Name","Substrate.Gene","Substrate.Mod","log2FC","p")
  KSData.dataset.abbrev = KSData.dataset.abbrev[order(KSData.dataset.abbrev$Kinase.Gene,
                                                      KSData.dataset.abbrev$Substrate.Gene, KSData.dataset.abbrev$Substrate.Mod,
                                                      KSData.dataset.abbrev$p), ]
  if (nrow(KSData.dataset.abbrev)>1) {
    KSData.dataset.abbrev = aggregate(log2FC ~ Kinase.Gene + Kinase.Name +
                                        Substrate.Gene + Substrate.Mod , data = KSData.dataset.abbrev,
                                      FUN = mean)
  }else{
    KSData.dataset.abbrev = KSData.dataset.abbrev
  }
  
  KSData.dataset.abbrev = KSData.dataset.abbrev[order(KSData.dataset.abbrev$Kinase.Gene),]
  
  kinase.list = as.vector(KSData.dataset.abbrev$Kinase.Gene)
  kinase.list = as.matrix(table(kinase.list))
  if (nrow(kinase.list)>1) {
    Mean.FC = aggregate(log2FC ~ Kinase.Gene, data = KSData.dataset.abbrev,FUN = mean)##log2FC的计算方式为所有激酶对应底物的log2FC均值
  }else{
    Mean.FC = KSData.dataset.abbrev[,c('Kinase.Gene','log2FC')]
  }
  
  Mean.FC = Mean.FC[order(Mean.FC[, 1]), ]
  Mean.FC$mS = Mean.FC[, 2]##mS为log2FC
  Mean.FC$Enrichment = Mean.FC$mS/abs(mean(new$log2FC, na.rm = T))%>%as.numeric()##Enrichment的计算方法为mS除以所有激酶的log2FC均值的绝对值
  Mean.FC$m = kinase.list %>%as.numeric()#统计激酶对应的底物数量
  Mean.FC$z.score = (((Mean.FC$mS - mean(new$log2FC, na.rm = T)) *
                        sqrt(Mean.FC$m))/sd(new$log2FC, na.rm = T))%>%as.numeric()##z.score打分的计算方式为：激酶log2FC减所有激酶的log2FC均值乘以激酶底物数量平方根，再除以所有激酶的log2FC的标准差
  Mean.FC$p.value = pnorm(-abs(Mean.FC$z.score))%>%as.numeric()##p.value是计算z.score打分绝对值负值的正态分布
  Mean.FC$FDR = p.adjust(Mean.FC$p.value, method = "fdr")%>%as.numeric()
  savexlsx(KSData.dataset.abbrev, file = paste0(path,"/Kinase-Substrate Links.xlsx"))
  savexlsx(Mean.FC[order(Mean.FC$Kinase.Gene), -which(colnames(Mean.FC)=="mS")],file = paste0(path,"/KSEA Kinase Scores.xlsx"))
  
  #z_score_bar
  KSEA_z_score_bar(savepath = savepath, incompare = incompare, 
              m.cutoff = m.cutoff,p.cutoff = p.cutoff ,fontfamily = "sans")
  #KSEA_bar
  KSEA_Bar(savepath = savepath, incompare = incompare, 
           m.cutoff = m.cutoff,p.cutoff = p.cutoff ,fontfamily = "sans")
  #KSEA_point
  KSEA_point(savepath = savepath, incompare = incompare, height = 10, width = 8,
             m.cutoff = m.cutoff,p.cutoff = p.cutoff ,fontfamily = "sans")
  #KSEA_Num
  KSEA_Number(savepath = savepath, incompare = incompare,
              m.cutoff = m.cutoff,p.cutoff = p.cutoff ,fontfamily = "sans")
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 14, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 9, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")

  # 此图参数
  parser$add_argument("-ic","--incompare",default = "A_B", help = "比较组名称")
  parser$add_argument("-s","--savepath",default = "KSEA_Activity/", help = "绘图存放分析文件夹路径，默认KSEA_Activity/")
  parser$add_argument("-ip","--inputpath",default = "./", help = "输入文件路径，默认./")
  parser$add_argument("-p","--phosfile",default = "差异表达矩阵(未筛选)-ptm.xlsx", help = "磷酸化位点表达矩阵默认差异表达矩阵(未筛选)-ptm.xlsx")
  parser$add_argument("-or","--species",default = "hsa", help = "物种英文名，hsa，mmu，rno,默认hsa")
  parser$add_argument("-m","--m.cutoff",default=2, type= "double",help="激酶对应底物的最小数量,默认2")
  parser$add_argument("-v","--p.cutoff",default=0.05, type= "double",help="富集p值,默认0.05")
  parser$add_argument("-fc","--fc",default=NULL, type= "double",help="此项目fc值,默认1.2")
  parser$add_argument("-pv","--p.value",default=NULL, type= "double",help="此项目p值,默认pvalue<0.05")
  args <- parser$parse_args()
  KSEA_Activity <- do.call(KSEA_Activity,args = args)
}
