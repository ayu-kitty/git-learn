#' 处理PD标记下机数据-单标
#' @export
pdfetch<-function(inputpath="./projectdata/",savepath="./rawdata/"){
  label<-read.xlsx(paste0(inputpath,"样品标记对照表.xlsx"))
  #蛋白表
  profil<-dir(inputpath,pattern = "*pro*.xlsx")
  pro<-readdata(paste0(inputpath,profil))
  if(file.exists("./background/Annotation_Data.xlsx")){
    anno<-readxlsx("./background/Annotation_Data.xlsx")[,-1]
    names(anno)[1]<-"Accession"
    pro[,"Gene Name"]<-left_join(pro,anno,by="Accession")$Gene_Name
  }else pro[,"Gene Name"]<-NA
  newpro <- select(pro,c("Accession","Gene Name","Description","Coverage [%]","# Peptides","# PSMs","# Unique Peptides","# AAs","MW [kDa]","calc. pI","Score Sequest HT: Sequest HT",grep("^Abundance:",colnames(pro),value=T)))
  
  names(newpro)[grep("^Abundance:",colnames(newpro))]<-label[,1]
  if("Mix"%in%label[,1]){
	  newpro[,grep("Mix",label[,1],invert = T,value = T)]<-newpro[,grep("Mix",label[,1],invert = T,value = T)]/newpro[,"Mix"]
	  newpro<-newpro[,-grep("Mix",names(newpro))]
  }
  savexlsx(newpro,paste0(savepath,"Protein quantitation.xlsx"))
  #肽段表
  pepfil<-dir(inputpath,pattern = "*pep*.xlsx") 
  pep<-readdata(paste0(inputpath,pepfil))
  newpep <- pep[,c("Annotated Sequence","Modifications","Qvality PEP","Qvality q-value","# Protein Groups","# Proteins","# PSMs","Master Protein Accessions","Positions in Master Proteins","# Missed Cleavages","Theo. MH+ [Da]",grep("^Abundance:",colnames(pep),value=T),"Quan Info",grep("^Found in Sample:",colnames(pep),value=T),"Confidence (by Search Engine): Sequest HT","Percolator q-Value (by Search Engine): Sequest HT","Percolator PEP (by Search Engine): Sequest HT","XCorr (by Search Engine): Sequest HT")]
  names(newpep)[grep("^Abundance:",colnames(newpep))]<-label[,1]
  savexlsx(newpep,paste0(savepath,"Peptides.xlsx"))
  ppm<-as.data.frame(matrix(nrow = 1,ncol=5))
  colnames(ppm)<-c("Items","二级谱图总数(MS/MS)","有效谱图数(PSM)","鉴定蛋白数(Protein Groups)","肽段数(Peptides)")
  readdata(paste0(inputpath,"MSMS.xlsx"))->ms
  
  ppm[1,1]<-"FDR小于0.01"
  ppm[1,4]<-nrow(newpro)
  ppm[1,5]<-nrow(newpep)
  ppm[1,2]<-nrow(ms)
  ppm[1,3]<-sum(pep$'# PSMs')
  savexlsx(ppm,paste0(savepath,"Identified_number.xlsx"))
  files<-dir(inputpath)
  mvfil<-setdiff(files,c(profil,pepfil,"MSMS.xlsx"))%>% paste0("'",inputpath,.,"'")%>% paste0(.,collapse=" ")
  system(paste0("cp -r ",mvfil," ",savepath))
}
#' 处理PD标记下机数据-多标
#' @export
pdmultifetch<-function(inputpath="./projectdata/",savepath="./rawdata/"){
  label<-readdata(paste0(inputpath,"样品标记对照表.xlsx"))
  multi_group <- dir(path=inputpath,pattern="^Label set*")
  ppm<-as.data.frame(matrix(nrow = 1,ncol=5))
  crepro<-list()
  colnames(ppm)<-c("Items","二级谱图总数(MS/MS)","有效谱图数(PSM)","鉴定蛋白数(Protein Groups)","肽段数(Peptides)")
  for(i in 1:length(multi_group)){
    createdir(paste0(savepath,multi_group[i]))
    #蛋白表
    profil<-dir(paste0(inputpath,multi_group[i]),pattern = "*pro*.xlsx")
    pro<-readdata(paste0(inputpath,multi_group[i],"/",profil))
    if(file.exists("./background/Annotation_Data.xlsx")){
      anno<-readxlsx("./background/Annotation_Data.xlsx")[,-1]
      names(anno)[1]<-"Accession"
      pro[,"Gene Name"]<-left_join(pro,anno,by="Accession")$Gene_Name
    }else pro[,"Gene Name"]<-NA
    newpro <- select(pro,c("Accession","Gene Name","Description","Coverage [%]","# Peptides","# PSMs","# Unique Peptides","# AAs","MW [kDa]","calc. pI","Score Sequest HT: Sequest HT",grep("^Abundance:",colnames(pro),value=T)))
    
    names(newpro)[grep("^Abundance:",colnames(newpro))]<-label[label[,1]==multi_group[i],2]
	if("Mix"%in%label[label[,1]==multi_group[i],2]){
		  newpro[,grep("Mix",label[label[,1]==multi_group[i],2],invert = T,value = T)]<-newpro[,grep("Mix",label[label[,1]==multi_group[i],2],invert = T,value = T)]/newpro[,"Mix"]
		  newpro<-newpro[,-grep("Mix",names(newpro))]
	}
    filpro<-filter(newpro,`# Unique Peptides`>0,`Score Sequest HT: Sequest HT`>0)
    crepro[[i]]<-select(filpro,c("Accession",grep("Gene Name",names(filpro),value = T),"Description",all_of(grep("Mix",label[label[,1]==multi_group[i],2],invert = T,value = T))))
    savexlsx(newpro,paste0(savepath,multi_group[i],"/Protein quantitation.xlsx"))
    #肽段表
    pepfil<-dir(paste0(inputpath,multi_group[i]),pattern = "*pep*.xlsx")
    pep<-readdata(paste0(inputpath,multi_group[i],"/",pepfil))
    newpep <- pep[,c("Annotated Sequence","Modifications","Qvality PEP","Qvality q-value","# Protein Groups","# Proteins","# PSMs","Master Protein Accessions","Positions in Master Proteins","# Missed Cleavages","Theo. MH+ [Da]",grep("^Abundance:",colnames(pep),value=T),"Quan Info",grep("^Found in Sample:",colnames(pep),value=T),"Confidence (by Search Engine): Sequest HT","Percolator q-Value (by Search Engine): Sequest HT","Percolator PEP (by Search Engine): Sequest HT","XCorr (by Search Engine): Sequest HT")]
    names(newpep)[grep("^Abundance:",colnames(newpep))]<-label[label[,1]==multi_group[i],2]
    savexlsx(newpep,paste0(savepath,multi_group[i],"/Peptides.xlsx"))
    
    readdata(paste0(inputpath,multi_group[i],"/MSMS.xlsx"))->ms
    
    ppm[i,1]<-multi_group[i]
    ppm[i,4]<-nrow(newpro)
    ppm[i,5]<-nrow(newpep)
    ppm[i,2]<-nrow(ms)
    ppm[i,3]<-sum(pep$'# PSMs')
  }
  crepronew<-crepro[[1]]
  for(j in 2:length(crepro)){
    crepronew<-inner_join(crepronew,crepro[[j]], by = c('Accession',grep("Gene Name",names(crepro[[1]]),value = T),"Description")) 
  }
  write.xlsx(crepronew,paste0(savepath,"Protein quantitation.xlsx"))
  savexlsx(ppm,paste0(savepath,"Identified_number.xlsx"))
  mvfil<-setdiff(list.files(inputpath),list.dirs(inputpath,full.names = F,recursive = F))%>% paste0("'",inputpath,.,"'")%>% paste0(.,collapse=" ")
  system(paste0("cp -r ",mvfil," ",savepath))
}
