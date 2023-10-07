#' 处理Spectronaut DIA下机数据
#' @export
spfetch<-function(inputpath="./projectdata/",savepath="./rawdata/"){
  #肽段表
  pepfil<-dir(inputpath,pattern = "*pep*.xls")
  pep<-read.delim2(paste0(inputpath,pepfil),check.names = F)
  pep[pep=="NaN"]<-NA
  pep[pep=="Filtered"]<-NA  
  pep[,grep(".Quantity",names(pep))]<-sapply(pep[,grep(".Quantity",names(pep))], as.numeric)
  names(pep)<-names(pep)%>%gsub("PG.|EG.|PEP.","",.)%>%gsub("\\[[0-9]+\\] ","",.)%>%gsub(".Quantity","",.)
  pep<-pep %>% distinct(ModifiedSequence, .keep_all = T)
  savexlsx(pep,paste0(savepath,"Peptides.xlsx"))
  #统计unique肽段
  newpep<-pep[!duplicated(pep$StrippedSequence),]
  unipro<-as.data.frame(table(filter(newpep,toupper(newpep$IsProteinGroupSpecific)=="TRUE")$ProteinGroups))
  names(unipro)<-c("ProteinGroups","Unique peptides")
  #蛋白表
  profil<-dir(inputpath,pattern = "*pro*.xls")
  pro<-read.delim2(paste0(inputpath,profil),check.names = F)
  pro[pro=="NaN"]<-NA
  pro[pro=="Filtered"]<-NA  
  samplecol <-grep(".Quantity",names(pro))
  pro[,samplecol]<-sapply(pro[,samplecol], as.numeric)
  pro[,samplecol][pro[,samplecol]<1]<-NA
  names(pro)<-names(pro)%>%gsub("PG.","",.)%>%gsub("\\[[0-9]+\\] ","",.)%>%gsub(".Quantity","",.)
  pro<-left_join(pro,unipro,by="ProteinGroups")
  pro[,"Unique peptides"][is.na(pro[,"Unique peptides"])]<-0
  pro[,"Accession"]<-gsub(";.*","",pro$ProteinGroups)
  if(file.exists("./background/Annotation_Data.xlsx")){
    anno<-readxlsx("./background/Annotation_Data.xlsx")[,-1]
    names(anno)[1]<-"Accession"
    pro[,"Gene Name"]<-left_join(pro,anno,by="Accession")$Gene_Name
  }else pro[,"Gene Name"]<-NA
  
  pro[,"Description"]<-gsub(";.*","",pro$ProteinDescriptions)
  pro[,"MW [kDa]"]<-as.numeric(gsub(";.*","",pro$MolecularWeight))/1000
  pro[,"Qvalue"]<-as.numeric(pro[,"Qvalue"])
  pro[,"Cscore"]<-as.numeric(pro[,"Cscore"])
  pro[,"Peptides"]<-as.numeric(pro[,"NrOfStrippedSequencesIdentified (Experiment-wide)"])
  newpro<-select(pro,c("Accession","Gene Name","Description","Organisms","MW [kDa]","Qvalue","Cscore","Peptides","Unique peptides",all_of(samplecol)))
  savexlsx(newpro,paste0(savepath,"Protein quantitation.xlsx"))
  
  ppm<-as.data.frame(matrix(nrow = 1,ncol=3))
  colnames(ppm)<-c("Items","蛋白数(Protein Groups)","肽段数(Peptides)")
  ppm[1,1]<-"DIA"
  ppm[1,2]<-nrow(newpro)
  ppm[1,3]<-nrow(pep)
  #DDA表
  ddafil<-dir(inputpath,pattern = "*DDA*")
  if(length(ddafil)!=0){
    dda<-read.delim2(paste0(inputpath,ddafil),check.names = F)
    SP_dda <- dda %>% distinct(ModifiedPeptide, .keep_all = T)
    ll<-data.frame("before"=c("BGSInferenceId","ModifiedPeptide","StrippedPeptide","IsProteotypic","ProteinDescription","Genes","Organisms"),
                   "after"=c("ProteinIDs","ModifiedPeptide","StrippedPeptide","IsProteotypic","ProteinDescription","Gene Name","Organisms"))
    rownames(ll)<-ll[,1]
    seln<-intersect(ll[,1],colnames(SP_dda))
    SP_dda <- select(SP_dda,all_of(seln))
    colnames(SP_dda) <-ll[seln,2]
    irt<-length(grep("iRT",SP_dda[,1]))
    ppm[2,1]<-"DDA"
    ppm[2,2]<-length(unique(SP_dda[,1]))-irt
    ppm[2,3]<-nrow(SP_dda)-irt
    savexlsx(SP_dda,paste0(savepath,"DDA_library.xlsx"))
  }
  files<-dir(inputpath)
  savexlsx(ppm,paste0(savepath,"Identified_number.xlsx"))
  mvfil<-setdiff(files,c(profil,pepfil,ddafil))
  if(!is.na(mvfil[1])){
    mvfils<- paste0("'",inputpath,mvfil,"'")%>% paste0(.,collapse=" ")
    system(paste0("cp -r ",mvfils," ",savepath))
  }
}
