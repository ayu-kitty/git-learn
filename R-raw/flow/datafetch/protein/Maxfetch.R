#' 获取修饰类型函数
#'
#' @param colna 位点表列名 
#' @export
ptm<-function(colna = colna){
	if("Phospho (STY) site IDs"%in% colna){
		ptmn<-"Phospho (STY) site IDs"
	}else if("Acetyl (K) site IDs"%in% colna){
		ptmn<-"Acetyl (K) site IDs"
	}else if("GlyGly (K) site IDs"%in% colna){
		ptmn<-"GlyGly (K) site IDs"
	}else if("Deamidation 18O (N) site IDs"%in% colna){
                ptmn<-"Deamidation 18O (N) site IDs"
        }
	ptmn
}
#' 去除反库污染物
#' @export
delrc<-function(dat = dat,proname="Leading razor protein"){
  if(proname=="F"){
    rc<-setdiff(1:nrow(dat),c(grep("+",dat[,"Reverse"],fixed = T),grep("+",dat[,"Potential contaminant"],fixed = T)))
  }else{
    rc1<-setdiff(1:nrow(dat),c(grep("+",dat[,"Reverse"],fixed = T),grep("+",dat[,"Potential contaminant"],fixed = T)))
    rc2<-setdiff(1:nrow(dat),c(grep("CON__",dat[,proname]),grep("REV__",dat[,proname])))
    rc<-unique(rc1,rc2)
  }
  dat[rc,]
}

#' 处理Maxquant labelfree下机数据
#' @export
lbfetch<-function(inputpath="./projectdata/",savepath="./rawdata/"){
  pro<-read.delim2(paste0(inputpath,"proteinGroups.txt"),check.names = F)
  pro<-delrc(pro,proname="Majority protein IDs")
  pep<-read.delim2(paste0(inputpath,"peptides.txt"),check.names = F)
  pep<-delrc(pep)
  #pro
  procol<-c("Protein IDs","Majority protein IDs","Fasta headers","Peptides","Unique peptides","Sequence coverage [%]","Mol. weight [kDa]","Reverse","Potential contaminant")
  samplecolpro <-grep("LFQ intensity ",colnames(pro))#获取样本名列数
  colnames(pro)<-gsub("LFQ intensity *","",colnames(pro))  #替换样本名列名
  pro[,"Accession"]<- gsub(";.*","",pro$`Majority protein IDs`)
  pro[,samplecolpro]<-apply(pro[,samplecolpro], 2,as.numeric)
  pro[,"Sequence coverage [%]"]<-as.numeric(pro[,"Sequence coverage [%]"])
  pro[,"Mol. weight [kDa]"]<-as.numeric(pro[,"Mol. weight [kDa]"])
  pro[,"Qvalue]"]<-as.numeric(pro[,"Q-value"])
  pro[samplecolpro][pro[samplecolpro]==0]<-NA
  if(file.exists("./background/Annotation_Data.xlsx")){
    anno<-readxlsx("./background/Annotation_Data.xlsx")[,-1]
    names(anno)[1]<-"Accession"
    pro[,"Gene Name"]<-left_join(pro,anno,by="Accession")$Gene_Name
  }else pro[,"Gene Name"]<-NA
  prode <- select(pro,c("Accession","Gene Name",all_of(procol),all_of(samplecolpro)))
  
  savexlsx(prode,paste0(savepath,"Protein quantitation.xlsx"))
  #pep
  pepcol<-c("Sequence","Mass","Proteins","Leading razor protein","Charges","Score")
  samplecolpep <-grep("LFQ intensity ",colnames(pep))
  colnames(pep)<-gsub("LFQ intensity *","",colnames(pep))
  pep[,samplecolpep]<-apply(pep[,samplecolpep], 2,as.numeric)
  pep[samplecolpep][pep[samplecolpep]==0]<-NA
  pepde <- select(pep,c(all_of(pepcol),all_of(samplecolpep)))
  
  savexlsx(pepde,paste0(savepath,"Peptides.xlsx"))
  ppm<-as.data.frame(matrix(nrow = 1,ncol=3))
  colnames(ppm)<-c("Items","蛋白数(Protein Groups)","肽段数(Peptides)")
  ppm[1,1]<-"FDR小于0.01"
  ppm[1,2]<-nrow(prode)
  ppm[1,3]<-nrow(pepde)
  savexlsx(ppm,paste0(savepath,"Identified_number.xlsx"))
  files<-dir(inputpath)
  mvfil<-setdiff(files,c("proteinGroups.txt","peptides.txt"))
  if(!is.na(mvfil[1])){
    mvfils<-paste0("'",inputpath,mvfil,"'")%>% paste0(.,collapse=" ")
    system(paste0("cp -r ",mvfils," ",savepath))
  }
}
#' 处理Maxquant 磷酸化TMT下机数据
#' @export
ptfetch<-function(inputpath="./projectdata/",savepath="./rawdata/"){
  pro<-read.delim2(paste0(inputpath,"proteinGroups.txt"),check.names = F)
  pro<-delrc(pro,proname="Protein IDs")
  pep<-read.delim2(paste0(inputpath,"modificationSpecificPeptides.txt"),check.names = F)
  pep<-delrc(pep,proname="F")
  fil<-dir(inputpath,pattern="*Sites")
  site<-read.delim2(paste0(inputpath,fil),check.names = F)
  site<-delrc(site,proname = "Leading proteins")
  label<-read.xlsx(paste0(inputpath,"样品标记对照表.xlsx"))
  grep("delet",label[,1])->del
  if(length(del)>0){
    labeln<-label[-del,]
  }else labeln<-label 
  ptmn<-ptm(colnames(pro))
  ptmname<-gsub(" site IDs","",ptmn)
  #pro
  procol<-c("Protein IDs","Majority protein IDs","Fasta headers","Peptides","Unique peptides","Sequence coverage [%]","Mol. weight [kDa]","Reverse","Potential contaminant","Oxidation (M) site positions",paste0(ptmname," site positions"))
  pro<-pro[pro[,ptmn]!="",]
  gsub("Reporter intensity *","",colnames(pro)) %>% grep("^\\d",.)->samn
  colnames(pro)[samn]<-label[,1]
  pro[,label[,1]]<-apply(pro[,label[,1]], 2,as.numeric)
  pro[,"Sequence coverage [%]"]<-as.numeric(pro[,"Sequence coverage [%]"])
  pro[,"Mol. weight [kDa]"]<-as.numeric(pro[,"Mol. weight [kDa]"])
  prode<-select(pro,c(all_of(procol),all_of(labeln[,1])))
  prode[,all_of(labeln[,1])][prode[,all_of(labeln[,1])]==0]<-NA
  savexlsx(prode,paste0(savepath,"Protein quantitation.xlsx"))
  #pep
  pepcol<-c("Sequence","Modifications","Mass","Proteins","Oxidation (M)",ptmname,"Missed cleavages","Retention time","Charges","PEP","Delta score")
  pep<-pep[pep[,ptmn]!="",]
  gsub("Reporter intensity *","",colnames(pep)) %>% grep("^\\d*$",.)->sa
  colnames(pep)[sa]<-label[,1]
  pep[,label[,1]]<-apply(pep[,label[,1]], 2,as.numeric)
  pepde <- select(pep,c(all_of(pepcol),all_of(labeln[,1])))
  pepde[,all_of(labeln[,1])][pepde[,all_of(labeln[,1])]==0]<-NA
  savexlsx(pepde,paste0(savepath,"Peptides.xlsx"))
  #site 标记需用三价之和
  sitecol<-c("Proteins","Positions within proteins","Leading proteins","Protein","Fasta headers","Localization prob","PEP","Delta score","Amino acid","Sequence window","Peptide window coverage",paste0(ptmname," Probabilities"),	"Position in peptide","Charge","Mass error [ppm]","Reverse","Potential contaminant")
  samplecolsite<-str_extract(colnames(site),"Reporter intensity [0-9]*___[0-9]*") %>% na.omit() %>% as.character()
  site[,samplecolsite]<-apply(site[,samplecolsite], 2,as.numeric)
  tempsite<-site[,samplecolsite]
  ssam<-gsub("___[0-9]*","",samplecolsite)
  names(tempsite)<-ssam
  site[,label[,1]]<-sapply(split.default(tempsite, factor(names(tempsite),levels=unique(names(tempsite)))), rowSums)
  site[,"Accession"]<- gsub(";.*","",site$Proteins)
  site[,"Positions within proteins"]<- gsub(";.*","",site$`Positions within proteins`) %>% as.numeric()
  site[,"ID"]<- paste0(site$Accession,":",site$`Amino acid`,site$`Positions within proteins`)
  site[,"Localization prob"]<-as.numeric(site[,"Localization prob"])
  site[,"Delta score"]<-as.numeric(site[,"Delta score"])
  site[,"PEP"]<-as.numeric(site[,"PEP"])
  if(file.exists("./background/Annotation_Data.xlsx")){
    anno<-readxlsx("./background/Annotation_Data.xlsx")[,-1]
    names(anno)[1]<-"Accession"
    site[,"Gene Name"]<-left_join(site,anno,by="Accession")$Gene_Name
  }else site[,"Gene Name"]<-NA
  newsite <- select(site,c("ID","Accession","Gene Name",all_of(sitecol),all_of(labeln[,1])))
  
  newsite[,all_of(labeln[,1])][newsite[,all_of(labeln[,1])]==0]<-NA
  savexlsx(newsite,paste0(savepath,gsub(".txt",".xlsx",fil)))
  ppm<-as.data.frame(matrix(nrow = 1,ncol=4))
  colnames(ppm)<-c("Items","蛋白数(Protein Groups)","肽段数(Peptides)","位点数(Sites)")
  ppm[1,1]<-"FDR小于0.01"
  ppm[1,2]<-nrow(prode)
  ppm[1,3]<-nrow(pepde)
  ppm[1,4]<-nrow(newsite)
  savexlsx(ppm,paste0(savepath,"Identified_number.xlsx"))
  files<-dir(inputpath)
  mvfil<-setdiff(files,c("proteinGroups.txt","modificationSpecificPeptides.txt",fil))
  if(!is.na(mvfil[1])){
    mvfils<-paste0("'",inputpath,mvfil,"'")%>% paste0(.,collapse=" ")
    system(paste0("cp -r ",mvfils," ",savepath))
  }
}
#' 处理Maxquant 修饰labelfree下机数据
#' @export
plfetch<-function(inputpath="./projectdata/",savepath="./rawdata/"){
  pro<-read.delim2(paste0(inputpath,"proteinGroups.txt"),check.names = F)
  pro<-delrc(pro,proname="Protein IDs")
  pep<-read.delim2(paste0(inputpath,"modificationSpecificPeptides.txt"),check.names = F)
  pep<-delrc(pep,proname="F")
  fil<-dir(inputpath,pattern="*Sites")
  site<-read.delim2(paste0(inputpath,fil),check.names = F)
  site<-delrc(site,proname = "Leading proteins")
  ptmn<-ptm(colnames(pro))
  ptmname<-gsub(" site IDs","",ptmn)
  #pro
  procol<-c("Protein IDs","Majority protein IDs","Fasta headers","Peptides","Unique peptides","Sequence coverage [%]","Mol. weight [kDa]","Reverse","Potential contaminant","Oxidation (M) site positions",paste0(ptmname," site positions"))
  samplecolpro <-grep("LFQ intensity ",colnames(pro))#获取样本名列数
  pro[,samplecolpro]<-apply(pro[,samplecolpro], 2,as.numeric)
  pro[,"Sequence coverage [%]"]<-as.numeric(pro[,"Sequence coverage [%]"])
  pro[,"Mol. weight [kDa]"]<-as.numeric(pro[,"Mol. weight [kDa]"])
  colnames(pro)<-gsub("LFQ intensity *","",colnames(pro))
  pro<-pro[pro[,ptmn]!="",]
  pro[samplecolpro][pro[samplecolpro]==0]<-NA
  prode<-select(pro,c(all_of(procol),all_of(samplecolpro)))
  savexlsx(prode,paste0(savepath,"Protein quantitation.xlsx"))
  #pep
  pepcol<-c("Sequence","Modifications","Mass","Proteins","Oxidation (M)",ptmname,"Missed cleavages","Retention time","Charges","PEP","Delta score")
  samplecolpep <-grep("Intensity ",colnames(pep))
  colnames(pep)<-gsub("Intensity *","",colnames(pep))
  pep<-pep[pep[,ptmn]!="",]
  pep[,samplecolpep]<-apply(pep[,samplecolpep], 2,as.numeric)
  pep[samplecolpep][pep[samplecolpep]==0]<-NA
  pepde <- select(pep,c(all_of(pepcol),all_of(samplecolpep)))
  
  savexlsx(pepde,paste0(savepath,"Peptides.xlsx"))
  #site
  sitecol<-c("Proteins","Positions within proteins","Leading proteins","Protein","Fasta headers","Localization prob","PEP","Delta score","Amino acid","Sequence window","Peptide window coverage",paste0(ptmname," Probabilities"),	"Position in peptide","Charge","Mass error [ppm]","Reverse","Potential contaminant")
  samplecolsite <-grep("Intensity ",colnames(site))
  site[,"Accession"]<- gsub(";.*","",site$Proteins)
  site[,"Localization prob"]<-as.numeric(site[,"Localization prob"])
  site[,"Positions within proteins"]<- gsub(";.*","",site$`Positions within proteins`) %>% as.numeric()
  site[,"ID"]<- paste0(site$Accession,":",site$`Amino acid`,site$`Positions within proteins`)
  site[,"Delta score"]<-as.numeric(site[,"Delta score"])
  site[,"PEP"]<-as.numeric(site[,"PEP"])
  site[,samplecolsite]<-apply(site[,samplecolsite], 2,as.numeric)
  setdiff(grep("Intensity ",colnames(site)),grep("___",colnames(site)))->sas
  colnames(site)[sas]<-gsub("Intensity *","",colnames(site)[sas])
  if(file.exists("./background/Annotation_Data.xlsx")){
    anno<-readxlsx("./background/Annotation_Data.xlsx")[,-1]
    names(anno)[1]<-"Accession"
    site[,"Gene Name"]<-left_join(site,anno,by="Accession")$Gene_Name
  }else site[,"Gene Name"]<-NA
  newsite <- select(site,c("ID","Accession","Gene Name",all_of(sitecol),all_of(colnames(site)[sas])))
  newsite[,all_of(colnames(site)[sas])][newsite[,all_of(colnames(site)[sas])]==0]<-NA
  savexlsx(newsite,paste0(savepath,gsub(".txt",".xlsx",fil)))
  ppm<-as.data.frame(matrix(nrow = 1,ncol=4))
  colnames(ppm)<-c("Items","蛋白数(Protein Groups)","肽段数(Peptides)","位点数(Sites)")
  ppm[1,1]<-"FDR小于0.01"
  ppm[1,2]<-nrow(prode)
  ppm[1,3]<-nrow(pepde)
  ppm[1,4]<-nrow(newsite)
  savexlsx(ppm,paste0(savepath,"Identified_number.xlsx"))
  files<-dir(inputpath)
  mvfil<-setdiff(files,c("proteinGroups.txt","modificationSpecificPeptides.txt",fil))
  if(!is.na(mvfil[1])){
    mvfils<-paste0("'",inputpath,mvfil,"'")%>% paste0(.,collapse=" ")
    system(paste0("cp -r ",mvfils," ",savepath))
  }
}