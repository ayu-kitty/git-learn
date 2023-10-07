#!/opt/conda/bin/Rscript

#富集数据库类型
#' @export
getenrichtype<-function(type=tyte){
  daty<-data.frame("type"=c("G","K","I","R","W"),
                   "filn"=c("GO","KEGG","InterPro","Reactome","Wikipathways"),
                   "backn"=c("go","kegg","interpro","reactome","wikipathway")
  )
  rownames(daty)<-daty[,1]
  filn<-daty[type,2]
  backn<-daty[type,3]
  res<-list(filn,backn)
  names(res)<-c("filn","backn")
  return(res)
}

#KEGG背景整理函数
#' @export
funback<-function(bad,level){
  baid<-strsplit(bad[,2],",")%>%unlist()%>%unique()
  sub("[0-9].*","",baid[1])%>%gsub("ko",.,level[,1])->level[,1]
  babl<-level[which(level[,1] %in% baid),]
  term_pn<-sapply(1:length(babl[,1]),function(i){
    sum(str_detect(bad[,2],babl[i,1]))
  })
  babl[,5]<-term_pn
  babl[,6]<-unique(bad[,1])%>%length()
  return(babl)
}
#wikipathway/reactome背景整理函数
#' @export
funbackwr<-function(bad){
  baid<-strsplit(bad[,2],",")%>%unlist()%>%unique()
  babl<-cbind(unlist(strsplit(bad[,2],",")),trimws(unlist(strsplit(bad[,3],"[|]"))))%>%as.data.frame()%>%unique()
  term_pn<-sapply(1:length(babl[,1]),function(i){
    sum(str_detect(bad[,2],babl[i,1]))
  })
  babl[,3]<-term_pn
  babl[,4]<-unique(bad[,1])%>%length()
  return(babl)
}
#diff背景整理函数
#' @export
diffterm<-function(proid,bad,babl){
  diffb<-bad[which(bad[,1] %in% proid[,1]),] #提取差异蛋白对应go背景
  if(nrow(diffb)!=0){
    diffid<-strsplit(diffb[,2],",")%>%unlist()%>%unique()    #差异背景GOID
    diffbl<-babl[which(babl[,1] %in% diffid),] #
    #生成差异ID对应term的蛋白及个数
    term_pnd<-sapply(1:length(diffbl[,1]),function(i){
      sum(str_detect(diffb[,2],diffbl[i,1]))
    })   #个数
    term_pro<-sapply(1:length(diffbl[,1]),function(i){
      paste(diffb[str_detect(diffb[,2],diffbl[i,1]),4],collapse=",")
    })   #蛋白
    colll<-ncol(diffbl)
    diffbl[,colll+1]<-term_pnd
    diffbl[,colll+2]<-nrow(diffb)
    diffbl[,colll+3]<-term_pro
  }else diffbl<-as.data.frame(matrix(ncol = 10,nrow = 0))
  
  return(diffbl)
}

#GO富集函数
#' @export
goenrich<-function(db,con,outd){
  type<-c("Total","Up","Down")
  gob<-read.delim(paste0(db,"gene_go.backgroud.xls"), header=F, sep="\t", quote="")
  gob[,"id"]<-gsub(":.*","",gob[,1])
  gob<-select(gob,c(4,2,3,1))
  golevel<-read.delim(paste0(db,"category.xls"),header=F, sep="\t", quote="")
  golevel<-select(golevel,c(1,ncol(golevel)))
  godd<-cbind(unlist(strsplit(gob[,2],",")),unlist(strsplit(gob[,3],"[|]"))) %>% .[!duplicated(.),] %>% as.data.frame()
  gos<-strsplit(gob[,2],",")%>%unlist()%>%table()%>%as.data.frame()
  cate<-merge(godd,golevel,by="V1",all.x=T) %>% merge(.,gos,by.x="V1",by.y=".",all.x=T)
  cate[,5]<-nrow(gob)
  outdir<-paste0(outd,"GO/")
  createdir(outdir)
  prol<-readdata(con,header = T)
  prol[,1]<-gsub(":.*","",prol[,1])
  prol<-as.data.frame(prol)
  tre<-gsub("-diff-.*","",gsub(".*/","",con))
  createdir(paste0(outdir,tre))
  if(ncol(prol)==3){
    for(j in 1:3){
      if(type[j]=="Total"){
        prolt<-prol
      }else prolt<-prol[prol[,3]==type[j],]
      prolt<-prolt[!duplicated(prolt[,1]),]
      if(nrow(prolt)!=0 & length(intersect(prolt[,1],gob[,1]))!=0){
        enrich<-select(diffterm(prolt,gob,cate),c(1,3,2,6,7,4,5,8))
        colnames(enrich)<-c("id","Category","Term","ListHits","ListTotal","PopHits","PopTotal","Substances")
        enrich["p-value"] <- phyper(enrich[,"ListHits"]-1, enrich[,"PopHits"],
                                   enrich[,"PopTotal"]-enrich[,"PopHits"], enrich[,"ListTotal"], lower.tail=F)
        enrich["q-value"] <- p.adjust(enrich[,"p-value"], method="fdr")
        enrich["Enrichment_score"] <-
          (enrich["ListHits"]*enrich["PopTotal"])/(enrich["ListTotal"]*enrich["PopHits"])
        enrich <- enrich[order(enrich["p-value"]), ]
        enrich<-enrich[c(1,2,3,4,5,6,7,9,10,11,8)]
        write.table(enrich, paste0(outdir,tre, "/enrichment-go-",tre,"-",type[j],".xls" ), sep="\t", row.names=F, quote=F)
      }else{
        if(type[j]=="Total"){
          return("stop")
        }
      }
    }
  }else{
    prolt<-prol[!duplicated(prol[,1]),]%>%as.data.frame()
    if(nrow(prolt)!=0 & length(intersect(prolt[,1],gob[,1]))!=0){
      enrich<-select(diffterm(prolt,gob,cate),c(1,3,2,6,7,4,5,8))
      colnames(enrich)<-c("id","Category","Term","ListHits","ListTotal","PopHits","PopTotal","Substances")
      enrich["p-value"] <- phyper(enrich[,"ListHits"]-1, enrich[,"PopHits"],
                                 enrich[,"PopTotal"]-enrich[,"PopHits"], enrich[,"ListTotal"], lower.tail=F)
      enrich["q-value"] <- p.adjust(enrich[,"p-value"], method="fdr")
      enrich["Enrichment_score"] <-
        (enrich["ListHits"]*enrich["PopTotal"])/(enrich["ListTotal"]*enrich["PopHits"])
      enrich <- enrich[order(enrich["p-value"]), ]
      enrich<-enrich[c(1,2,3,4,5,6,7,9,10,11,8)]
      write.table(enrich, paste0(outdir,tre, "/enrichment-go-",tre,"-Total.xls" ), sep="\t", row.names=F, quote=F)
    }else{
      return("stop")
    }
  }
  
  return("run")
}

#KEGG富集函数
#' @export
keggenrich<-function(db,con,outd){
  type<-c("Total","Up","Down")
  bad<-read.delim(paste0(db,"gene_kegg.backgroud.xls"), header=F, sep="\t", quote="")
  bad[,"id"]<-gsub(":.*","",bad[,1])
  bad<-select(bad,c("id",2,3,1))
  data("kegglevel")
  outdir<-paste0(outd,"KEGG/")
  createdir(outdir)
  babl<-funback(bad,level)
  tre<-gsub("-diff-.*","",gsub(".*/","",con))
  
  prol<-readdata(con,header = T)
  prol[,1]<-gsub(":.*","",prol[,1])
  prol<-as.data.frame(prol)
  createdir(paste0(outdir,tre))
  if(ncol(prol)==3){
    for(j in 1:3){
      if(type[j]=="Total"){
        prolt<-prol
      }else prolt<-prol[prol[,3]==type[j],]
      prolt<-prolt[!duplicated(prolt[,1]),]
      if(nrow(prolt)!=0 & length(intersect(prolt[,1],bad[,1]))!=0){
        enrich<-select(diffterm(prolt,bad,babl),c(1,2,3,4,7,8,5,6,9))
        colnames(enrich)<-c("id","Classification_level1","Classification_level2","Term","ListHits","ListTotal","PopHits","PopTotal","Substances")
        enrich["p-value"] <- phyper(enrich[,"ListHits"]-1, enrich[,"PopHits"],
                                   enrich[,"PopTotal"]-enrich[,"PopHits"], enrich[,"ListTotal"], lower.tail=F)
        enrich["q-value"] <- p.adjust(enrich[,"p-value"], method="fdr")
        enrich["Enrichment_score"] <-
          (enrich["ListHits"]*enrich["PopTotal"])/(enrich["ListTotal"]*enrich["PopHits"])
        enrich <- enrich[order(enrich["p-value"]), ]
        enrich<-enrich[c(1,2,3,4,5,6,7,8,10,11,12,9)]
        write.table(enrich, paste0(outdir,tre, "/enrichment-kegg-",tre,"-",type[j],".xls" ), sep="\t", row.names=F, quote=F)
      }else{
        if(type[j]=="Total"){
          return("stop")
        }
      }
    }
  }else{
    prolt<-prol[!duplicated(prol[,1]),]%>%as.data.frame()
    if(nrow(prolt)!=0 & length(intersect(prolt[,1],bad[,1]))!=0){
      enrich<-select(diffterm(prolt,bad,babl),c(1,2,3,4,7,8,5,6,9))
      colnames(enrich)<-c("id","Classification_level1","Classification_level2","Term","ListHits","ListTotal","PopHits","PopTotal","Substances")
      enrich["p-value"] <- phyper(enrich[,"ListHits"]-1, enrich[,"PopHits"],
                                 enrich[,"PopTotal"]-enrich[,"PopHits"], enrich[,"ListTotal"], lower.tail=F)
      enrich["q-value"] <- p.adjust(enrich[,"p-value"], method="fdr")
      enrich["Enrichment_score"] <-
        (enrich["ListHits"]*enrich["PopTotal"])/(enrich["ListTotal"]*enrich["PopHits"])
      enrich <- enrich[order(enrich["p-value"]), ]
      enrich<-enrich[c(1,2,3,4,5,6,7,8,10,11,12,9)]
      write.table(enrich, paste0(outdir,tre, "/enrichment-kegg-",tre,"-Total.xls" ), sep="\t", row.names=F, quote=F)
      
    }else{
      return("stop")
    }
    
  }
  return("run")
}

#wikipathways/reactome/interpro富集函数
#' @export
enrich<-function(db,con,enclass,outd){
  type<-c("Total","Up","Down")
  bad<-read.delim(paste0(db,"gene_",getenrichtype(enclass)$backn,".backgroud.xls"), header=F, sep="\t", quote="")
  bad[,"id"]<-gsub(":.*","",bad[,1])
  bad<-select(bad,c(4,2,3,1))
  babl<-funbackwr(bad)
  outdir<-paste0(outd,getenrichtype(enclass)$filn,"/")
  tre<-gsub("-diff-.*","",gsub(".*/","",con))
  createdir(outdir)
  prol<-readdata(con,header = T)
  prol[,1]<-gsub(":.*","",prol[,1])
  prol<-as.data.frame(prol)
  createdir(paste0(outdir,tre))
  if(ncol(prol)==3){
    for(j in 1:3){
      if(type[j]=="Total"){
        prolt<-prol
      }else prolt<-prol[prol[,3]==type[j],]
      prolt<-prolt[!duplicated(prolt[,1]),]
      if(nrow(prolt)!=0 & length(intersect(prolt[,1],bad[,1]))!=0){
        enrich<-select(diffterm(prolt,bad,babl),c(1,2,5,6,3,4,7))
        colnames(enrich)<-c("id","Term","ListHits","ListTotal","PopHits","PopTotal","Substances")
        enrich["p-value"] <- phyper(enrich[,"ListHits"]-1, enrich[,"PopHits"],
                                   enrich[,"PopTotal"]-enrich[,"PopHits"], enrich[,"ListTotal"], lower.tail=F)
        enrich["q-value"] <- p.adjust(enrich[,"p-value"], method="fdr")
        enrich["Enrichment_score"] <-
          (enrich["ListHits"]*enrich["PopTotal"])/(enrich["ListTotal"]*enrich["PopHits"])
        enrich <- enrich[order(enrich["p-value"]), ]
        enrich<-enrich[c(1,2,3,4,5,6,8,9,10,7)]
        if(enclass!="I"){
          enrich[,"Term"]<-paste0("(",enrich[,"id"],")",enrich[,"Term"])
        }
        write.table(enrich, paste0(outdir,tre, "/enrichment-",getenrichtype(enclass)$filn,"-",tre,"-",type[j],".xls" ), sep="\t", row.names=F, quote=F)
      }else{
        if(type[j]=="Total"){
          return("stop")
        }
      }
    }
  }else{
    prolt<-prol[!duplicated(prol[,1]),]%>%as.data.frame()
    if(nrow(prolt)!=0 & length(intersect(prolt[,1],bad[,1]))!=0){
      enrich<-select(diffterm(prolt,bad,babl),c(1,2,5,6,3,4,7))
      colnames(enrich)<-c("id","Term","ListHits","ListTotal","PopHits","PopTotal","Substances")
      enrich["p-value"] <- phyper(enrich[,"ListHits"]-1, enrich[,"PopHits"],
                                 enrich[,"PopTotal"]-enrich[,"PopHits"], enrich[,"ListTotal"], lower.tail=F)
      enrich["q-value"] <- p.adjust(enrich[,"p-value"], method="fdr")
      enrich["Enrichment_score"] <-
        (enrich["ListHits"]*enrich["PopTotal"])/(enrich["ListTotal"]*enrich["PopHits"])
      enrich <- enrich[order(enrich["p-value"]), ]
      enrich<-enrich[c(1,2,3,4,5,6,8,9,10,7)]
      if(enclass!="I"){
        enrich[,"Term"]<-paste0("(",enrich[,"id"],")",enrich[,"Term"])
      }
      write.table(enrich, paste0(outdir,tre, "/enrichment-",getenrichtype(enclass)$filn,"-",tre,"-Total.xls" ), sep="\t", row.names=F, quote=F)
      
    }else{
      return("stop")
    }
    
  }
  return("run")
}
