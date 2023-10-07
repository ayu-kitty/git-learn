
#' 空代与空转联合分析
#'
#' @param samplename 样本名
#' @param mode 正负离子模式
#' @param datapath 数据路径
#' @param keggpath kegg路径
#' @param species 物种
#' @param savapath 保存路径
#' @param meta 逻辑，是否仅对注释代谢物进行分析
#' @param missvalue 缺失值筛选标准
#' @param updowncolor 上下调颜色
#' @param prange 相关p值筛选标准
#' @param corrange 相关性筛选标准
#'
#' @export
TransAndMetaboanalyst <- function(samplename,
                                  mode = "neg",
                                  datapath = "./sample/union/uniondata/",
                                  keggpath = meta::databasepath(path = "/kegg/kegg2/"),
                                  species = "mmu",
                                  savapath = "./sample/union/result/",
                                  meta = T,
                                  missvalue = 0.5,
                                  updowncolor = c("red","blue"),
                                  prange = 0.05,
                                  corrange = 0.5){
  
  library("pathview")
  library("meta")
  library("dplyr")
  wd <- getwd()
  
  sample <- read.table(file = paste(datapath,samplename,mode,"samplename.txt",sep = "/"),
                       header = T,check.names = F,stringsAsFactors = F)
  sample <- sample[,c("transname","metaname")]
  
  
  metadata <- read.table(file = paste(datapath,samplename,mode,"metadata.txt",sep = "/"),
                         header = T,check.names = F,stringsAsFactors = F)
  if(meta){
    metadata <- metadata[!is.na(metadata$Metabolites),]
  }
  transdata <- read.table(file = paste(datapath,samplename,mode,"transdata.txt",sep = "/"),
                          header = T,check.names = F,stringsAsFactors = F)
  
  metadata <- metadata[apply(metadata[,colnames(metadata) %in% sample$metaname], 1, function(x){sum(x == 0)/length(x)}) < missvalue,]
  metadata <- metadata[!duplicated(metadata$Metabolites),]
  
  transdata <- transdata[apply(transdata[,colnames(transdata) %in% sample$transname], 1, function(x){sum(x == 0)/length(x)}) < missvalue,]
  transdata <- transdata[!duplicated(transdata$Accession),]
  
  # 相关性分析
  meta::setwdfile(paste0(savapath,"/1.correlation/",samplename),force = T)
  cormetadata <- metadata[,sample$metaname]
  row.names(cormetadata) <- metadata$mz
  cortransdata <- transdata[,sample$transname]
  row.names(cortransdata) <- transdata$Accession
  
  cor <- meta::auto_cornetworkto(data = cormetadata,
                                 y = t(cortransdata),
                                 name = "cornetwork",
                                 prange = prange,
                                 corrange = corrange,
                                 yname = "Trans",
                                 nummax = 10000,
                                 dealname = F,
                                 allana = F,num=3)
  
  metadata <- metadata[metadata$mz %in% cor$data$name,]
  transdata <- transdata[transdata$Accession %in% cor$data$linkname,]
  
  cormetadata <- metadata[,sample$metaname]
  row.names(cormetadata) <- metadata$mz
  cortransdata <- transdata[,sample$transname]
  row.names(cortransdata) <- transdata$Accession
  
  cor <- auto_correlation(data = cormetadata,
                          y = t(cortransdata),
                          name = "T_M",
                          insig = "label_sig",
                          sig.level = c(0.001,0.01,0.05),
                          pch.col = "black",
                          order = "original",
                          nummax = 10000,
                          fnumber=2,
                          height = ifelse(dim(cormetadata)[1]>9,dim(cormetadata)[1]/3,6),
                          width = ifelse(dim(cortransdata)[1]>9,dim(cortransdata)[1]/3,6))
  
  colnames(cormetadata) <- colnames(cortransdata)
  allcordata <- rbind(cormetadata,cortransdata)
  cor <- auto_correlation(data = allcordata,
                          name = "TM",
                          insig = "label_sig",
                          sig.level = c(0.001,0.01,0.05),
                          pch.col = "black",
                          order = "hclust",
                          mode = "full",
                          fnumber=2,
                          height = ifelse(dim(allcordata)[1]>9,dim(allcoradata)[1]/3,6),
                          width = ifelse(dim(allcordata)[1]>9,dim(allcoradata)[1]/3,6))
  
  
  savexlsx5(data = cortransdata,name ="数据矩阵.xlsx",sheet = "转录")
  savexlsx5(data = cormetadata,name ="数据矩阵.xlsx",sheet = "代谢")
  
  setwd(wd)
  
  # 数据处理
  if (sum(!is.na(transdata$Accession))>1 & sum(!is.na(metadata$KEGG))>1){
    geneiddata <- read.table(file = paste0(keggpath,"/",species,"/data/gene"),
                             check.names = F,
                             stringsAsFactors = F,
                             header = T,
                             sep = "\t")
    geneiddata <- geneiddata[!is.na(geneiddata$`gene name`),]
    geneiddata <- geneiddata[!is.na(geneiddata$`gene id`),]
    geneiddata <- geneiddata[!duplicated(geneiddata$`gene id`),]
    geneiddata <- geneiddata[geneiddata$`gene name` %in% transdata$`Gene Name`,]
    colnames(geneiddata)[1] <- "Gene ID"
    transdata2 <- transdata[,c("Accession","Gene Name")]
    transdata2[,"log2(FC)"] <- 0
    transdata2 <- merge(x = transdata2,y = geneiddata[,c(1,2)],by.x = "Gene Name",by.y = "gene name",all = F)
    transdata2[,"up/down"] <- "up"
    transdata2[transdata2[,"log2(FC)"]  < 0,"up/down"] <- "down"
    transdata2[transdata2[,"log2(FC)"]  == 0,"up/down"] <- "none"
    
    keggdata <- metadata[,c("Metabolites","KEGG")]
    keggdata[,"log2(FC)"] <- 0
    keggdata <- keggdata[!is.na(keggdata[,"KEGG"]),]
    keggdata <- keggdata[!duplicated(keggdata[,"KEGG"]),]
    keggdata[,"up/down"] <- "up"
    keggdata[keggdata[,"log2(FC)"]  < 0,"up/down"] <- "down"
    keggdata[keggdata[,"log2(FC)"]  == 0,"up/down"] <- "none"
    
    # 通路信息整理
    pathwaydata <- read.table(file = paste0(keggpath,"/",species,"/data/entry"),
                              check.names = F,
                              stringsAsFactors = F,
                              header = T,
                              sep = "\t")
    pathwaydata_gene <- pathwaydata[pathwaydata$name %in% transdata2$`Gene ID`,]
    pathwaydata_cpd <- pathwaydata[pathwaydata$name %in% keggdata[,"KEGG"],]
    
    pathwaydata_gene1 <- pathwaydata_gene[!duplicated(pathwaydata_gene$pathwayname),]
    pathwaydata_cpd1 <- pathwaydata_cpd[!duplicated(pathwaydata_cpd$pathwayname),]
    
    movepathway <- intersect(pathwaydata_gene1$pathwayname,pathwaydata_cpd1$pathwayname)
    pathwaydata_gene <- pathwaydata_gene[pathwaydata_gene$pathwayname %in% movepathway,]
    pathwaydata_cpd <- pathwaydata_cpd[pathwaydata_cpd$pathwayname %in% movepathway,]
    pathwaydata_all <- rbind(pathwaydata_cpd,pathwaydata_gene)
    pathwaydata_all2 <- pathwaydata[pathwaydata$pathwayname %in% movepathway,]
    pathwaydata_all2 <- pathwaydata_all2[pathwaydata_all2$name %in% movepathway,]
    pathwaydata_all <- rbind(pathwaydata_all,pathwaydata_all2)
    
    meta::setwdfile(paste0(savapath,"/2.pathway_map/",samplename),force = T)
    if (length(movepathway)>0){
      # keggmap转移
      file.copy(from = paste0(keggpath,"/",species,"/",movepathway,".html"),
                to = "./",
                overwrite = T)
      file.copy(from = paste0(keggpath,"/",species,"/",movepathway,".xml"),
                to = "./",
                overwrite = T)
      file.copy(from = paste0(keggpath,"/",species,"/",movepathway,".png"),
                to = "./",
                overwrite = T)
      
      # keggmap上色
      gene.data <- transdata2$`log2(FC)`
      names(gene.data) <- transdata2$`Gene ID`
      meta.data <- keggdata$`log2(FC)`
      names(meta.data) <- keggdata[,"KEGG"]
      for ( i in 1:length(movepathway)) {
        pathview(gene.data = gene.data,
                 cpd.data = meta.data,
                 gene.idtype = "kegg",
                 pathway.id = movepathway[i],
                 species = species,
                 kegg.dir = "./",
                 keys.align = "y",
                 kegg.native = T,
                 key.pos = "topright",
                 same.layer=T,
                 # new.signature=F,
                 match.data=T,
                 limit = list(gene = c(-1,1), cpd = c(-1,1)),
                 bins = list(gene = 2, cpd = 2),
                 discrete = list(gene = F, cpd = F),
                 low = list(gene = updowncolor[2], cpd = updowncolor[2]),
                 mid = list(gene = updowncolor[1], cpd = updowncolor[1]),
                 high = list(gene = updowncolor[1], cpd = updowncolor[1]))
      }
      
      # 移除xml及修改png
      unlink(paste0(movepathway,".xml"),force = T)
      unlink(paste0(movepathway,".png"),force = T)
      file.rename(from = paste0(movepathway,".pathview.png"),to = paste0(movepathway,".png"))
      
      # 数据处理后保存
      
      pathwayinfo <- read.table(file = paste0(keggpath,"/",species,"/data/info"),
                                check.names = F,
                                stringsAsFactors = F,
                                header = T,
                                sep = "\t")
      pathwayinfo <- pathwayinfo[,c("name","title","level1","level2")]
      pathwaydata_gene <- pathwaydata_gene[,c("name","pathwayname")]
      pathwaydata_gene <- merge(x = transdata2,y = pathwaydata_gene,by.x = "Gene ID",by.y = "name")
      pathwaydata_gene <- merge(x = pathwaydata_gene,y = pathwayinfo,by.x = "pathwayname",by.y = "name")
      names(pathwaydata_gene)[1] <- "Pathway"
      
      pathwaydata_cpd <- pathwaydata_cpd[,c("name","pathwayname")]
      pathwaydata_cpd <- merge(x = keggdata,y = pathwaydata_cpd,by.x = "KEGG",by.y = "name")
      pathwaydata_cpd <- merge(x = pathwaydata_cpd,y = pathwayinfo,by.x = "pathwayname",by.y = "name")
      names(pathwaydata_cpd)[1] <- "Pathway"
      
      savexlsx1(data = transdata2,name ="ID.xlsx",sheet = "转录")
      savexlsx1(data = keggdata,name ="ID.xlsx",sheet = "代谢")
      savexlsx1(data = pathwaydata_gene,name ="Pathway.xlsx",sheet = "转录")
      savexlsx1(data = pathwaydata_cpd,name ="Pathway.xlsx",sheet = "代谢")
      
      setwd(wd)
      
      # 3.KGML分析
      meta::setwdfile(paste0(savapath,"/3.kgml/",samplename),force = T)
      
      # 节点数据整理
      nodedatatrans <- transdata2[,c("Gene ID","Accession","log2(FC)","up/down")]
      colnames(nodedatatrans)[1:2] <- c("Node","Name")
      nodedatatrans[,"type"] <- "gene"
      nodedatametabo <- keggdata[,c("KEGG","Metabolites","log2(FC)","up/down")]
      colnames(nodedatametabo)[1:2] <- c("Node","Name")
      nodedatametabo[,"type"] <- "compound"
      nodedata <- rbind(nodedatatrans,nodedatametabo)
      nodedatapathway <- pathwayinfo[pathwayinfo$name %in% movepathway,c("name","title")]
      colnames(nodedatapathway)[1:2] <- c("Node","Name")
      nodedatapathway[,"log2(FC)"] <- 0
      nodedatapathway[,"up/down"] <- "map"
      # nodedatapathway[,"up/down"] <- whto(a = pathwayinfo[,c("name","level1")],b = nodedatapathway$Node)
      nodedatapathway[,"type"] <- "map"
      nodedata <- rbind(nodedata,nodedatapathway)
      
      # 链接数据整理
      relationdata <- read.table(file = paste0(keggpath,"/",species,"/data/relation"),
                                 check.names = F,
                                 stringsAsFactors = F,
                                 header = T,
                                 sep = "\t")
      entry1 <- pathwaydata_all[,c("id","pathwayname","name")]
      colnames(entry1) <- c("entry1","pathwayname","from")
      relationdata <- merge(x = relationdata,y = entry1,by = c("entry1","pathwayname"))
      
      entry2 <- pathwaydata_all[,c("id","pathwayname","name")]
      colnames(entry2) <- c("entry2","pathwayname","to")
      relationdata <- merge(x = relationdata,y = entry2,by = c("entry2","pathwayname"))
      
      relationdata <- relationdata[,c("from","to","type")]
      relationdata2 <- pathwaydata_all[,c("name","pathwayname","type")]
      colnames(relationdata2) <- c("from","to","type")
      relationdata2[,"type"] <- "maplink"
      relationdata <- rbind(relationdata,relationdata2)
      relationdata <- relationdata[!duplicated(relationdata[,c("from","to")]),]
      relationdata <- relationdata[relationdata$from != relationdata$to,]
      
      relationdata2 <- relationdata
      relationdata2[,"ID"] <- 1:dim(relationdata2)[1]
      relationdata2[,"group"] <- "up"
      relationdata3 <- relationdata2[,c("to","from","type","ID")]
      colnames(relationdata3) <- c("from","to","type","ID")
      relationdata3[,"group"] <- "down"
      relationdata4 <- rbind(relationdata2,relationdata3)
      relationdata4 <- relationdata4[!duplicated(relationdata4[,c("from","to")]),]
      relationdata4 <- relationdata4[relationdata4$group == "down",]
      relationdata2 <- relationdata2[relationdata2$ID %in% relationdata4$ID,]
      relationdata <- relationdata2[,c("from","to","type")]
      relationdata <- relationdata[relationdata$from != relationdata$to,]
      
      for ( i in 1:dim(nodedata)[1]) {
        nodedata[i,"Degree"] <- dim(relationdata[(relationdata[,"from"] %in% nodedata$Node[i]) | (relationdata[,"to"] %in% nodedata$Node[i]),])[1]
      }
      nodedata <- nodedata[nodedata$Degree !=0,]
      nodedata <- nodedata[!duplicated(nodedata[,1]),]
      
      networkdata <- list(relationdata = relationdata,
                          nodedata = nodedata)
      
      allkeggnetwork(data = networkdata,
                     shape = c("compound" = 21,
                               "gene" = 24,
                               "map" = 22),
                     shapelabels = c("compound","gene","map"),
                     width = 10,height = 10)
      
      savexlsx1(data = nodedata,name ="Network.xlsx",sheet = "node")
      savexlsx1(data = relationdata,name ="Network.xlsx",sheet = "link")
      
      setwd(wd)
    }else{
      message(paste0("由于样本：", samplename, " 空代和空转没有共同pathway，所以没有KEGG富集分析、KEGG网络图分析。"))
      write(paste0("由于样本：", samplename, " 空代和空转没有共同pathway，所以没有KEGG富集分析、KEGG网络图分析。"), "说明.txt", append = T)
      setwd(wd)
    }
  }else{
    message(paste0("由于样本：", samplename," 空代和空转没有共同pathway，所以没有KEGG富集分析、KEGG网络图分析。"))
    write(paste0("由于样本：", samplename," 空代和空转没有共同pathway，所以没有KEGG富集分析、KEGG网络图分析。"), "说明.txt", append = T)
  }
}

#' 空代与空转联合分析
#'
#' @param ZoneFile 含有聚类或微区标签的barcode或坐标（x，y）表格的路径
#' @param Collable 标签列名
#' @param subset 聚类或微区中的子集，默认为NULL，即所有子区域
#' @param cluster 趋势分析时，自定义簇的个数，默认为NULL
#' @param num 满足相关性筛选的最少feature个数
#' 
#' @export
TransAndMetaboanalyst1 <-  function(samplename,
                                    mode = "neg",
                                    datapath = "./sample/union/uniondata/",
                                    keggpath = meta::databasepath(path = "/kegg/kegg2/"),
                                    species = "mmu",
                                    savapath = "./sample/union/result/",
                                    meta = T,
                                    ZoneFile = "./sample/union/transdata/",
                                    subset = NULL,
                                    Collable="clusters",
                                    cluster = NULL,
                                    missvalue = 0.5,
                                    updowncolor = c("red","blue"),
                                    prange = 0.05, num = 6,
                                    corrange = 0.5){
  
  library("pathview")
  library("meta")
  wd <- getwd()
  
  # data
  sample <- read.table(file = paste(datapath,samplename,mode,"samplename.txt",sep = "/"), header = T,check.names = F,stringsAsFactors = F)
  sample <- sample[,c("transname","metaname")]
  metadata0 <- read.table(file = paste(datapath,samplename,mode,"metadata.txt",sep = "/"), header = T,check.names = F,stringsAsFactors = F)
  if(meta){metadata0 <- metadata0[!is.na(metadata0$Metabolites),]}
  transdata0 <- read.table(file = paste(datapath,samplename,mode,"transdata.txt",sep = "/"), header = T,check.names = F,stringsAsFactors = F)
  
  # barcode of subset
  subset_barcode <- SelectZone(ZoneFile=ZoneFile, Select=subset, SampleName = samplename, Collable=Collable)
  for (sb in names(subset_barcode)){
    base::print(sb)
    transname_subset <- base::intersect(subset_barcode[[sb]], colnames(transdata0))
    metaname_subset <- sample[sample$transname %in% transname_subset,"metaname"]
    if (length(metaname_subset)>1){
      
      metadata <- metadata0[apply(metadata0[,colnames(metadata0) %in% metaname_subset], 1, function(x){sum(as.numeric(x) == 0, na.rm = T)/length(x)}) < missvalue,
                            c(colnames(metadata0)[1:17], metaname_subset)
      ]
      metadata <- metadata[!duplicated(metadata$Metabolites),]
      transdata <- transdata0[apply(transdata0[,colnames(transdata0) %in% subset_barcode[[sb]]], 1, function(x){sum(as.numeric(x) == 0, na.rm = T)/length(x)}) < missvalue,
                              c(colnames(transdata0)[1:2], transname_subset)
      ]
      transdata <- transdata[!duplicated(transdata$Accession),]
      
      cormetadata <- metadata[, metaname_subset] 
      row.names(cormetadata) <- metadata$Metabolites
      cortransdata <- transdata[, transname_subset] 
      row.names(cortransdata) <- transdata$Accession
      
      if (dim(cormetadata)[2]>10&dim(cortransdata)[2]>10){
        
        # 趋势分析
        meta::setwdfile(file.path(savapath,"1.trendAnalysis",samplename, sb), force = T)
        cormetadata_ <- cormetadata
        colnames(cormetadata_) <- colnames(cortransdata)
        trenddata <- rbind(cortransdata, cormetadata_)
        
        cluster1 <- cluster %||% min(as.integer(nrow(trenddata)**0.5), 9)
        auto_stemanalyst1(data=trenddata, cluster = cluster1)
        setwd(wd)
        
        # 相关性分析
        meta::setwdfile(file.path(savapath,"2.correlation",samplename, sb), force = T)
        cor <- meta::auto_cornetworkto(data = cormetadata,
                                       y = t(cortransdata),
                                       name = "cornetwork",
                                       prange = prange,
                                       corrange = corrange,
                                       yname = "Trans",
                                       nummax = 10000,
                                       num = num,
                                       dealname = F,
                                       allana = F)
        
        if (!is.null(cor$data)){
          
          # filter by cor and pvalue
          metadata <- metadata[metadata$Metabolites %in% cor$data$name,]
          transdata <- transdata[transdata$Accession %in% cor$data$linkname,]
          
          cormetadata <- metadata[, metaname_subset]
          row.names(cormetadata) <- metadata$Metabolites
          cortransdata <- transdata[, transname_subset]
          row.names(cortransdata) <- transdata$Accession
          
          cor <- auto_correlation(data = cormetadata,
                                  y = t(cortransdata),
                                  name = "T_M",
                                  insig = "label_sig",
                                  sig.level = c(0.001,0.01,0.05),
                                  pch.col = "black",
                                  order = "original",
                                  nummax = 10000,
                                  height = min(dim(cormetadata)[1]/5, 80),
                                  width = min(dim(cortransdata)[1]/5, 80))
          
          colnames(cormetadata) <- colnames(cortransdata)
          allcordata <- rbind(cormetadata,cortransdata)
          cor <- auto_correlation(data = allcordata,
                                  name = "TM",
                                  insig = "label_sig",
                                  sig.level = c(0.001,0.01,0.05),
                                  pch.col = "black",
                                  order = "hclust",
                                  mode = "full",
                                  height = min(dim(allcordata)[1]/10, 80),
                                  width = min(dim(allcordata)[1]/10, 80))
          
          
          meta::savexlsx5(data = cortransdata,name ="数据矩阵.xlsx",sheet = "转录")
          meta::savexlsx5(data = cormetadata,name ="数据矩阵.xlsx",sheet = "代谢")
          
          setwd(wd)
          
          # 数据处理
          if (sum(!is.na(transdata$Accession))>1 & sum(!is.na(metadata$KEGG))>1){
            geneiddata <- read.table(file = paste0(keggpath,"/",species,"/data/gene"),
                                     check.names = F,
                                     stringsAsFactors = F,
                                     header = T,
                                     sep = "\t")
            geneiddata <- geneiddata[!is.na(geneiddata$`gene name`),]
            geneiddata <- geneiddata[!is.na(geneiddata$`gene id`),]
            geneiddata <- geneiddata[!duplicated(geneiddata$`gene id`),]
            geneiddata <- geneiddata[geneiddata$`gene name` %in% transdata$`Gene Name`,]
            colnames(geneiddata)[1] <- "Gene ID"
            transdata2 <- transdata[,c("Accession","Gene Name")]
            transdata2[,"log2(FC)"] <- 0
            transdata2 <- merge(x = transdata2,y = geneiddata[,c(1,2)],by.x = "Gene Name",by.y = "gene name",all = F)
            transdata2[,"up/down"] <- "up"
            transdata2[transdata2[,"log2(FC)"]  < 0,"up/down"] <- "down"
            transdata2[transdata2[,"log2(FC)"]  == 0,"up/down"] <- "none"
            
            keggdata <- metadata[,c("Metabolites","KEGG")]
            keggdata[,"log2(FC)"] <- 0
            keggdata <- keggdata[!is.na(keggdata[,"KEGG"]),]
            keggdata <- keggdata[!duplicated(keggdata[,"KEGG"]),]
            keggdata[,"up/down"] <- "up"
            keggdata[keggdata[,"log2(FC)"]  < 0,"up/down"] <- "down"
            keggdata[keggdata[,"log2(FC)"]  == 0,"up/down"] <- "none"
            
            # 通路信息整理
            pathwaydata <- read.table(file = paste0(keggpath,"/",species,"/data/entry"),
                                      check.names = F,
                                      stringsAsFactors = F,
                                      header = T,
                                      sep = "\t")
            pathwaydata_gene <- pathwaydata[pathwaydata$name %in% transdata2$`Gene ID`,]
            pathwaydata_cpd <- pathwaydata[pathwaydata$name %in% keggdata[,"KEGG"],]
            
            pathwaydata_gene1 <- pathwaydata_gene[!duplicated(pathwaydata_gene$pathwayname),]
            pathwaydata_cpd1 <- pathwaydata_cpd[!duplicated(pathwaydata_cpd$pathwayname),]
            
            movepathway <- intersect(pathwaydata_gene1$pathwayname,pathwaydata_cpd1$pathwayname)
            pathwaydata_gene <- pathwaydata_gene[pathwaydata_gene$pathwayname %in% movepathway,]
            pathwaydata_cpd <- pathwaydata_cpd[pathwaydata_cpd$pathwayname %in% movepathway,]
            pathwaydata_all <- rbind(pathwaydata_cpd,pathwaydata_gene)
            pathwaydata_all2 <- pathwaydata[pathwaydata$pathwayname %in% movepathway,]
            pathwaydata_all2 <- pathwaydata_all2[pathwaydata_all2$name %in% movepathway,]
            pathwaydata_all <- rbind(pathwaydata_all,pathwaydata_all2)
            
            meta::setwdfile(file.path(savapath,"3.pathway_map",samplename, sb),force = T)
            if (length(movepathway)>0){
              # keggmap转移
              file.copy(from = paste0(keggpath,"/",species,"/",movepathway,".html"),
                        to = "./",
                        overwrite = T)
              file.copy(from = paste0(keggpath,"/",species,"/",movepathway,".xml"),
                        to = "./",
                        overwrite = T)
              file.copy(from = paste0(keggpath,"/",species,"/",movepathway,".png"),
                        to = "./",
                        overwrite = T)
              
              # keggmap上色
              gene.data <- transdata2$`log2(FC)`
              names(gene.data) <- transdata2$`Gene ID`
              meta.data <- keggdata$`log2(FC)`
              names(meta.data) <- keggdata[,"KEGG"]
              for ( i in 1:length(movepathway)) {
                pathview(gene.data = gene.data,
                         cpd.data = meta.data,
                         gene.idtype = "kegg",
                         pathway.id = movepathway[i],
                         species = species,
                         kegg.dir = "./",
                         keys.align = "y",
                         kegg.native = T,
                         key.pos = "topright",
                         same.layer=T,
                         # new.signature=F,
                         match.data=T,
                         limit = list(gene = c(-1,1), cpd = c(-1,1)),
                         bins = list(gene = 2, cpd = 2),
                         discrete = list(gene = F, cpd = F),
                         low = list(gene = updowncolor[2], cpd = updowncolor[2]),
                         mid = list(gene = updowncolor[1], cpd = updowncolor[1]),
                         high = list(gene = updowncolor[1], cpd = updowncolor[1]))
              }
              
              # 移除xml及修改png
              unlink(paste0(movepathway,".xml"),force = T)
              unlink(paste0(movepathway,".png"),force = T)
              file.rename(from = paste0(movepathway,".pathview.png"),to = paste0(movepathway,".png"))
              
              # 数据处理后保存
              pathwayinfo <- read.table(file = paste0(keggpath,"/",species,"/data/info"),
                                        check.names = F,
                                        stringsAsFactors = F,
                                        header = T,
                                        sep = "\t")
              pathwayinfo <- pathwayinfo[,c("name","title","level1","level2")]
              pathwaydata_gene <- pathwaydata_gene[,c("name","pathwayname")]
              pathwaydata_gene <- merge(x = transdata2,y = pathwaydata_gene,by.x = "Gene ID",by.y = "name")
              pathwaydata_gene <- merge(x = pathwaydata_gene,y = pathwayinfo,by.x = "pathwayname",by.y = "name")
              names(pathwaydata_gene)[1] <- "Pathway"
              
              pathwaydata_cpd <- pathwaydata_cpd[,c("name","pathwayname")]
              pathwaydata_cpd <- merge(x = keggdata,y = pathwaydata_cpd,by.x = "KEGG",by.y = "name")
              pathwaydata_cpd <- merge(x = pathwaydata_cpd,y = pathwayinfo,by.x = "pathwayname",by.y = "name")
              names(pathwaydata_cpd)[1] <- "Pathway"
              
              meta::savexlsx1(data = transdata2,name ="ID.xlsx",sheet = "转录")
              meta::savexlsx1(data = keggdata,name ="ID.xlsx",sheet = "代谢")
              meta::savexlsx1(data = pathwaydata_gene,name ="Pathway.xlsx",sheet = "转录")
              meta::savexlsx1(data = pathwaydata_cpd,name ="Pathway.xlsx",sheet = "代谢")
              
              setwd(wd)
              
              # 3.KGML分析
              meta::setwdfile(file.path(savapath,"4.kgml",samplename, sb),force = T)
              
              # 节点数据整理
              nodedatatrans <- transdata2[,c("Gene ID","Accession","log2(FC)","up/down")]
              colnames(nodedatatrans)[1:2] <- c("Node","Name")
              nodedatatrans[,"type"] <- "gene"
              nodedatametabo <- keggdata[,c("KEGG","Metabolites","log2(FC)","up/down")]
              colnames(nodedatametabo)[1:2] <- c("Node","Name")
              nodedatametabo[,"type"] <- "compound"
              nodedata <- rbind(nodedatatrans,nodedatametabo)
              nodedatapathway <- pathwayinfo[pathwayinfo$name %in% movepathway,c("name","title")]
              colnames(nodedatapathway)[1:2] <- c("Node","Name")
              nodedatapathway[,"log2(FC)"] <- 0
              nodedatapathway[,"up/down"] <- "map"
              # nodedatapathway[,"up/down"] <- whto(a = pathwayinfo[,c("name","level1")],b = nodedatapathway$Node)
              nodedatapathway[,"type"] <- "map"
              nodedata <- rbind(nodedata,nodedatapathway)
              
              # 链接数据整理
              relationdata <- read.table(file = paste0(keggpath,"/",species,"/data/relation"),
                                         check.names = F,
                                         stringsAsFactors = F,
                                         header = T,
                                         sep = "\t")
              entry1 <- pathwaydata_all[,c("id","pathwayname","name")]
              colnames(entry1) <- c("entry1","pathwayname","from")
              relationdata <- merge(x = relationdata,y = entry1,by = c("entry1","pathwayname"))
              
              entry2 <- pathwaydata_all[,c("id","pathwayname","name")]
              colnames(entry2) <- c("entry2","pathwayname","to")
              relationdata <- merge(x = relationdata,y = entry2,by = c("entry2","pathwayname"))
              
              relationdata <- relationdata[,c("from","to","type")]
              relationdata2 <- pathwaydata_all[,c("name","pathwayname","type")]
              colnames(relationdata2) <- c("from","to","type")
              relationdata2[,"type"] <- "maplink"
              relationdata <- rbind(relationdata,relationdata2)
              relationdata <- relationdata[!duplicated(relationdata[,c("from","to")]),]
              relationdata <- relationdata[relationdata$from != relationdata$to,]
              
              relationdata2 <- relationdata
              relationdata2[,"ID"] <- 1:dim(relationdata2)[1]
              relationdata2[,"group"] <- "up"
              relationdata3 <- relationdata2[,c("to","from","type","ID")]
              colnames(relationdata3) <- c("from","to","type","ID")
              relationdata3[,"group"] <- "down"
              relationdata4 <- rbind(relationdata2,relationdata3)
              relationdata4 <- relationdata4[!duplicated(relationdata4[,c("from","to")]),]
              relationdata4 <- relationdata4[relationdata4$group == "down",]
              relationdata2 <- relationdata2[relationdata2$ID %in% relationdata4$ID,]
              relationdata <- relationdata2[,c("from","to","type")]
              relationdata <- relationdata[relationdata$from != relationdata$to,]
              
              for ( i in 1:dim(nodedata)[1]) {
                nodedata[i,"Degree"] <- dim(relationdata[(relationdata[,"from"] %in% nodedata$Node[i]) | (relationdata[,"to"] %in% nodedata$Node[i]),])[1]
              }
              nodedata <- nodedata[nodedata$Degree !=0,]
              nodedata <- nodedata[!duplicated(nodedata[,1]),]
              
              networkdata <- list(relationdata = relationdata,
                                  nodedata = nodedata)
              
              meta::allkeggnetwork(data = networkdata,
                                   shape = c("compound" = 21,
                                             "gene" = 24,
                                             "map" = 22),
                                   shapelabels = c("compound","gene","map"),
                                   width = 10,height = 10)
              
              meta::savexlsx1(data = nodedata,name ="Network.xlsx",sheet = "node")
              meta::savexlsx1(data = relationdata,name ="Network.xlsx",sheet = "link")
              
              setwd(wd)
              
            }else{
              message(paste0("由于样本：", samplename, "-",Collable,"-",sb," 空代和空转没有共同pathway，所以没有KEGG富集分析、KEGG网络图分析。"))
              write(paste0("由于样本：", samplename, "-",Collable,"-",sb," 空代和空转没有共同pathway，所以没有KEGG富集分析、KEGG网络图分析。"), "说明.txt", append = T)
              setwd(wd)
            }
          }else{
            message(paste0("由于样本：", samplename, "-",Collable,"-",sb," 空代和空转没有共同pathway，所以没有KEGG富集分析、KEGG网络图分析。"))
            write(paste0("由于样本：", samplename, "-",Collable,"-",sb," 空代和空转没有共同pathway，所以没有KEGG富集分析、KEGG网络图分析。"), "说明.txt", append = T)
          }
        }else{
          unlink("说明.txt", force = T, recursive = T)
          write(paste0("由于样本：", samplename, "-",Collable,"-",sb,"相关性分析中，相关性>",
                       corrange, ",显著性<", prange, "标准筛选后的关联少于",
                       num, "个, 所以没有相关性热图和网络图，也没有KEGG网络图分析。"), "说明.txt", append = T)
        }
        setwd(wd)
      }else{
        write(paste0("由于样本：", samplename, "-",Collable,"-",sb,", spot点数小于10, 所以没有空间趋势分析、相关性分析、KEGG网络图分析。"), "说明.txt", append = T)
      }
    }
  }
}

#' 选择区域
#' 
#' @param ZoneFile 含有聚类或微区标签的barcode或坐标（x，y）表格的路径
#' @param Collable 标签列名
#' @param select 聚类或微区中的子集，默认为NULL，即所有子区域
#'
#' @export
SelectZone <- function(ZoneFile, SampleName= "LN3", Collable="clusters", Select=NULL){
  
  library(dplyr)
  library(magrittr)
  library(rlang)
  
  filename <- file.path(ZoneFile, SampleName, "clusters_infor.csv")
  barcode_info <- read.table(filename, sep = ",", header = T,check.names = F,stringsAsFactors = F)
  barcode_info_sam <- dplyr::filter(barcode_info, sampleid==SampleName)
  
  # select
  Select <- Select %||% unique(barcode_info_sam[, Collable])
  Select %<>% as.character() %>% unique() %>% .[order(.)]
  if (!all(Select %in% barcode_info_sam[, Collable])) stop("selecting set must be all in clusters, please check it ~ ")
  result <- lapply(Select, function(x){paste0(SampleName,".", gsub("-[0-9]*","", barcode_info_sam[barcode_info_sam[, Collable]==x, "Barcode"]))})
  names(result) <- Select
  result
}

#' @export
auto_stemanalyst1 <- function(data,
                              circle_num = 5,
                              max_cluster = 50,
                              cluster = NULL,
                              filter_std = F,
                              min.std = 0,
                              standardise = T,
                              method = "cmeans", # "cmeans" "ufcl"
                              ...) {
  
  suppressMessages(library("lmbio"))
  suppressMessages(library("Mfuzz"))
  suppressMessages(library("dplyr"))
  set.seed(123)
  
  data1 <- data
  colnames(data1) <- paste0(colnames(data1), ".raw")
  colone_value <- data1[, 1]
  if (!standardise) {
    print("数据进行log2归一化处理")
    value_log <- function(x, na.rm = FALSE) (log2(x / colone_value))
    data2 <- dplyr::mutate_all(.tbl = data1, .funs = value_log)
    colnames(data2) <- paste0(colnames(data), ".Expression Change(log2(v(i)/v(0))")
    row.names(data2) <- rownames(data1)
  } else {
    print("数据进行Z-score标准化处理...")
    data2 <- t(scale(t(data1))) # 对row进行Z-score标准化
    colnames(data2) <- paste0(colnames(data), ".Expression Abundance(Z-score)")
  }
  
  print("热图总图保存")
  data3 <- as.data.frame(data2)
  names(data3) <- colnames(data)
  # meta::auto_heatmap(data3, 
  #                    mapname = "Summary-heatmap",
  #                    ...)
  
  data4 <- cbind(Feature = row.names(data2),
                 data2,
                 data1)
  
  savexlsx1(data = data4,
            filename = "分析原始数据.xlsx",
            sheet = "分析原始数据")
  
  # 需要matrix数据类型
  count_matrix <- data.matrix(data2)
  eset <- new("ExpressionSet", exprs = count_matrix)
  
  # 根据标准差去除样本间差异太小的基因
  if (!filter_std) {
    eset <- eset
  } else {
    eset <- Mfuzz::filter.std(eset, min.std = min.std)
  }
  #  评估出最佳的m值
  # m <- max(Mfuzz::mestimate(eset), log(ncol(count_matrix)**0.33, 5))
  m <- Mfuzz::mestimate(eset)
  
  # change m and c
  base::print(length(Mfuzz::mfuzz(eset, c = cluster, m = m, ...)$size))
  while(length(Mfuzz::mfuzz(eset, c = cluster, m = m, ...)$size)==0) {
    m <- m+0.2
  }
  while(length(Mfuzz::mfuzz(eset, c = cluster, m = m, ...)$size)!=cluster) {
    cluster <- length(Mfuzz::mfuzz(eset, c = cluster, m = m, ...)$size)
  }
  
  c <- cluster
  cl <- Mfuzz::mfuzz(eset, c = c, m = m, ...)
  
  # # # 聚类总图
  # n1 <- ceiling(sqrt(c + 1))
  # n2 <- ceiling((c + 1) / n1)
  # 
  # plotfile(mapname = "Cluster-summary",
  #          height = n1 * 4,
  #          width = n2 * 4,
  #          ...)
  # time.labels <- if (ncol(data)>10) rep(" ", ncol(data))else colnames(data)
  # Mfuzz::mfuzz.plot(eset, cl,
  #                   mfrow = c(n2, n1),
  #                   new.window = FALSE,
  #                   time.labels = time.labels)
  # 
  # mfuzzColorBar <- function (col, horizontal = FALSE, ...){
  #   require(marray) || stop("Library marray is required")
  #   if (missing(col)) {
  #     col <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700",
  #              "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00",
  #              "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00",
  #              "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40",
  #              "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7",
  #              "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF",
  #              "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF",
  #              "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF",
  #              "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF",
  #              "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF",
  #              "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078",
  #              "#FF0060", "#FF0048", "#FF0030", "#FF0018")
  #   }else if (length(col) > 1) {
  #   }else if (col == "fancy") {
  #     fancy.blue <- c(c(255:0), rep(0, length(c(255:0))),
  #                     rep(0, length(c(255:150))))
  #     fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
  #     fancy.red <- c(c(0:255), rep(255, length(c(255:0))),
  #                    c(255:150))
  #     col <- rgb(b = fancy.blue/255, g = fancy.green/255,
  #                r = fancy.red/255)
  #   }
  #   par(mar = c(5, 2, 4, 3) + 0.1)
  #   maColorBar(seq(0, 1, 0.01), col = col, horizontal = FALSE,
  #              k = 11, ...)
  # }
  # 
  # suppressWarnings(mfuzzColorBar(main = "Membership", cex.main = 1)) # 画lengend，默认dpi为300
  # plotsave()
  
  # # 保存每个cluster中的element个数
  cl.size <- data.frame("Cluster number" = paste0("Cluster ", seq(c)),
                        "Cluster amount" = cl$size, check.names = F)
  
  savexlsx1(data = cl.size,
            filename = "Cluster Statistics.xlsx",
            sheet = "Cluster Statistics")
  
  
  # 查看基因和cluster之间的membership
  membership <- cl$membership
  membership1 <- data.frame(membership)
  names(membership1) <- paste0("Cluster ", seq(c))
  
  membership1 <- cbind(Feature = row.names(membership1),
                       membership1)
  
  savexlsx1(data = membership1,
            filename = "Membership.xlsx",
            sheet = "Membership")
  
  # 提取cluster下的element,并保存成文件
  for (i in seq(c)) {
    cluster_i <- data.frame("Cluster number" = cl$cluster[cl$cluster == i],
                            check.names = F)
    cluster_i[, 1] <- paste0("Cluster ", cluster_i[, 1])
    print(paste0("保存Cluster ", i, "的热图"))
    data_i <- data3[row.names(cluster_i), , drop = F] # 筛选出该Cluster的metabolites
    # 画热图
    auto_heatmap(data = data_i,
                 mapname = paste0("Cluster ", i, "-heatmap"),
                 ...)
    
    cluster_i <- cbind(cluster_i, data4[row.names(cluster_i), ]) # 合并
    # 提取每个cluster的element的membership
    mem.ship <- data.frame("Membership" = membership[row.names(cluster_i), i])
    cluster_i <- cbind(cluster_i, mem.ship) # 合并
    # 保存分析结果
    print(paste0("保存Cluster ", i, "的结果"))
    savexlsx1(data = cluster_i,
              filename = paste0("Cluster ", i, ".xlsx"),
              sheet = paste0("Cluster ", i))
    
    # 查看属于cluster cores的基因list. cluster cores为membership > 0.7的metabolites
    acore.list <- Mfuzz::acore(eset, cl, min.acore = 0.7)
    acore.list[[i]] <- dplyr::rename(acore.list[[i]], Feature = NAME)
    acore.list[[i]] <- dplyr::rename(acore.list[[i]], Membership = MEM.SHIP)
    
    savexlsx1(data = acore.list[[i]],
              filename = paste0("Cluster ", i, ".xlsx"),
              sheet = paste0("Cluster ", i, " cores"))
    
    # 画图，线的颜色以各个代谢物在该cluster的membership来定
    choosedata <- data.frame(t(select(cluster_i, contains("Expression")))) # 抽提画图信息（时间轴，表达量）
    h <- length(rownames(choosedata))
    choosedata <- reshape2::melt(choosedata) # 将宽数据框变为长数据框， 便于画图
    choosedata["Samples"] <- colnames(data)
    choosedata["Membership"] <- rep(cluster_i$Membership, each = h, times = 1) # 添加membership，指示颜色
    
    cl_draw1(data = data,
             choosedata = choosedata,
             i = i,
             standardise = standardise,
             ...)
  }
}

# 画cluster的折线图
#' @export
cl_draw1 <- function(data,
                     choosedata,
                     i,
                     standardise,
                     width = 9,
                     height = 7,
                     xlab = "Spots",
                     ...) {
  suppressMessages(library("ggplot2"))
  theme_set(theme_classic())
  p <- ggplot(data = choosedata, aes(x = factor(choosedata$Samples, levels = colnames(data)),
                                     y = value)) +
    geom_line(aes(color = Membership, group = Membership)) +
    labs(title = paste0("Cluster ", i),
         subtitle = paste("contains ", length(choosedata$Membership) / length(colnames(data)), " Features" , sep = ""),
         x = xlab) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5), 
          axis.text.x = element_blank(), axis.ticks = element_blank(), 
          legend.text.align = 0.5)
  
  col <- c("#FFFF66", "#FFFF33", "#FFFF00", "#33FF66", "#33CC33", "#33CC00",
           "#0099FF", "#0066FF", "#6633CC", "#663399", "#663366", "#990099",
           "#FF0066", "#FF0033", "#CC0000")
  
  if (!standardise) {
    p <- p + scale_color_gradientn(colours = col, breaks = seq(0, 1, 0.1)) +
      scale_y_continuous(name = "Expression Change (log2(v(i)/v(0))")
  } else {
    p <- p + scale_color_gradientn(colours = col, breaks = seq(0, 1, 0.1)) +
      scale_y_continuous(name = "Expression Abundance (Z-score)")
  }
  
  # 保存图片
  print(paste("保存Cluster ", i, "的聚类图", sep = ""))
  ggplotsave(plot = p,
             mapname = paste0("Cluster ", i),
             width = width, height = height,
             ...)
}



