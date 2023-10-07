#!/opt/conda/bin/Rscript

#' @export
progmt <- function(inputpath = "../background/",inputfile = "F"){
  if("data.frame" %in% class(inputfile)){
    omty <- TRUE
  }else if("character" %in% class(inputfile)){
    if(inputfile == "F"){
      omty <- FALSE
    }else{omty <- TRUE}
  }
  pacman::p_load(dplyr,openxlsx,stringr)
  if(!file.exists(paste0(inputpath,"/go.gmt"))){
    if(file.exists(paste0(inputpath,"/gene_go.backgroud.xls"))){
      go<-read.delim(paste0(inputpath,"/gene_go.backgroud.xls"), sep="\t", header=F, quote="")
      go[,2]<-as.character(go[,2])
      go[,3]<-as.character(go[,3])
      cbind(unlist(strsplit(go[,3],"[|]")),unlist(strsplit(go[,2],",")))->pat
      sapply(go[,2],function(x){length(str_extract_all(x,",")[[1]])+1})%>%as.numeric()->pn
      cbind(pat,rep(go[,1],pn))->patp
      cbind(patp,apply(patp,1,function(x){(paste0(x[1],"(",x[2],")"))}))->patp
      unique(patp[,4])->pathway
      
      sink(paste0(inputpath,"/go.gmt"))
      for(i in 1:length(pathway)){
        cat(cat(pathway[i],"",patp[patp[,4]==pathway[i],3],sep="\t"),sep="\n")
      }
      sink()
    }
  }
  if(omty){
    if(file.exists(paste0(inputpath,"/background.xlsx"))){
      back <- readdata(paste0(inputpath,"/background.xlsx"),sheet=2)
      meta <- readdata(inputfile)
      
      if("Metabolites" %in% colnames(meta) & "KEGG" %in% colnames(meta)){
        mk <- select(meta,c("Metabolites","KEGG"))
        if(any(grepl(pattern = ";",mk$Metabolites))){
          meta_2 <- mk[0,]
          meta_1 <- mk
          for ( j in 1:dim(meta_1)[1]) {
            meta_1_1 <- unlist(strsplit(split = ";\n",x = meta_1[j,1]))
            meta_1_1 <- unlist(strsplit(split = "; ",x = meta_1_1))
            
            kegg_1_1 <- unlist(strsplit(split = ";\n",x = meta_1[j,2]))
            kegg_1_1 <- unlist(strsplit(split = "; ",x = kegg_1_1))
            
            if(length(meta_1_1) == 0){
              next
            }
            for ( k in 1:length(meta_1_1)) {
              if(meta_1_1[k] == ""){
                next
              }
              meta_1_1_1 <- meta_1[j,,drop = F]
              meta_1_1_1[1,1] <- meta_1_1[k]
              meta_1_1_1[1,2] <- kegg_1_1[k]
              meta_2 <- rbind(meta_2,meta_1_1_1)
              # break
            }
          }
          mk <- meta_2
        }
      }else if("Metabolites" %in% colnames(meta)){
        mk <- select(meta,c("Metabolites"))
        if(any(grepl(pattern = ";",mk$Metabolites))){
          mk <- data.frame("Metabolites"=unlist(strsplit(x = mk$Metabolites,split = ";\n")))
          mk <- data.frame("Metabolites"=unlist(strsplit(x = mk$Metabolites,split = "; ")))
        }
        mk <- getmetainfo(data = mk,idlist = "Metabolites",needlist = "KEGG")
      }else{
        stop("无Metabolites列")
      }
      
      mk <- na.omit(mk)
      mk <- mk[mk[,2] != "",]
      mk <- mk[!duplicated(mk[,1]),]
      # kegg <- left_join(mk,back,by="KEGG")
      # kegg <- kegg[,-2]
      kegg <- merge(x = mk,y = back,by="KEGG",all.y = T)
      kegg[is.na(kegg[,2]),"Metabolites"] <- kegg[is.na(kegg[,2]),"KEGG"]
      kegg <- kegg[,-1]
      kegg <- na.omit(kegg)
      kegg[,2] <- as.character(kegg[,2])
      kegg[,3] <- as.character(kegg[,3])
      cbind(unlist(strsplit(kegg[,3],"[|]")),unlist(strsplit(kegg[,2],","))) -> pat
      sapply(kegg[,2],function(x){length(str_extract_all(x,",")[[1]])+1})%>%as.numeric() -> pn
      cbind(pat,rep(kegg[,1],pn)) -> patp
      cbind(patp,apply(patp,1,function(x){(paste0(x[1],"(",x[2],")"))})) -> patp
      unique(patp[,4]) -> pathway
      
      sink(paste0(inputpath,"/kegg.gmt"))
      for(i in 1:length(pathway)){
        cat(cat(pathway[i],"",patp[patp[,4]==pathway[i],3],sep="\t"),sep="\n")
      }
      sink()
    }
  }else if(!file.exists(paste0(inputpath,"/kegg.gmt"))){
    if(file.exists(paste0(inputpath,"/gene_kegg.backgroud.xls"))){
      kegg<-read.delim(paste0(inputpath,"/gene_kegg.backgroud.xls"), sep="\t", header=F, quote="")
      kegg[,2]<-as.character(kegg[,2])
      kegg[,3]<-as.character(kegg[,3])
      cbind(unlist(strsplit(kegg[,3],"[|]")),unlist(strsplit(kegg[,2],",")))->pat
      sapply(kegg[,2],function(x){length(str_extract_all(x,",")[[1]])+1})%>%as.numeric()->pn
      cbind(pat,rep(kegg[,1],pn))->patp
      cbind(patp,apply(patp,1,function(x){(paste0(x[1],"(",x[2],")"))}))->patp
      unique(patp[,4])->pathway
      
      sink(paste0(inputpath,"/kegg.gmt"))
      for(i in 1:length(pathway)){
        cat(cat(pathway[i],"",patp[patp[,4]==pathway[i],3],sep="\t"),sep="\n")
      }
      sink()
    }
  }
  
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-ip","--inputpath",type="character", default="../background/", help="背景文件存放路径，默认../background/", metavar="character")
  parser$add_argument("-if","--inputfile",type="character", default="F", help="代谢使用，数据矩阵文件.xlsx", metavar="character")
  
  args <- parser$parse_args()
  progmt <- do.call(progmt,args = args)
}