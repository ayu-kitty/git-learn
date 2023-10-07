#!/opt/conda/bin/Rscript

#' 空转barcode替换像素点并整理分组信息
#' 
#' @param workpath 文件读取路径
#' @param mode 需要处理的正负离子模式
#' @param savepath 文件保存路径
#' @param savename 文件保存名字
#' @param list 需要处理的比较组
#' @param ionmerge 正负离子是否需要合并
#' @param barcodelist 转录提供的亚型barcode列表
#' @param barcodefile 用于获取实际barcode对照的file,neg,pos,all均可,推荐有all用all
#'
#' @export
barcodetosp_dataextraction <- function(workpath = "/data/nas/176/空代-空代实验组-jpp/项目数据/2022/DLM20227913/空转空代联合-0804/售后1/2.uniondata",
                                       mode = c("neg","pos"),
                                       savepath = "./",
                                       savename = "barcodemeta",
                                       list = c(),
                                       ionmerge = T,
                                       barcodelist = "./6个样本4种细胞亚型Kiss1阳性barcode信息.xlsx",
                                       barcodefile = "barcodemeta-all.txt"){
  
  
  barcodetosp_merge(workpath = workpath,
                    mode = mode,
                    savepath = savepath,
                    savename = savename,
                    list = list,
                    ionmerge = ionmerge)
  
  barcodetosp_samplegroup(workpath = workpath,
                          savepath = savepath,
                          list = list,
                          barcodelist = barcodelist,
                          barcodefile = barcodefile)
}


#' 空转barcode替换像素点
#'
#' @param workpath 文件读取路径
#' @param mode 需要处理的正负离子模式
#' @param savepath 文件保存路径
#' @param savename 文件保存名字
#' @param list 需要处理的比较组
#' @param ionmerge 正负离子是否需要合并
#'
#' @export
barcodetosp_merge <- function(workpath = "/data/nas/176/空代-空代实验组-jpp/项目数据/2022/DLM20227913/空转空代联合-0804/售后1/2.uniondata",
                              mode = c("neg","pos"),
                              savepath = "./",
                              savename = "barcodemeta",
                              list = c(),
                              ionmerge = T){
  # 进入工作路径
  setwd(workpath)
  # 获取所有需要处理的比较组
  if (is.null(list)){
    list <- list.files("./",full.names = T,pattern = "_")
  }else{list = list}
  
  
  # 比较组循环
  for (k in 1:length(list)){
    setwd(list[k])
    # 正负离子循环
    for (i in mode){
      # 读取每个比较组每个正负模式下的samplename和metadata文件
      samplelist <- lmbio::readdata(paste0("./",i,"/samplename.txt"))
      metaata <- lmbio::readdata(paste0("./",i,"/metadata.txt"))
      metadata2 <- metaata
      for (j in 1:length(samplelist$metaname)){
        col <- colnames(metadata2)
        col <- gsub(samplelist$metaname[j],samplelist$transname[j],col)
        colnames(metadata2) <- col
      }
      # 保存barcode对应替换后的矩阵
      savetxt(metadata2,filename = paste0(savepath,savename,"-",i,".txt"))
    }
    
    # 是否需要正负离子合并
    if (ionmerge){
      # 正负合并
      neg <- lmbio::readdata(filename = "barcodemeta-neg.txt")
      neg$ionmode <- "neg"
      pos <- lmbio::readdata(filename = "barcodemeta-pos.txt")
      pos$ionmode <- "pos"
      neg <- t(neg)
      pos <- t(pos)
      
      all <- merge(neg,pos,by="row.names")
      
      all <- t(all)
      
      all <- data.frame(all)
      colnames(all) <- all[1,]
      all <- all[-1,]
      
      all <- all[,c(grep("Adduct",colnames(all)):length(all),grep("_",colnames(all)))]
      
      savetxt(data = all,filename = paste0(savepath,savename,"-all.txt"))
      
    }
    setwd(workpath)
  }
}




#' 根据亚型获取分组情况
#'
#' @param workpath 文件读取路径
#' @param savepath 文件保存路径
#' @param list 需要处理的比较组
#' @param barcodelist 转录提供的亚型barcode列表
#' @param barcodefile 用于获取实际barcode对照的file,neg,pos,all均可,推荐有all用all
#'
#' @export
barcodetosp_samplegroup <- function(workpath = "/data/nas/176/空代-空代实验组-jpp/项目数据/2022/DLM20227913/空转空代联合-0804/售后1/2.uniondata",
                                    savepath = "./",
                                    list = c(),
                                    barcodelist = "./6个样本4种细胞亚型Kiss1阳性barcode信息.xlsx",
                                    barcodefile = "barcodemeta-all.txt"){
  
  setwd(workpath)
  if (is.null(list)){
    list <- list.files("./",pattern = "_")
  }else{list = list}
  
  
  
  for (i in 1:length(list)){
    #barcode处理,获得所有有信息列表中的barcode
    group1 <- lmbio::readxlsx(filename = barcodelist,sheet = list[i])
    setwd(list[i])
    #group2 <- c(group1$Astrocytes,group1$Microglia,group1$Neuron,group1$`Oligo+OPC`)
    group3 <- group1[!is.na(group1) & group1 != ""]
    
    #给barcode去重
    group4 <- unique(group3)
    
    #获取所有的barcode,每个all中的allname1
    allname <- colnames(lmbio::readdata(paste0("./",barcodefile)))
    allname1 <- allname[grep("_",allname)]
    allname1 <- gsub("\\.","\\-",allname1)
    
    #获取other组difgroup2
    difgroup <- setdiff(allname1,group4)
    difgroup2 <- data.frame(barcode = gsub("\\.","\\-",difgroup),group = "Other")
    
    for (j  in 1:length(colnames(group1))){
      areaname <- colnames(group1)[j]
      
      print(paste0("~分析组:",list[i],"-第",j,"个-区域名:",areaname))
      
      # 去目标亚型和数据集的交集
      area <- group1[,areaname][!is.na(group1[,areaname]) & group1[,areaname] != ""]
      area2 <- intersect(area,allname1)
      
      if (length(area2)>9){
        areagroup <- data.frame(barcode = area2,group = areaname)
      }else{
        areagroup <- NULL
      }
      
      # 非此次比较选区
      more1 <- union(area2,difgroup)
      more2 <- setdiff(allname1,more1)
      
      
      if (length(more2)>0){
        moregroup <- data.frame(barcode = more2,group = "more")
      }else{
        moregroup <- NULL
      }
      
      
      endgroup <- rbind(areagroup,difgroup2,moregroup)
      savexlsx(endgroup,filename = paste0(savepath,"分组情况.xlsx"),sheet = areaname)
      
    }
    
    setwd("../")
    
  }
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-wp","--workpath",default = "./2.uniondata", help = "文件读取路径")
  parser$add_argument("-m","--mode",default = c("neg","pos"), help = "需要处理的正负离子模式")
  parser$add_argument("-sp","--savepath",default = "./", help = "文件保存路径")
  parser$add_argument("-sn","--savename",default = "barcodemeta", help = "文件保存名字")
  parser$add_argument("-l","--list",default = c(), help = "需要处理的比较组")
  parser$add_argument("-i","--ionmerge",default = T, help = "正负离子是否需要合并",action ='store_false')
  parser$add_argument("-bl","--barcodelist",default = "转录提供的亚型barcode列表")
  parser$add_argument("-bf","--barcodefile",default = "barcodemeta-all.txt", help = "用于获取实际barcode对照的file,neg,pos,all均可,推荐有all用all")
  
  args <- parser$parse_args()
  
  writeinfo()
  
  mulargs <- do.call(what = barcodetosp_dataextraction,args = args)
  
  writeinfo(endtime = T)
  
}

