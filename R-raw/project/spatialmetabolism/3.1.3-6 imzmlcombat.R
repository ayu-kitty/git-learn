#!/opt/conda/bin/Rscript

#' 空代批次校正处理
#'
#' @param mass.range mz范围
#' @param resolution 分辨率
#' @param units 分辨率单位
#' @param imzmlpath 数据路径
#' @param savepath 保存路径
#' @param moderange 正负离子模式
#' @param batchtype 批次校正方法，harmony、anchors、combat
#' @param normalize 逻辑，是否标准化
#'
#' @export
Combatspace <- function(mass.range = NULL,
                        resolution = 5,
                        units = "ppm",
                        imzmlpath = "./sample/final/",
                        savepath = "./sample/adjustdata/",
                        moderange = c("neg","pos"),
                        batchtype = "anchors",
                        normalize = F,
                        sampletobatch = NULL,
                        batchlist = "slidename",
                        combatmod = NULL,
                        dealminus = "fill0",
                        qcname = "qc_data",
                        qcweight = 0.6,
                        groupweight = 0.6,
                        sampleweight = 0.9,
                        ...){
  suppressMessages(library("Seurat"))
  suppressMessages(library("Cardinal"))
  suppressMessages(library("harmony"))
  suppressMessages(library("sva"))
  
  print("---------批次校正中---------")
  
  for (mode in moderange) {
    
    filename <- list.files(path = imzmlpath,
                           pattern = paste0("-", mode, ".imzML$"),
                           full.names = F,
                           recursive = T)
    
    samplename <- gsub(pattern = paste0("-", mode, ".imzML"),
                       replacement = "",
                       x = filename)
    if(length(samplename) == 0){
      print(paste0(imzmlpath,"未找到",mode,"模式的数据"))
      next
    }
    
    if(any(grepl(pattern = paste0("^",qcname,".*"),x = samplename)) & is.null(combatmod)){
      combatmod <- data.frame("modlist" = "sample",
                              "samplename" = samplename,
                              stringsAsFactors = F)
      combatmod[grepl(pattern = paste0("^",qcname,".*"),x = samplename),"modlist"] <- "qc"
      print("~查看感兴趣分组")
      base::print(combatmod)
    }
    
    print("数据读取")
    i <- 1
    print(samplename[i])
    mse <- readMSIData(paste0(imzmlpath, samplename[i], "-", mode, ".imzML"),
                       attach.only = F,
                       mass.range = mass.range, resolution = resolution, units = units)
    spectradata <- as.data.frame(iData(mse, "intensity"))
    names(spectradata) <- paste(samplename[i],
                                mode, coord(mse)$x, coord(mse)$y,
                                sep = "-")
    row.names(spectradata) <- format(mz(mse), nsmall = 5, trim = T)
    
    sampleinfo <- readRDS(file = paste0(imzmlpath, samplename[i], "-", mode, ".rds"))
    
    if(!is.null(sampletobatch)){
      if(!all(colnames(sampletobatch) %in% c("samplename","batch"))){
        stop("sampletobatch参数的格式不对")
      }
      
      sampleinfo <- merge(sampleinfo,sampletobatch,by = "samplename",all.x = T,sort = F)
      if(any(is.na(sampleinfo$batch))){
        stop("sampletobatch参数对部分样本未提供批次组别")
      }
      batchlist <- "batch"
    }
    if(!is.null(combatmod)){
      if(!all(colnames(combatmod) %in% c("samplename","modlist"))){
        stop("combatmod参数的格式不对")
      }
      
      sampleinfo <- merge(sampleinfo,combatmod,by = "samplename",all.x = T,sort = F)
      if(any(is.na(sampleinfo$modlist))){
        stop("combatmod参数对部分样本未提供预期分组组别")
      }
    }
    
    row.names(sampleinfo) <- paste(sampleinfo$samplename,
                                   sampleinfo$mode,
                                   sampleinfo$x,
                                   sampleinfo$y,
                                   sep = "-")
    
    Seuratdata <- CreateSeuratObject(counts = spectradata,meta.data = sampleinfo,assay = "Meta")
    allSeuratdata <- Seuratdata
    
    if(length(samplename) > 1){
      for ( i in 2:length(samplename)) {
        print(samplename[i])
        mse <- readMSIData(paste0(imzmlpath, samplename[i], "-", mode, ".imzML"),
                           attach.only = F,
                           mass.range = mass.range, resolution = resolution, units = units)
        spectradata <- as.data.frame(iData(mse, "intensity"))
        names(spectradata) <- paste(samplename[i],
                                    mode, coord(mse)$x, coord(mse)$y,
                                    sep = "-")
        row.names(spectradata) <- format(mz(mse), nsmall = 5, trim = T)
        
        sampleinfo <- readRDS(file = paste0(imzmlpath, samplename[i], "-", mode, ".rds"))
        
        if(!is.null(sampletobatch)){
          sampleinfo <- merge(sampleinfo,sampletobatch,by = "samplename",all.x = T,sort = F)
          if(any(is.na(sampleinfo$batch))){
            stop("sampletobatch参数对部分样本未提供批次组别")
          }
        }
        if(!is.null(combatmod)){
          sampleinfo <- merge(sampleinfo,combatmod,by = "samplename",all.x = T,sort = F)
          if(any(is.na(sampleinfo$modlist))){
            stop("combatmod参数对部分样本未提供预期分组组别")
          }
        }
        
        row.names(sampleinfo) <- paste(sampleinfo$samplename,
                                       sampleinfo$mode,
                                       sampleinfo$x,
                                       sampleinfo$y,
                                       sep = "-")
        
        Seuratdata <- CreateSeuratObject(counts = spectradata,meta.data = sampleinfo,assay = "Meta")
        allSeuratdata <- merge(allSeuratdata,Seuratdata)
      }
    }
    
    print("特征提取")
    allSeuratdata <- FindVariableFeatures(allSeuratdata,verbose = F,nfeatures = 100000)
    
    if(normalize){
      print("归一化处理")
      allSeuratdata <- NormalizeData(allSeuratdata)
    }
    
    if(length(unique(allSeuratdata@meta.data[,batchlist]))>1){
      print("批次校正处理")
      if(batchtype == "harmony"){
        middledata <- GetAssayData(object = allSeuratdata)
        print("数据格式转化中")
        # alldata <- as.matrix(middledata)
        # alldata <- Matrix::as.matrix(middledata)
        alldata <- as_matrix(middledata)
        alldata[alldata == 0] <- 0.000001
        print("数据格式转化完成")
        meta_data <- allSeuratdata@meta.data
        # if("modlist" %in% colnames(meta_data)){
        #   dealdata <- HarmonyMatrix(data_mat = alldata,
        #                             meta_data =  meta_data,
        #                             vars_use = c(batchlist,"modlist"),
        #                             reference_values = "qc",
        #                             do_pca = FALSE)
        # }else{
          dealdata <- HarmonyMatrix(data_mat = alldata,
                                    meta_data =  meta_data,
                                    vars_use = batchlist,
                                    do_pca = FALSE)
        # }
      }else if(batchtype == "anchors"){
        meta_data <- allSeuratdata@meta.data
        allSeuratdatalist <- SplitObject(allSeuratdata,split.by = batchlist)
        
        features <- SelectIntegrationFeatures(object.list = allSeuratdatalist, nfeatures = 500)
        
        allSeuratdatalist <- lapply(X = allSeuratdatalist, FUN = function(x) {
          x <- ScaleData(x, features = features, verbose = FALSE)
          x <- RunPCA(x, features = features, verbose = FALSE)
        })
        
        # allSeuratdatalist2 <<- allSeuratdatalist
        allSeuratdataanchor <- FindIntegrationAnchors(allSeuratdatalist,
                                                      anchor.features = 100000,
                                                      k.anchor = 20,
                                                      reduction = "rpca")
        allSeuratdataanchor@anchor.features <- row.names(allSeuratdatalist[[1]])
        # allSeuratdataanchor2 <<- allSeuratdataanchor
        if("modlist" %in% colnames(meta_data)){
          allSeuratdataanchordata <- allSeuratdataanchor@anchors
          allSeuratdataanchordata <- allSeuratdataanchordata[allSeuratdataanchordata$score > min(c(qcweight,groupweight,sampleweight)),]
          i <- 1
          while (i <= dim(allSeuratdataanchordata)[1]) {
            if((allSeuratdataanchor@object.list[[allSeuratdataanchordata[i,"dataset1"]]]@meta.data[allSeuratdataanchordata[i,"cell1"],"modlist"] == "qc" &
                allSeuratdataanchor@object.list[[allSeuratdataanchordata[i,"dataset1"]]]@meta.data[allSeuratdataanchordata[i,"cell1"],"modlist"] == allSeuratdataanchor@object.list[[allSeuratdataanchordata[i,"dataset2"]]]@meta.data[allSeuratdataanchordata[i,"cell2"],"modlist"] &
                allSeuratdataanchordata[i,"score"] >= qcweight) |
               (allSeuratdataanchor@object.list[[allSeuratdataanchordata[i,"dataset1"]]]@meta.data[allSeuratdataanchordata[i,"cell1"],"modlist"] != "qc" &
                allSeuratdataanchor@object.list[[allSeuratdataanchordata[i,"dataset2"]]]@meta.data[allSeuratdataanchordata[i,"cell2"],"modlist"] != "qc" &
                allSeuratdataanchordata[i,"score"] >= sampleweight)|
               (allSeuratdataanchor@object.list[[allSeuratdataanchordata[i,"dataset1"]]]@meta.data[allSeuratdataanchordata[i,"cell1"],"modlist"] != "qc" &
                allSeuratdataanchor@object.list[[allSeuratdataanchordata[i,"dataset1"]]]@meta.data[allSeuratdataanchordata[i,"cell1"],"modlist"] == allSeuratdataanchor@object.list[[allSeuratdataanchordata[i,"dataset2"]]]@meta.data[allSeuratdataanchordata[i,"cell2"],"modlist"] &
                allSeuratdataanchordata[i,"score"] >= groupweight)){
              
              if(allSeuratdataanchor@object.list[[allSeuratdataanchordata[i,"dataset1"]]]@meta.data[allSeuratdataanchordata[i,"cell1"],"modlist"] == "qc"){
                allSeuratdataanchordata[i,"type"] <- "qc"
              }else{
                allSeuratdataanchordata[i,"type"] <- "sample"
              }
              
            }else{
              allSeuratdataanchordata[i,"type"] <- "del"
            }
            i <- i +1
          }
          
          allSeuratdataanchordata <- allSeuratdataanchordata[allSeuratdataanchordata$type!="del",]
          
          allSeuratdataanchor@anchors <- allSeuratdataanchordata[,-dim(allSeuratdataanchordata)[2]]
        }else{
          allSeuratdataanchordata <- allSeuratdataanchor@anchors
          allSeuratdataanchordata <- allSeuratdataanchordata[allSeuratdataanchordata$score > min(c(qcweight,groupweight,sampleweight)),]
          allSeuratdataanchor@anchors <- allSeuratdataanchordata
        }
        
        allSeuratdatanew <- IntegrateData(allSeuratdataanchor,normalization.method = "LogNormalize")
        
        middledata <- GetAssayData(object = allSeuratdatanew,assay = "integrated")
        # middledata2 <<- middledata
        dealdata <- as_matrix(middledata)
        dealdata <- dealdata[order(as.numeric(row.names(dealdata))),]
        
      }else if(batchtype == "combat"){
        middledata <- GetAssayData(object = allSeuratdata)
        print("数据格式转化中")
        # alldata <- as.matrix(middledata)
        # alldata <- Matrix::as.matrix(middledata)
        alldata <- as_matrix(middledata)
        print("数据格式转化完成")
        meta_data <- allSeuratdata@meta.data
        
        if("modlist" %in% colnames(meta_data)){
          mod <- model.matrix(~as.factor(modlist), data=meta_data)
        }else{
          mod <- NULL
        }
        
        dealdata <- ComBat(dat = alldata,
                           batch = meta_data[,batchlist],
                           mod = mod,
                           par.prior = TRUE,
                           prior.plots = FALSE)
        
      }else if(batchtype == "none"){
        print("~不进行批次矫正")
        
        dealdata <- as.matrix(GetAssayData(object = allSeuratdata))
        meta_data <- allSeuratdata@meta.data
      }else{
        print("请输入正确的批次矫正方法")
      }
    }else{
      print("~上机玻片数为1，不进行批次矫正")
      
      dealdata <- as.matrix(GetAssayData(object = allSeuratdata))
      meta_data <- allSeuratdata@meta.data
    }
    
    # [x] 将批次校正中小于0的数据赋值
    if(!normalize){
      if(min(dealdata) < 0){
        if(dealminus == "fill0"){
          dealdata[dealdata < 0] <- 0
        }else if(dealminus == "addmin"){
          dealdata <- dealdata-min(dealdata)
        }else if(dealminus == "fillraw"){
          middledata <- GetAssayData(object = allSeuratdata)
          alldata <- as_matrix(middledata)
          dealdata[dealdata < 0] <- alldata[dealdata < 0]
          dealdata[alldata == 0] <- 0
        }
      }
    }
    
    print("数据保存")
    for ( i in 1:length(samplename)){
      print(samplename[i])
      sample_meta_data <- meta_data[meta_data$samplename == samplename[i],]
      sampledata <- dealdata[,rownames(sample_meta_data)]
      
      mse <- readMSIData(paste0(imzmlpath, samplename[i], "-", mode, ".imzML"),
                         attach.only = F,
                         mass.range = mass.range, resolution = resolution, units = units)
      spectra(mse) <- sampledata
      
      # 数据保存
      filename <- paste0(savepath,"/",samplename[i],"-",mode,".imzML")
      filename_info <- paste0(savepath, samplename[i], "-", mode, ".rds")
      
      saveimzml(data = mse,filename = filename)
      saverds(data = sample_meta_data,filename = filename_info)
      
    }
  }
}

#' @export
as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

#' @export
splitdata <- function(data) {
  id <- strsplit(as.character(data[2]), ",")[[1]]
  return(as.data.frame(cbind(data[1], id)))
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-i","--imzmlpath",default = "./sample/final/", help = "imzml原始文件路径,默认./sample/final/")
  parser$add_argument("-m","--massrange",default = NULL,type= "double",help = "mz范围,默认NULL",
                      nargs = 2,dest = "mass.range")
  parser$add_argument("-re","--resolution",default = 5, type= "double",help = "分辨率,默认5")
  parser$add_argument("-u","--units",default = "ppm",help = "分辨率单位,默认ppm")
  parser$add_argument("-mr","--moderange",default = c("neg","pos"),help = "正负离子模式",nargs = "+",
                      choices = c("neg","pos"))
  
  parser$add_argument("-s","--savepath",default = "./sample/adjustdata/", help = "imzml数据保存路径,默认./sample/adjustdata/")
  
  parser$add_argument("-nl","--normalize",default = F, help = "是否标准化",action='store_true')
  parser$add_argument("-bt","--batchtype",default = "anchors", help = "批次校正方式,包含harmony、anchors、combat、none")
  parser$add_argument("-bl","--batchlist",default = "slidename",help = "批次来源,包含slidename、samplename")
  parser$add_argument("-sb","--sampletobatch",default = NULL,help = "根据样本提供批次,如1:A1,B1,C1 2:A2,B2,C2",
                      nargs = "+")
  parser$add_argument("-cm","--combatmod",default = NULL,help = "根据样本提供感兴趣分组,如A:A1,A2 B:B1,B2 C:C1,C2",
                      nargs = "+")
  parser$add_argument("-dm","--dealminus",default = "fillraw",help = "填充负数的方式,fill0(填充0)、addmin(加上最小值的绝对值)、none(不填充),fillraw(填充原始值)")
  parser$add_argument("-qw","--qcweight",default = 0.6, type= "double",help = "qc权重,越小qc权重越大,默认0.4")
  parser$add_argument("-gw","--groupweight",default = 0.6, type= "double",help = "分组权重,越小分组权重越大,默认0.4")
  parser$add_argument("-sw","--sampleweight",default = 0.9, type= "double",help = "样本权重,越小样本权重越大,默认0.8")
  
  args <- parser$parse_args()
  
  if(!is.null(args$sampletobatch)){
    txt <- args$sampletobatch
    txt <- strsplit(txt,split = ":")
    txt <- Reduce(rbind,txt)
    txt <- as.data.frame(txt)
    txt <- apply(txt, 1, splitdata)
    txt <- Reduce(rbind, txt)
    colnames(txt) <- c("batch","samplename")
    args$sampletobatch <- txt
    print("~查看批次组别")
    print(txt)
    if(any(duplicated(txt$samplename))){
      stop("样本有重复分组")
    }
  }
  
  if(!is.null(args$combatmod)){
    txt <- args$combatmod
    txt <- strsplit(txt,split = ":")
    txt <- Reduce(rbind,txt)
    txt <- as.data.frame(txt)
    txt <- apply(txt, 1, splitdata)
    txt <- Reduce(rbind, txt)
    colnames(txt) <- c("modlist","samplename")
    args$combatmod <- txt
    print("~查看感兴趣分组")
    print(txt)
    if(any(duplicated(txt$samplename))){
      stop("样本有重复分组")
    }
  }
  
  writeinfo()
  createdir(filename = args$savepath,linkdir = T)
  
  result <- do.call(what = Combatspace,args = args)
  
  writeinfo(endtime = T)
  
}
