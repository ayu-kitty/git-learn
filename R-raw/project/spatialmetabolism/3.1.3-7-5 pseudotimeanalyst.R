#!/opt/conda/bin/Rscript

#' pseudotimeanalyst
#'
#' 拟时序分析
#'
#' @param flitermethod 挑选mz方法，包含monocle,Seurat,filtermz,none
#' @param filtermz 挑选mz方法为filtermz时，提供需要的mz
#' @param reduction_method 降维方法
#' @param samplename 样本名称
#' @param mode 离子模式
#' @param rdspath sscc中间rds文件路径
#' @param savepath 保存路径
#' @param ... 见[monocle::orderCells()]
#'
#' @export
pseudotimeanalyst <- function(flitermethod = "monocle", # monocle,Seurat,filtermz,none
                              filtermz = NULL,
                              reduction_method = "DDRTree",
                              samplename = "1462",
                              mode = "neg",
                              rdspath = "./sample/cluster/sscc/",
                              savepath = "./sample/cluster/pseudotime/",
                              norm_method = "none",
                              selectcluster = NULL,
                              maxpoint = 20000,
                              filename = paste0(rdspath,samplename,"-",mode,"-data.rds"),
                              savedata = T,
                              ...){

  # 加载monocle
  suppressMessages(library("monocle"))

  # 读取sample/cluster/sscc目录下rds文件的数据
  # filename = paste0(rdspath,samplename,"-",mode,"-data.rds")
  clusterdata <- readdata(filename = filename)

  # 获取表达矩阵，并将表达矩阵转化为matrix格式
  cellData <- as.matrix(clusterdata$spectradata)
  # 像素点对应信息
  pdata <- clusterdata$plotdata
  rownames(pdata) <- pdata$Name

  if(dim(cellData)[2] > maxpoint){
    print(paste0("由于像素点数大于",maxpoint,",减少像素点至",maxpoint))
    selectname <- sample(colnames(cellData),maxpoint)
    cellData <- cellData[,colnames(cellData) %in% selectname]
    pdata <- pdata[pdata$Name %in% selectname,]
  }
  newcellData <- cellData
  rownames(cellData) <- paste0("mz:",rownames(cellData))

  # mz对应信息
  fdata <- data.frame(gene_short_name = rownames(cellData),row.names = rownames(cellData))

  if(is.null(selectcluster)){
  }else{
    pdata <- pdata[pdata$Cluster %in% selectcluster,]
    cellData <- cellData[,pdata$Name]
  }

  print(paste0("像素点数量:",dim(cellData)[2]))

  # 将pdata和fdata格式转化为AnnotatedDataFrame
  newpdata <- new("AnnotatedDataFrame",data = pdata)
  newfdata <- new("AnnotatedDataFrame",data = fdata)
  # 创建CellDataSet格式数据
  print("创建CellDataSet格式数据")
  cds <- newCellDataSet(cellData = as(cellData,"sparseMatrix"),
                        phenoData = newpdata,
                        featureData = newfdata,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())

  # 统计
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- detectGenes(cds,min_expr = 0.1)

  # monocle高变基因筛选
  if(flitermethod == "monocle"){
    disp_table <- dispersionTable(cds)
    mzrange <- subset(disp_table,mean_expression >= 0.1 & dispersion_empirical >= 1*dispersion_fit)$gene_id
  }else if(flitermethod == "Seurat"){
    suppressMessages(library("Seurat"))
    Seuratdata <- CreateSeuratObject(counts = cellData,meta.data = pdata,assay = "Meta")
    Seuratdata <- FindVariableFeatures(Seuratdata,verbose = F)
    mzrange <-VariableFeatures(Seuratdata)
  }else if(flitermethod == "filtermz"){
    if(is.null(filtermz)){
      mzrange <- rownames(cellData)
    }else{
      mzrange <- paste0("mz:",format(filtermz,nsmall = 5,trim = T))
    }
  }else if(flitermethod == "none"){
    mzrange <- rownames(cellData)
  }else{
    print("mz筛选方法选择不对，默认使用所有mz")
    mzrange <- rownames(cellData)
  }
  print("高变mz筛选")
  print(paste0("高变mz数量:",length(mzrange)))
  cds <- setOrderingFilter(cds,ordering_genes = mzrange)
  print("降维分析")
  cds <- reduceDimension(cds,
                         max_components = 2,
                         reduction_method = reduction_method,
                         norm_method = norm_method)
  print("拟时序分析")
  cds <- orderCells(cds,...)

  alldata <- list(spectradata = newcellData,
                  cds = cds,
                  pdata = pData(cds),
                  fdata = fData(cds),
                  plotdata = pData(cds))
  
  if(savedata){
    cdsfilename <- paste0(savepath,samplename,"-",mode,"-data.rds")
    saverds(data = alldata,filename = cdsfilename)
  }else{
    return(alldata)
  }
}

#' mulpseudotimeanalyst
#'
#' 批量拟时序分析
#'
#' @param rdspath 数据路径
#' @param moderange 正负离子模式
#' @param anamoudle 分析模块
#' @param wantsample 想处理的样本名
#' @param delsample 不想处理的样本名
#' @param ... 见anamoudle函数
#'
#' @export
mulpseudotimeanalyst <- function(rdspath = "./sample/cluster/sscc/",
                                 moderange = c("neg", "pos"),
                                 anamoudle = pseudotimeanalyst,
                                 wantsample = NULL,
                                 delsample = c("bg_data","qc_data"),
                                 ...){
  for (mode in moderange) {

    # 获取.imzML文件绝对路径
    filename <- list.files(path = rdspath,
                           pattern = paste0("-", mode, "-data.rds$"),
                           full.names = F,
                           recursive = T)
    if(length(filename) == 0){
      print(paste0("在",rdspath,"目录下未找到", mode, "模式的rds文件"))
      next
    }

    # 获取样品名
    samplename <- gsub(pattern = paste0("-", mode, "-data.rds"),
                       replacement = "",
                       x = filename)

    if(!is.null(wantsample)){
      samplename <- samplename[samplename %in% wantsample]
    }

    if(!is.null(delsample)){
      for(delsample2 in delsample) {
        samplename <- samplename[!grepl(pattern = paste0("^",delsample2,".*"),x = samplename)]
      }
    }

    if (length(samplename) > 0) {

      for ( i in 1:length(samplename)) {

        anamoudle(samplename = samplename[i],
                  mode = mode,
                  rdspath = rdspath,
                  ...)

        gc(verbose = T)
      }

    } else {
      print(paste0("在",rdspath,"目录下未找到", mode, "模式的rds文件"))
    }
  }
}

