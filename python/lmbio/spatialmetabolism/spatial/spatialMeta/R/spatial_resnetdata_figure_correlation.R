spatial_resnetdata_figcor <- function(imzmlfile,
                                      resnetdata,
                                      mode = ifelse(grepl(pattern = "-neg\\.",basename(imzmlfile)),"neg","pos"),
                                      samplename = gsub(pattern = paste0("-",mode,"\\..*$"),replacement = "",x = basename(imzmlfile)),
                                      selectmz = NULL,
                                      savepath = "./sample/fig/result/resnet/",
                                      corfilter = 0.7,
                                      corfiltertype = "+",
                                      pfilter = 0.05,
                                      adjust = "BH",
                                      cormethod = "pearson",
                                      use = "pairwise",
                                      xname = "mz-x",
                                      yname = "mz-y",
                                      ...){
  library(lmbio)
  resnet  <-  resnetdata

  if(is.null(selectmz)){
    x <- resnet
    y <- NULL
  }else{
    # select mz
    x <- resnet
    x <- x[selectmz,,drop=F]
    # rest mz
    y <- resnet
    y <- resnet[row.names(y) != selectmz,,drop=F]
  }
  # calc corrdata
  corrdata <- corrcal(x = x,
                      y = y,
                      corfilter = corfilter,
                      corfiltertype = corfiltertype,
                      pfilter = pfilter,
                      adjust = adjust,
                      method = cormethod,
                      use = use,
                      xname = xname,
                      yname = yname)

  # filter by p value & r2
  linkdata <- corrdata$filterdata$linkdata
  linkdata2 <- linkdata
  # global list
  linkdata2 <<- linkdata2

  if(dim(linkdata)[1] == 0){
    return(corrdata)
  }

  # 绘制总图
  if(is.null(selectmz)){
    # x,y
    fitmz <- c(linkdata[,1],linkdata[,2])
    fitmz <- unique(fitmz)
    x <- resnet
    x <- x[fitmz,,drop=F]
    y <- NULL
    write.table(x = x,
                file = paste0(savepath,"/","rawdata-",samplename,"-",mode,".xls"),
                sep = "\t",
                row.names = T)
  } else {
    x <- resnet
    x <- x[selectmz,,drop=F]

    fitmz <- linkdata[,2]
    fitmz <- unique(fitmz)

    y <- resnet
    y <- y[fitmz,,drop=F]
    # select mz
    write.table(x = x,
                file = paste0(savepath,"/","rawdata-",samplename,"-",mode,"-select.xls"),
                sep = "\t",
                row.names = T)
    # high corr y
    write.table(x = y,
                file = paste0(savepath,"/","rawdata-",samplename,"-",mode,"-y.xls"),
                sep = "\t",
                row.names = T)

  }
  # 计算相关性，并绘图 network + correlation
  # corrdatamz <- datatocorr(x = x,
  #                          y = y,
  #                          networkmapname = paste0("corrnetwork-",samplename,"-",mode),
  #                          plotmapname = paste0("corrplot-",samplename,"-",mode),
  #                          savepath = savepath,
  #                          corfilter = corfilter,
  #                          corfiltertype = corfiltertype,
  #                          pfilter = pfilter,
  #                          adjust = adjust,
  #                          cormethod = cormethod,
  #                          use = use,
  #                          xname = xname,
  #                          yname = yname,
  #                          ...)
  #
  # write.table(x = corrdatamz$framedata$linkdata,
  #             file = paste0(savepath,"/","corrplot-",samplename,"-",mode,".xls"),
  #             sep = "\t",
  #             row.names = F)
  #
  # write.table(x = corrdatamz$filterdata$linkdata,
  #             file = paste0(savepath,"/","corrnetwork-",samplename,"-",mode,".xls"),
  #             sep = "\t",
  #             row.names = F)

  # HE img
  fig <- row.names(resnet)[length(row.names(resnet))]

  linkdata3 <- linkdata2[linkdata2[,1] == fig | linkdata2[,2] == fig,]
  fitmz <- c(linkdata3[,1],linkdata3[,2])
  fitmz <- unique(fitmz)

  x <- resnet[fitmz,,drop=F]
  y <- NULL
  createdir(paste0(savepath,"/",samplename,"-",mode,"/",fig))

  write.table(x = x,
              file = paste0(savepath,"/",samplename,"-",mode,"/",fig,"/","rawdata.xls"),
              sep = "\t",
              row.names = T)

  corrdatamz <- datatocorr(x = x,
                           y = y,
                           networkmapname = "corrnetwork",
                           plotmapname = "corrplot",
                           savepath = paste0(savepath,"/",samplename,"-",mode,"/",fig,"/"),
                           corfilter = corfilter,
                           corfiltertype = corfiltertype,
                           pfilter = pfilter,
                           adjust = adjust,
                           cormethod = cormethod,
                           use = use,
                           xname = xname,
                           yname = yname,
                           ...)

  write.table(x = corrdatamz$framedata$linkdata,
              file = paste0(savepath,"/",samplename,"-",mode,"/",fig,"/","corrplot.xls"),
              sep = "\t",
              row.names = F)

  write.table(x = corrdatamz$filterdata$linkdata,
              file = paste0(savepath,"/",samplename,"-",mode,"/",fig,"/","corrnetwork.xls"),
              sep = "\t",
              row.names = F)

  filtermz <- fitmz[-which(fitmz == fig)]

  imzmlimage(filename = imzmlfile,
             mapname = paste0("Intensity-",samplename,"-",mode),
             savepath = paste0(savepath,"/",samplename,"-",mode,"/",fig,"/"),
             mapmz = T,
             filtermz = filtermz)

  if (is.null(selectmz)){
    # do nothing
  }
  else {
    linkdata3 <- linkdata2[linkdata2[,1] == fig | linkdata2[,2] == fig,]
    # selectMZ
    linkdata4 <- linkdata3[(linkdata3[,1] %in% selectmz | linkdata3[,2] %in% selectmz),]
    fitmz <- c(linkdata4[,1],linkdata4[,2])
    fitmz <- unique(fitmz)
    x <- resnet[fitmz,,drop=F]
    y <- NULL

    corrdatamz <- datatocorr(x = x,
                         y = y,
                         networkmapname = "corrnetwork",
                         plotmapname = "corrplot",
                         savepath = paste0(savepath,"/",samplename,"-",mode,"-select/"),
                         corfilter = corfilter,
                         corfiltertype = corfiltertype,
                         pfilter = pfilter,
                         adjust = adjust,
                         cormethod = cormethod,
                         use = use,
                         xname = xname,
                         yname = yname,
                         ...)

      createdir(paste0(savepath,"/",samplename,"-",mode,"-select/"))

      write.table(x = corrdatamz$framedata$linkdata,
				file = paste0(savepath,"/",samplename,"-",mode,"-select/","corrplot.xls"),
				sep = "\t",
				row.names = F)


	  write.table(x = corrdatamz$filterdata$linkdata,
				file = paste0(savepath,"/",samplename,"-",mode,"-select/","corrnetwork.xls"),
				sep = "\t",
				row.names = F)

      filtermz <- fitmz[-which(fitmz == fig)]

      imzmlimage(filename = imzmlfile,
                 mapname = paste0("Intensity-",samplename,"-",mode),
                 savepath = paste0(savepath,"/",samplename,"-",mode,"-select/"),
                 mapmz = T,
                 filtermz = filtermz)

  }


}

#' @export
datatocorr<- function(x,
                      y = NULL,
                      ci = FALSE,
                      xname = "Featurex",
                      yname = "Featurey",
                      corfilter = 0.95,
                      corfiltertype = "+-",
                      pfilter = 0.05,
                      adjust = "none",
                      transx = T,
                      transy = T,
                      cormethod = "pearson",
                      use = "pairwise",
                      networkmapmoudle = corrnetwork,
                      networkmapname = "corrnetwork",
                      plotmapmoudle = corrplot,
                      plotmapname = "corrplot",
                      savepath = "./",
                      networkparams = list(algorithm = "nicely"),
                      plotparams = list()){

  corrdata <- corrcal(x = x,
                      y = y,
                      ci = ci,
                      xname = xname,
                      yname = yname,
                      corfilter = corfilter,
                      corfiltertype = corfiltertype,
                      pfilter = pfilter,
                      adjust = adjust,
                      transx = transx,
                      transy = transy,
                      method = cormethod,
                      use = use)

  if(is.null(corrdata)){
    write.csv(x = "The data is empty, no correlation analysis is performed",
              file = paste0(savepath,"/","desc.txt"),
              append = T)
    return(corrdata)
  }

  params <- list(corr = corrdata$rawcor$r,
                 p.mat = corrdata$rawcor$p.adj,
                 mapname = plotmapname,
                 savepath = savepath)
  params <- c(params,plotparams)
  plotdata <- do.call(what = plotmapmoudle,args = params)
  corrdata$plotdata <- plotdata

  if(dim(corrdata$filterdata$linkdata)[1] == 0){
    write(x = "After correlation>",corfilter,",significant<",pfilter,",no matching results",
          file = paste0(savepath,"/","desc.txt"),
          append = T)
    return(corrdata)
  }

  params <- list(linkdata = corrdata$filterdata$linkdata,
                 nodedata = corrdata$filterdata$nodedata,
                 mapname = networkmapname,
                 savepath = savepath)
  params <- c(params,networkparams)
  networkdata <- do.call(what = networkmapmoudle,args = params)
  corrdata$networkdata <- networkdata

  return(corrdata)

}