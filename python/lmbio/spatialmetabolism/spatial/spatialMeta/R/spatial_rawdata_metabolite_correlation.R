#!/opt/conda/bin/Rscript
#' @export
spatial_rawdata_metacor <- function(imzmlfile,
								   df = df,
                                   mode = ifelse(grepl(pattern = "-neg\\.",basename(imzmlfile)),"neg","pos"),
                                   samplename = gsub(pattern = paste0("-",mode,"\\..*$"),replacement = "",x = basename(imzmlfile)),
                                   selectmz = NULL,
                                   savepath = "./sample/fig/result/meta/",
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
  # dataframe spectradata
  spectradata <- df

  if(is.null(selectmz)){
    x <- spectradata
    y <- NULL
  }else{
    x <- spectradata
    x <- x[selectmz,,drop=F]
    y <- spectradata
    y <- y[row.names(y) != selectmz,,drop=F]
  }

  # 计算相关性 lmbio::corrcal
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

  linkdata <- corrdata$filterdata$linkdata
  linkdata2 <- linkdata
  linkdata2 <<- linkdata2

  if(dim(linkdata)[1] == 0){
    return(corrdata)
  }

  # -----------------------------plot total pic-------------------------------------
  if(is.null(selectmz)){
    fitmz <- c(linkdata[,1],linkdata[,2])
    fitmz <- unique(fitmz)
    x <- spectradata
    x <- x[fitmz,,drop=F]
    y <- NULL
    write.table(x = x,
				file = paste0(savepath,"/","rawdata-",samplename,"-",mode,".xls"),
				sep = "\t",
				row.names = T)
  }else{
    x <- spectradata
    x <- x[selectmz,,drop=F]
    fitmz <- linkdata[,2]
    fitmz <- unique(fitmz)
    y <- spectradata
    y <- y[fitmz,,drop=F]
    write.table(x = x,
				file = paste0(savepath,"/","rawdata-",samplename,"-",mode,"-select.xls"),
				sep = "\t",
				row.names = T)
    write.table(x = y,
				file = paste0(savepath,"/","rawdata-",samplename,"-",mode,"-y.xls"),
				sep = "\t",
				row.names = T)
  }

  corrdatamz <- datatocorr(x = x,
                           y = y,
                           networkmapname = paste0("corrnetwork-",samplename,"-",mode),
                           plotmapname = paste0("corrplot-",samplename,"-",mode),
                           savepath = savepath,
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
				file = paste0(savepath,"/","corrplot-",samplename,"-",mode,".xls"),
				sep = "\t",
				row.names = F)

  write.table(x = corrdatamz$filterdata$linkdata,
				file = paste0(savepath,"/","corrnetwork-",samplename,"-",mode,".xls"),
				sep = "\t",
				row.names = F)


  # export relations
  # ------------------------------------------------------------------
  if(is.null(selectmz)){
    while (dim(linkdata2)[1] > 0) {
    linknum <- 0
    linkdata3 <- linkdata2[linkdata2[,1] == linkdata2[1,1] | linkdata2[,2] == linkdata2[1,1],]
    while (dim(linkdata3)[1] != linknum) {
      linknum <- dim(linkdata3)[1]
      fitmz <- c(linkdata3[,1],linkdata3[,2])
	  
      fitmz <- unique(fitmz)

      linkdata3 <- linkdata2[linkdata2[,1] %in% fitmz | linkdata2[,2] %in% fitmz,]
    }
    x <- spectradata[fitmz,,drop=F]
    y <- NULL

    createdir(paste0(savepath,"/",samplename,"-",mode,"/",fitmz[1]))
    write.table(x = x,
              file = paste0(savepath,"/",samplename,"-",mode,"/",fitmz[1],"/rawdata.xls"),
              sep = "\t",
              row.names = T)

    corrdatamz <- datatocorr(x = x,
                             y = y,
                             networkmapname = "corrnetwork",
                             plotmapname = "corrplot",
                             savepath = paste0(savepath,"/",samplename,"-",mode,"/",fitmz[1],"/"),
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
              file = paste0(savepath,"/",samplename,"-",mode,"/",fitmz[1],"/","corrplot.xls"),
              sep = "\t",
              row.names = F)

    write.table(x = corrdatamz$filterdata$linkdata,
              file = paste0(savepath,"/",samplename,"-",mode,"/",fitmz[1],"/","corrnetwork.xls"),
              sep = "\t",
              row.names = F)

    imzmlimage(filename = imzmlfile,
               mapname = paste0("Intensity-",samplename,"-",mode),
               savepath = paste0(savepath,"/",samplename,"-",mode,"/",fitmz[1],"/"),
               mapmz = T,
               filtermz = fitmz)


    linkdata2 <- linkdata2[!(linkdata2[,1] %in% fitmz | linkdata2[,2] %in% fitmz),]
    }
  } else{
    # selectMz
    for (mz in selectmz) {
      linkdata3 <- linkdata2[linkdata2[,1] == mz | linkdata2[,2] == mz,]
      fitmz <- c(linkdata3[,1],linkdata3[,2])
      fitmz <- unique(fitmz)
      x <- spectradata[fitmz,,drop=F]
      y <- NULL
      createdir(paste0(savepath,"/",samplename,"-",mode,"/",mz))
      write.table(x = x,
				file = paste0(savepath,"/",samplename,"-",mode,"/",mz,"/","rawdata.xls"),
				sep = "\t",
				row.names = T)

      corrdatamz <- datatocorr(x = x,
                               y = y,
                               networkmapname = "corrnetwork",
                               plotmapname = "corrplot",
                               savepath = paste0(savepath,"/",samplename,"-",mode,"-select/",mz,"/"),
                               corfilter = corfilter,
                               corfiltertype = corfiltertype,
                               pfilter = pfilter,
                               adjust = adjust,
                               cormethod = cormethod,
                               use = use,
                               xname = xname,
                               yname = yname,
                               ...)
      createdir(paste0(savepath,"/",samplename,"-",mode,"-select/",mz))
      write.table(x = corrdatamz$framedata$linkdata,
				file = paste0(savepath,"/",samplename,"-",mode,"-select/",mz,"/","corrplot.xls"),
				sep = "\t",
				row.names = F)


	  write.table(x = corrdatamz$filterdata$linkdata,
				file = paste0(savepath,"/",samplename,"-",mode,"-select/",mz,"/","corrnetwork.xls"),
				sep = "\t",
				row.names = F)

      imzmlimage(filename = imzmlfile,
                 mapname = paste0("Intensity-",samplename,"-",mode),
                 savepath = paste0(savepath,"/",samplename,"-",mode,"-select/",mz,"/"),
                 mapmz = T,
                 filtermz = fitmz)
    }

  }

  return(corrdata)
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
	write.csv(x = "Data is empty, stop analysis!",
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
