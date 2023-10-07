#!/opt/conda/bin/Rscript

#' 图片保存
#' 
#' @param savepath 保存路径
#' @param mapname 图片名称
#' @param imagetype 图片格式
#' @param height 图片长度
#' @param width 图片宽度
#' @param dpi 图片分辨率
#' @param family 字体
#'
#' @export 
plotfile <- function(savepath = "./",
                     mapname = "test",
                     imagetype = c("png", "pdf"),
                     height = 6,
                     width = 5,
                     dpi = 300,
                     family = "sans",
                     units = "in"){
  suppressMessages(library("pryr"))
  # graphics.off()
  
  height2 <- convert_to_inches(x = height,units = units)
  width2 <- convert_to_inches(x = width,units = units)
  
  if(units == "px"){dpi <- "auto"}
  
  Cairo::Cairo(type = "raster",
               units = units,
               width = width2, 
               height = height2,
               family = family,
               bg = "white",
               dpi = dpi)
  
  if(units == "px"){dpi <- NA}
  
  tmpfun <- make_function(args = list(savepath = savepath,
                                      mapname = mapname,
                                      imagetype = imagetype,
                                      height = height,
                                      width = width,
                                      compression = "zip",
                                      dpi = dpi,
                                      units = units,
                                      family = family), 
                          body = body(plotsave))
  
  if("package:lmbio" %in% search()){
    detach("package:lmbio")
    assignInNamespace(x = "plotsave",value = tmpfun,ns = asNamespace("lmbio"))
    suppressMessages(library("lmbio"))
  }else{
    assignInNamespace(x = "plotsave",value = tmpfun,ns = asNamespace("lmbio"))
  }
}


#' 图片保存
#' 
#' @param savepath 保存路径
#' @param mapname 图片名称
#' @param imagetype 图片格式
#' @param height 图片长度
#' @param width 图片宽度
#' @param dpi 图片分辨率
#' @param compression tiff格式压缩模式
#' @param family 字体
#'
#' @export 
plotsave <- function(savepath = "./",
                     mapname = "test",
                     imagetype = c("png", "pdf"),
                     height = 6,
                     width = 5,
                     compression = "zip",
                     units = "in",
                     dpi = ifelse(units == "px",NA,300),
                     family = "sans"){
  
  height2 <- convert_to_inches(x = height,units = units)
  width2 <- convert_to_inches(x = width,units = units)
  
  for (type in imagetype) {
    if (is.na(type)) {
      devopt <- getOption("device")
      
      if(is.function(devopt)){
        if(length(dev.list()) > 1){
          print(dev.cur())
          print(dev.list())
          
          dev.copy(which = dev.prev())
          dev.off(dev.list()[names(dev.list()) == "Cairo"][1])
          
          print(dev.list())
        }else{
          dev.copy(device = devopt)
        }
        
      }else{
        if("RStudioGD" %in% getOption("device")){
          num <- dev.list()[names(dev.list()) == "RStudioGD"]
          if(length(num) > 0){
            dev.off(num) 
          }
          dev.copy(device = RStudioGD)
        }
      }
      
      return()
    } else {
      filename <- paste0(savepath,"/",mapname,".",type)
      createdir(filename = savepath)
      # print(paste0(filename,"保存中"))
      if (type == "tiff") {
        dev.copy(device = tiff,
                 filename = filename,
                 width = width2, height = height2,
                 units = units, pointsize = 12, res = dpi,
                 family = family, compression = "zip")
        dev.off()
        Sys.chmod(paths = filename,mode = "0777",use_umask = F)
      } else if (type == "pdf") {
        dev.copy(device = pdf,
                 file = filename,
                 width = width, height = height, 
                 family = family)
        dev.off()
        Sys.chmod(paths = filename,mode = "0777",use_umask = F)
      } else if (type == "jpg") {
        dev.copy(device = jpeg,
                 filename = filename,
                 width = width2, height = height2, 
                 units = units,pointsize = 12, res = dpi, 
                 family = family)
        dev.off()
        Sys.chmod(paths = filename,mode = "0777",use_umask = F)
      } else if (type == "png") {
        dev.copy(device = png,
                 filename = filename,
                 width = width2, height = height2,
                 units = units, pointsize = 12, res = dpi, 
                 family = family)
        dev.off()
        Sys.chmod(paths = filename,mode = "0777",use_umask = F)
      }
      
    }
  }
  
  graphics.off() 
}

#' @export
convert_to_inches <- function(x, units) {
  x <- switch(units, 
              px = x*72,
              `in` = x, 
              cm = x/2.54, 
              mm = x/2.54/10)
  return(x)
}
