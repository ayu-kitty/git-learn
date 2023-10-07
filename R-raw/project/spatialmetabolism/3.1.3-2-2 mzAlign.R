
#' @export
mulmzAlign3 <- function(object,...){
  
  group <- as.character(run(object))
  group <- unique(group)
  i <- 1
  base::print("进入质量轴校正")
  base::print(group[i])
  object1 <- subsetPixels(object,run(object) == group[i])
  object1 <- mzAlign3(object = object1,...)
  
  if (length(group) > 1) {
    
    for ( i in 2:length(group)) {
      base::print(group[i])
      object2 <- subsetPixels(object,run(object) == group[i])
      object2 <- mzAlign3(object = object2,...)
      runtry <- try({
        object1 <- BiocGenerics::cbind(object1,object2)
        gc(reset = T)
      })
      if("try-error" %in% class(runtry)){
        base::print(mz(object2)[!(mz(object2) %in% mz(object1))])
        base::print(mz(object1)[!(mz(object1) %in% mz(object2))])
        stop("质量轴未对上")
      }
    }
    
  }
  
  return(object1)
}

#' 分段质量轴矫正
#'
#' @param obeject mse数据
#' @param number 每number个像素点拆分进行分段矫正,默认1000 
#' 
#' @export
mzAlign3 <- function(object,number = 1000,mapname = run(object)[1],savepath = "./",...){
  
  mse <- object
  # 分段数计算
  num <- floor(dim(mse)[2]/number)
  
  if(num == 0){
    num <- 1
  }
  
  mse3 <- NULL
  
  for ( i in 1:num) {
    # 开始像素点计算
    startnum <- number*(i-1)+1
    
    # 结束像素点计算
    endnum <- number*i
    
    # 如果结束像素点不是最后一个像素点,则赋值为最后一个像素点
    if(endnum+number > dim(mse)[2]){
      endnum <- dim(mse)[2]
    }
    
    base::print(paste0("第",i,"段质量轴校正"))
    # 分段mse
    mse2 <- mse[,startnum:endnum]
    # 进行分段质量轴矫正
    mse2 <- part_mzAlign3(object = mse2,mapname = paste0(mapname,"-",i),savepath = savepath,...)
    
    # 矫正后数据合并
    if(is.null(mse3)){
      mse3 <- mse2
    }else{
      runtry <- try({
        mse3 <- BiocGenerics::cbind(mse3,mse2)
        gc(reset = T)
      })
      if("try-error" %in% class(runtry)){
        # mse3 <<- mse3
        # mse2 <<- mse2
        base::print(mz(mse3)[!(mz(mse3) %in% mz(mse2))])
        base::print(mz(mse2)[!(mz(mse2) %in% mz(mse3))])
        stop("质量轴未对上")
      }
    }
  }
  
  return(mse3)
}

#' 分段质量轴矫正
#'
#' @param obeject mse数据
#' @param ref 归一化参考mz,一般是内标,例如554,556
#' @param refmz 归一化mz列表
#' @param tolerance 峰矫正参数
#' @param units 单位ppm or mz
#' @export
part_mzAlign3 <- function(object, ref, refmz = NULL,
                          tolerance = NA, units = c("ppm", "mz"),
                          maxintensity = 400,
                          mapname = run(object)[1],
                          savepath = "./"){
  
  tol <- switch("ppm",
                ppm = c("relative" = unname(200) * 1e-6),
                mz = c("absolute" = unname(0.1)))
  tol2 <- switch("ppm",
                 ppm = c("relative" = unname(100) * 1e-6),
                 mz = c("absolute" = unname(0.1)))
  tol.ref <- switch(names(tol),
                    relative = "key",
                    absolute = "none")
  
  mse2 <- object
  resolution(mse2) <- c(ppm=4)
  
  # 获取所有mz
  mz <- mz(mse2)
  # 参考mz获取
  mz.ref <- ref[order(ref)]
  
  # 获得mz.ref ± 0.2 mz
  newmz <- NULL
  for ( i in 1:length(mz.ref)) {
    diffmz <- 200/1000000*mz.ref[i]
    
    newmz2 <- mz[(mz<= (mz.ref[i]+diffmz))&(mz>=(mz.ref[i]-diffmz))]
    newmz <- c(newmz,newmz2)
  }
  newmz <- newmz[order(newmz)]
  newmz <- newmz[!duplicated(newmz)]
  
  # 参考newmz 的mse
  newmse <- subsetFeatures(x = mse2,mz = newmz)
  
  # extractpixel <- pixelApply(newmse,function(x){!any(x[(newmz < 556.5) & (newmz > 554)] > 2000)})
  # if(sum(extractpixel) > 10){
  #   base::print("删除背景点进行质量轴计算")
  #   newmse <- subsetPixels(newmse,extractpixel)
  # }
  
  x <- featureApply(newmse, function(x){unname(quantile(x,probs=c(0.95)))})
  # x <- featureApply(newmse, mean)
  # x <- featureApply(newmse, median)
  x[x < 50] <- 0
  
  for ( i in 1:length(mz.ref)) {
    diffmz <- 200/1000000*mz.ref[i]
    
    if(all(x[newmz >= (mz.ref[i]+diffmz) | newmz <= (mz.ref[i]-diffmz)] < maxintensity)){
      
    }else{
      x[(x < maxintensity) & (newmz >= (mz.ref[i]+diffmz) | newmz <= (mz.ref[i]-diffmz))] <- 0
    }
  }
  
  max.test <- matter::locmax(x,halfWindow = 10)
  mz.test <- newmz[max.test]
  # 二分查找找到目标离子(内标)
  i1 <- matter::bsearch(mz.ref[mz.ref <= 400], mz.test, tol=tol, tol.ref=tol.ref)
  i2 <- matter::bsearch(mz.ref[mz.ref > 400], mz.test, tol=tol2, tol.ref=tol.ref)
  i <- c(i1,i2)
  
  found <- !is.na(i)
  if ( sum(found) < 1 ) {
    base::print(paste0("实际mz:",mz.test))
    # base::print(paste0("mz:",paste(newmz,collapse = ",")))
    base::print(paste0("强度:",paste(x,collapse = ",")))
    stop("未找到内标离子")
  }
  
  # 目标离子确认
  mz.ref <- mz.ref[found]
  i <- i[found]
  mz.test <- mz.test[i]
  base::print(paste0("实际mz:",paste(format(mz.test, nsmall = 5, trim = T),collapse = ";")))
  base::print(paste0("参考mz:",paste(format(mz.ref, nsmall = 5, trim = T),collapse = ";")))
  
  mse2 <- object
  mz <- mz(mse2)
  
  # 实际mz与参考mz的差值计算
  if(length(mz.ref) == 1){
    diff <- mz.ref - mz.test
    dmz <- diff[1]/mz.test[1]*mz
  }else{
    diff <- mz.ref - mz.test
    relativediff <- diff/mz.test
    
    base::print(paste0("mz差异:",paste(format(diff, nsmall = 5, trim = T,digits = 5),collapse = ";")))
    base::print(paste0("ppm差异:",paste(format(relativediff*1e6, nsmall = 1, trim = T,digits = 1),collapse = ";")))
    
    diff <- c(diff[1]/mz.test[1]*mz[1], diff, diff[length(diff)]/mz.test[length(mz.test)]*mz[length(mz)])
    mz.test2 <- c(mz[1], mz.test, mz[length(mz)])
    relativediff <- c(relativediff[1],relativediff,relativediff[length(relativediff)])
    
    span <- 1
    control <- loess.control()
    shift <- suppressWarnings(loess(diff ~ mz.test2, span=span, control=control))
    dmz <- predict(shift, mz)
    # shift <- suppressWarnings(loess(relativediff ~ mz.test2, span=span, control=control))
    # rdmz <- predict(shift, mz)
    # dmz <- mz*rdmz
  }
  
  plotdata <- data.frame(mz = mz,dmz = dmz)
  p <- ggplot2::ggplot(data = plotdata,mapping = ggplot2::aes(x = mz,y = dmz))+
    ggplot2::geom_point()+ggplot2::theme_bw()+
    ggplot2::geom_point(color="red",inherit.aes = F,data = data.frame(diff = diff,mz = mz.test2),mapping = ggplot2::aes(x = mz,y = diff))
  ggplotsave(plot = p,savepath = savepath,mapname = mapname,imagetype = c("png"))

  if(!is.null(refmz)){
    # 归一化
    tol <- switch(units[1],
                  ppm = c("relative" = unname(tolerance) * 1e-6),
                  mz = c("absolute" = unname(tolerance)))
    tol.ref <- switch(names(tol),
                      relative = "key",
                      absolute = "none")
    foundmz <- matter::bsearch(refmz, mz+dmz, tol=tol, tol.ref=tol.ref)
    if(any(is.na(foundmz))){
      stop("未找到对应的mz")
    }
    refmz2 <- mz[foundmz]
    mse2 <- mse2 %>%
      peakBin(ref = refmz2,tolerance = tolerance, units = units,type= "height") %>%
      process()
    mz(mse2) <- refmz
  }else{
    mse2 <- Cardinal::pull(mse2, as.matrix = TRUE)
    # spectra(mse2) <- as.matrix(spectra(mse2))
    mz(mse2) <- mz(mse2)+dmz
    
    if(sum(mz < min(mz(mse2))) > 0){
      for ( i in sum(mz < min(mz(mse2))):1) {
        mse3 <- mse2[1,]
        mz(mse3) <- mz[i]
        mse2 <- BiocGenerics::rbind(mse3,mse2)
        pData(mse2) <- pData(mse2)[,!duplicated(colnames(pData(mse2)))]
      }
    }
    
    if(sum(mz > max(mz(mse2))) > 0){
      for ( i in sum(mz > max(mz(mse2))):1) {
        mse3 <- mse2[dim(mse2)[1],]
        mz(mse3) <- mz[length(mz)-i+1]
        mse2 <- BiocGenerics::rbind(mse2,mse3)
        pData(mse2) <- pData(mse2)[,!duplicated(colnames(pData(mse2)))]
      }
    }
    
    extractmz <- matter::bsearch(key = mz,values = mz(mse2),nearest = T)
    if(any(duplicated(extractmz))){
      
      nummz <- which(duplicated(extractmz))
      mse3 <- subsetFeatures(x = mse2,mz = mz[1:(nummz[1]-1)])
      mz(mse3) <- mz[1:(nummz[1]-1)]
      
      for ( i in 1:length(nummz)){
        # base::print(nummz[i])
        # base::print(mz[c((nummz[i])-1,(nummz[i]))])
        if(i == length(nummz)){
          mse4 <- subsetFeatures(x = mse2,mz = mz[(nummz[i]):length(mz)])
          mz(mse4) <- mz[(nummz[i]):length(mz)]
          mse3 <- BiocGenerics::rbind(mse3,mse4)
          pData(mse3) <- pData(mse3)[,!duplicated(colnames(pData(mse3)))]
        }else{
          mse4 <- subsetFeatures(x = mse2,mz = mz[(nummz[i]):(nummz[i+1]-1)])
          mz(mse4) <- mz[(nummz[i]):(nummz[i+1]-1)]
          mse3 <- BiocGenerics::rbind(mse3,mse4)
          pData(mse3) <- pData(mse3)[,!duplicated(colnames(pData(mse3)))]
        }
      }
      
      mse2 <- mse3
    }else{
      mse2 <- subsetFeatures(x = mse2,mz = mz)
      mz(mse2) <- mz
    }
    
  }

  return(mse2)
}
