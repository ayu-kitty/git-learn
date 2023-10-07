#' @export
lmcal <- function(x,
                  y,
                  lmdatalist = list(lm(y~x),
                                    lm(y~x,weights = 1/x),
                                    lm(y~x,weights = 1/x^2),
                                    lm(y~x,weights = x),
                                    lm(y~x,weights = x^2)),
                  clrealrentionime,
                  c17realrentionime,
                  clrentionime,
                  c17rentionime,
                  cldiffrentionime,
                  c17diffrentionime,
                  errorstop = F){
  if(length(lmdatalist)==0){
    stop("校正效果不佳，请人工审核")
  }
  lmdata <- lmdatalist[[1]]

  clrentionime2 <- predict(lmdata,data.frame(x=c(clrealrentionime)))
  c17rentionime2 <- predict(lmdata,data.frame(x=c(c17realrentionime)))
  print(paste0("C17校正后保留时间：",c17rentionime2))
  print(paste0("2-氯苯丙氨酸校正后保留时间：",clrentionime2))
  
  if(!is.null(clrentionime) & !is.null(clrentionime)){
    
    if((abs(clrentionime2-clrentionime) < cldiffrentionime)&(abs(c17rentionime2-c17rentionime) < c17diffrentionime)){
      
    }else{
      print("校正效果不佳，更换模型运算")
      lmdata <- lmcal(x = x,
                      y = y,
                      lmdatalist = lmdatalist[-1],
                      clrealrentionime = clrealrentionime,
                      c17realrentionime = c17realrentionime,
                      clrentionime = clrentionime,
                      c17rentionime = c17rentionime,
                      cldiffrentionime = cldiffrentionime,
                      c17diffrentionime = c17diffrentionime)
    }
    
  }else if(!is.null(c17rentionime)){
    
    if((abs(c17rentionime2-c17rentionime) < c17diffrentionime)){
      
    }else{
      print("校正效果不佳，更换模型运算")
      lmdata <- lmcal(x = x,
                      y = y,
                      lmdatalist = lmdatalist[-1],
                      clrealrentionime = clrealrentionime,
                      c17realrentionime = c17realrentionime,
                      clrentionime = clrentionime,
                      c17rentionime = c17rentionime,
                      cldiffrentionime = cldiffrentionime,
                      c17diffrentionime = c17diffrentionime)
    }
    
  }else if(!is.null(clrentionime)){
    
    if((abs(clrentionime2-clrentionime) < cldiffrentionime)){
      
    }else{
      print("校正效果不佳，更换模型运算")
      lmdata <- lmcal(x = x,
                      y = y,
                      lmdatalist = lmdatalist[-1],
                      clrealrentionime = clrealrentionime,
                      c17realrentionime = c17realrentionime,
                      clrentionime = clrentionime,
                      c17rentionime = c17rentionime,
                      cldiffrentionime = cldiffrentionime,
                      c17diffrentionime = c17diffrentionime)
    }
  }
  
  return(lmdata)
}
