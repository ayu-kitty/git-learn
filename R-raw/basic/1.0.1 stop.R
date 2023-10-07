#! /opt/conda/bin/Rscript

#' 报错优化
#'
#' @param ... 
#'
#' @export
stop <- function(...,call. = FALSE){
  base::stop("oeError: ",call. = call.,...)
}

#' 警告优化
#'
#' @param ... 
#'
#' @export
warning <- function(...,call. = FALSE,immediate. = T){
  base::warning("oeWarning: ",...,call. = FALSE,immediate. = immediate.)
}

#' @export
loginfodata <- NULL

#' 打印优化
#'
#' @param ... 
#'
#' @export
print <- function(x,...){
  if(is.character(x)){
    if(any(grepl(pattern = "^~",x = x))){
      # if(exists("loginfodata")){
      #   loginfodata <- get("loginfodata",pos = 1)
      # }else{
      #   loginfodata <- NULL
      # }
      # assign(x = "loginfodata",value = c(loginfodata,x),envir = .GlobalEnv)
      # loginfodata <<- c(loginfodata,x)
      x <- gsub(pattern = "^~",replacement = "",x = x)
      base::cat(paste0("oeInfo: ",x),sep = "\n",...)
    }else if(any(grepl(pattern = "^#",x = x))){
      # if(exists("loginfodata")){
      #   loginfodata <- get("loginfodata",pos = 1)
      # }else{
      #   loginfodata <- NULL
      # }
      # assign(x = "loginfodata",value = c(loginfodata,x),envir = .GlobalEnv)
      # loginfodata <<- c(loginfodata,x)
      x <- gsub(pattern = "^#",replacement = "",x = x)
      base::cat(x,sep = "\n",...)
    }else{
      # if(exists("loginfodata")){
      #   loginfodata <- get("loginfodata",pos = 1)
      # }else{
      #   loginfodata <- NULL
      # }
      # assign(x = "loginfodata",value = c(loginfodata,x),envir = .GlobalEnv)
      # loginfodata <<- c(loginfodata,x)
    }
  }else if(is.numeric(x)){
    x <- as.character(x)
    print(x,...)
  }else{
    base::print(x,...)
  }
}
