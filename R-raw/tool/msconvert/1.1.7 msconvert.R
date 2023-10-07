#!/opt/conda/bin/Rscript

#' msconvert
#'
#' 将质谱数据转格式为mzml
#'
#' @param mzwd 质谱数据路径，默认当前路径
#' @param path 输出路径,默认/public/mzdata
#' @param raw 质谱数据格式，默认wiff
#' @param intern 在R中运行，默认FALSE
#' @param ignore.stdout 逻辑,是否忽略输出，默认FALSE
#' @param ignore.stderr 逻辑，是否忽略输出，默认FALSE
#' @param wait 逻辑，是否等待运行完成，默认FALSE
#' @param ... 见[system()]
#'
#' @export
msconvert <- function(mzwd = "./",
                      path = "/public/mzdata",
                      raw = "wiff",
                      intern = FALSE,
                      ignore.stdout = FALSE,
                      ignore.stderr = FALSE,
                      wait = F,
                      ...) {
  wd <- getwd()
  
  convertfile <- list.files(path = mzwd,pattern = paste0("\\.",raw,"$"),full.names = T)
  
  suppressWarnings(library(foreach))
  suppressWarnings(library(doParallel))
  registerDoParallel(cores=5)
  result <- foreach(j=seq_len(length(convertfile))) %dopar% {
    
    system(paste0(packagepath(path = "command/tool_msconvert.sh"),
                  " -o ", path,
                  " --ignoreUnknownInstrumentError '",
                  convertfile[j],"'"),
           intern = intern,
           ignore.stdout = ignore.stdout,
           ignore.stderr = ignore.stderr,
           wait = wait,
           ...)
    
  }
  
  return()
}

wdmsconvert <- function(datawd = "./",
                        mzmlwd = "./",
                        wait = F){
  
  if (dir.exists(datawd)) {
    massfile <- c(list.files(path = datawd, pattern = "\\.D$"),
                  list.files(path = datawd, pattern = "\\.wiff$"),
                  list.files(path = datawd, pattern = "\\.raw$"))
    
    if (length(massfile) == 0) {
      warning(paste0("未发现 ", datawd, " 目录下质谱数据"), immediate. = T)
      
      mzmlfile <- list.files(path = mzmlwd, pattern = "\\.mzML$")
      
      if (length(mzmlfile) == 0) {
        warning(paste0("未发现 ", mzmlwd, " 目录下mzML数据"), immediate. = T)
      } else {
        
      }
    } else {
      datamode <- gsub(pattern = ".*\\.", replacement = "", massfile[1])
      mzmlfile <- list.files(path = mzmlwd, pattern = "\\.mzML$")
      if (length(mzmlfile) == 0) {
        warning(paste0(datawd, " 下质谱数据未转格式，现进行转格式流程,转格式结果将保存到 ", mzmlwd), immediate. = T)
        
        msconvert(mzwd = datawd,
                  path = mzmlwd,
                  raw = datamode,
                  wait = wait)
        
      } else {
        warning(paste0(datawd, " 下质谱数据有转格式"), immediate. = T)
      }
    }
  } else {
    warning(paste0("未发现 ", datawd, " 目录"), immediate. = T)
  }
}

#' mulmsconvert
#'
#' 根据项目进行质谱数据转格式为mzml
#'
#' @param project_id 项目编号
#' @param wait 逻辑，是否等待运行完成，默认FALSE
#'
#' @export
mulmsconvert <- function(project_id = "test",
                         wait = F) {
  if (!grepl(pattern = "^[0-9A-Za-z-]*$", x = project_id)) {
    stop("项目编号不符合规则")
  }
  
  setwd(projectpath())
  
  if (dir.exists(project_id)) {
    setwd(project_id)
    wdmsconvert(datawd = "raw/质谱数据/LCMS/pos",
                mzmlwd = "raw/质谱数据/mzml/LCMS/pos",
                wait = wait)
    wdmsconvert(datawd = "raw/质谱数据/LCMS/neg",
                mzmlwd = "raw/质谱数据/mzml/LCMS/neg",
                wait = wait)
    wdmsconvert(datawd = "raw/质谱数据/GCMS",
                mzmlwd = "raw/质谱数据/mzml/GCMS",
                wait = wait)
  } else {
    stop("项目目录不存在，请先生成项目目录")
  }
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-dw","--datawd",default = "./", help = "数据路径")
  parser$add_argument("-mw","--mzmlwd",default = "./", help = "输出路径")
  parser$add_argument("-w","--wait",default = F, help = "是否放后台运行",action='store_true')
  
  args <- parser$parse_args()
  
  result <- do.call(what = wdmsconvert,args = args) 
  
}
