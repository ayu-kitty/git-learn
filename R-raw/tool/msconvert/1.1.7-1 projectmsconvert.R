#!/opt/conda/bin/Rscript

#' projectmsconvert
#'
#' 根据项目进行质谱数据转格式为mzml
#'
#' @param project_id 项目编号
#' @param wait 逻辑，是否等待运行完成，默认FALSE
#'
#' @export
projectmsconvert <- function(projectid = "test",
                             projectpath = "/public/lumingos/project/",
                             wait = F) {
  if (!grepl(pattern = "^[0-9A-Za-z-]*$", x = project_id)) {
    stop("项目编号不符合规则")
  }
  
  setwd(projectpath)
  
  if (dir.exists(projectid)) {
    setwd(projectid)
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
  
  parser$add_argument("-pi","--projectid",default = "test", help = "项目编号")
  parser$add_argument("-pp","--projectpath",default = "/public/lumingos/project/", help = "项目路径")
  parser$add_argument("-w","--wait",default = F, help = "是否放后台运行",action='store_true')
  
  args <- parser$parse_args()
  
  result <- do.call(what = projectmsconvert,args = args) 
  
}
