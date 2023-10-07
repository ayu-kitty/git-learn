#!/opt/conda/bin/Rscript

#' LCMSMS出报告函数
#' @param cloud 
#' @param ... 
#' @export
lmm<-function(cloud="T",...){
  fils<-dir("./")
  if(!"table.xlsx" %in% fils){
    for(fil in fils){
      setwd(fil)
      ta<-readxlsx("table.xlsx")
      if(!is.na(grep("LM",ta[3,2])[1])){
        mode="LM"
      }else mode ="OE"
      createdir("report")
      system("cp -r result report/")
      result <- system(paste0("python /data/hstore3/public/propip/report/LCMS/Report.py -m ",mode," -c ",cloud," > project_report_stderr.log"))
      if(result == "stop"){
        print(paste0("~",fil,"报告生成失败，请核查文件！"))
        return()
      }else print(paste0("~",fil,"报告已出具！"))
      
      setwd("..")
    }
  }else{
    ta<-readxlsx("table.xlsx")
    if(!is.na(grep("LM",ta[3,2])[1])){
      mode="LM"
    }else mode ="OE"
    createdir("report")
    system("cp -r result report/")
    system(paste0("python /data/hstore3/public/propip/report/LCMS/Report.py -m ",mode," -c ",cloud," > project_report_stderr.log"))
  }

}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  parser$add_argument("-c","--cloud", default = "T",help = "报告是否云交付，默认为T")
  args <- parser$parse_args()
  lmmreport <- do.call(lmm,args = args)
}
