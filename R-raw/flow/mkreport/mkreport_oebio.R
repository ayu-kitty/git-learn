#!/opt/conda/bin/Rscript

#' 代谢项目记录
#' @export
metalog <- function(savepath="./"){
  suppressMessages(library("dplyr"))
  prouser <- c("葛淑婷","付双","许以成","丁龙海","康倩倩","潘静","李阿雨")
  names(prouser) <- c("shuting.ge","shuang.fu","yicheng.xu","longhai.ding","qianqian.kang","jing.pan","ayu.li")
  if(file.exists(paste0(savepath,"/项目执行记录.txt"))){
    log <- readdata(paste0(savepath,"/项目执行记录.txt"),header=F)
  }else{
    log <- data.frame(id=c(1,2,3))
  }
  
  if(file.exists("../内部分析单.xlsx")){
    if("差异比较信息" %in% getsheetname("../内部分析单.xlsx")){
      comn <- readxlsx("../内部分析单.xlsx",sheet="差异比较信息") %>% nrow()
    }else comn=0
    oin <- readxlsx("../内部分析单.xlsx",sheet="分析基本信息")
    proid <- oin[oin[,1]=="项目编号",2]
    if(is.na(grep(proid,log[,1])[1])){
      file_user <- prouser[Sys.info()["user"]]  # 获取文件的所有者用户名
      complete_time <- file.info(dir(pattern = ".*_report.html"))$mtime %>% gsub(" .*","",.) %>% gsub("-","/",.)
      sn <- readxlsx("../内部分析单.xlsx",sheet="样本基本信息")
      samp<-nrow(filter(sn,分组!="QC"))
      samtype <- paste0(oin[oin[,1]=="样本物种",2],oin[oin[,1]=="样本类型",2])
      protype <- oin[oin[,1]=="项目类别",2]
      tmp <- c(proid,samtype,protype,file_user,samp,comn,complete_time,oin[oin[,1]=="映射物种",2])
      savetxt(data = paste(tmp,collapse = "\t"),
              filename = paste0(savepath,"/项目执行记录.txt"),append = T)
    }
  }
}

#' mkreport_oebio
#'
#' 生成分析报告
#'
#' @param data 项目信息
#' @param ...
#'
#' @export
mkreport_oebionew <- function(...) {
  
  load("raw.RData")
  mkreport_oebio(data = data,...)
  
}

#' mkreport_oebio
#'
#' 生成分析报告
#'
#' @param data 项目信息
#' @param ...
#'
#' @export
mkreport_oebio <- function(data,
                           ...) {
  UseMethod("mkreport_oebio")
}

#' mkreport_oebio.default
#'
#' 生成分析报告
#'
#' @param data 项目信息
#' @param mode 公司信息模板，默认LM
#' @param mode2 分析类型
#' @param type 项目类型
#' @param reportwd 项目路径
#' @param exreport 逻辑，是否拷贝试验报告
#' @param getoebio 逻辑，是否生成分析报告
#' @param check 逻辑，是否检查报告
#' @param checkname 检查后保存的检查文件名
#'
#' @export
mkreport_oebio.default <- function(data,
                                   mode = "auto",
                                   confirmfile = "内部分析单.xlsx",
                                   getoebio = T,
                                   CLOUD = F,
                                   check = T,
                                   checkname = "检查文件.xlsx",
                                   ...) {
  
  base::print(getwd())
  pacman::p_load(dplyr,stringr)
  oin<-readxlsx(confirmfile,sheet="分析基本信息")
  type<-oin[oin[,1]=="项目类别",2]
  type <- gsub(pattern = "非靶向代谢",replacement = "全谱代谢",x = type)
  type <- gsub(pattern = "-EMDB",replacement = "",x = type)
  type <- gsub(pattern = "-PMDB",replacement = "",x = type)
  reportwd <- oin[oin[,1]=="项目报告",2]
  if(mode=="auto"){
    ty<-str_extract_all(reportwd,pattern = 'LM|QD|OE|SG|HY')[[1]][1]
    if(is.na(ty)){
      mode<-"LM"
    }else mode<-ty
  }
  setwd(reportwd)
  
  if (getoebio) {
    if (oin[oin[,1]=="项目类型",2] == "精准靶向代谢"){
      oebio <- "/data/hstore4/database/mould/report/Oebio/TLCMS/TargetReport.py"
      system(paste("/opt/conda/bin/python3", oebio, '-m', mode ,'-cl',ifelse(CLOUD,"T","F"), sep=" "), intern = T) 
    }else{
      oebio <- paste0("/data/hstore4/database/mould/report/Oebio/",type,"/Report.py")
      system(paste("/opt/conda/bin/python3", oebio, '-m', mode ,'-cl',ifelse(CLOUD,"T","F"), sep=" "), intern = T)
      loctime=format(Sys.Date(),format="%Y%m")
      if(!is.na(dir("./",pattern = ".*html$")[1]) & CLOUD=="T"){
        metalog(savepath = "/data/hstore4/project/")
        htfile<-dir("./",pattern = ".*html$")[1]
        system(paste0("cp ",htfile," /data/nas/177/代谢/项目报告检查空间/预报告文件夹/",loctime))
        #NAS存储
        topath=gsub("/data/hstore4/project",paste0("/data/nas/173/代谢生信部/1.Projects/",format(Sys.Date(),format="%Y")),dirname(getwd()))
        createdir(topath)
        system(paste0("cp -r ../class* ../内部分析单*.xlsx ../分析确认单.xlsx ../*结题报告 ",topath))
      }
    }
  }
  else {
    stop(paste0(reportwd,"未生成分析报告！"))
  }
  
  setwd("../")
  
  # if (check) {
  #   checkreport(data = data, checkname = checkname)
  # }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  parser$add_argument("-mo","--mode",default = "auto", help = "公司模板,如LM;OE;QD等")
  parser$add_argument("-pl","--getoebio",default = T, help = "是否获取oebio报告")
  parser$add_argument("-cl","--CLOUD",default = F,help = "是否云交付")
  #parser$add_argument("-ck","--checkname",default = F,help = "checkout")
  
  args <- parser$parse_args()
  
  result <- do.call(what = mkreport_oebionew,args = args) 
  
}
