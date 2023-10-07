#' mkexreport_oebio
#'
#' 生成实验报告
#'
#' @param project_id 项目编号
#' @param analysis_id 分析编号
#' @param excopy 逻辑，是否生成试验报告填写模板
#' @param exmd 逻辑，是否生成实验报告模板
#'
#' @export
mkexreport_oebio <- function(project_id,
                             excopy = F,
                             exmd = T,
                             type = "非靶向代谢-LCMS",
                             mode = "LM") {
  setwd(projectpath())
  
  if (dir.exists(project_id)) {
    setwd(project_id)
  } else {
    stop(paste0("无", project_id, "目录"))
  }

  setwd("raw")
  data <- list()
  data$info$basic$项目类型 <- type
  if (excopy) {
    exfilecopyoebio(data)
  }
  if (exmd) {
    exreportoebio(data,mode = mode)
  }
}


exfilecopyoebio <- function(data,...) {
  UseMethod("exfilecopyoebio")
}


exfilecopyoebio.default <- function(data,
                                    type = data$info$basic$项目类型) {
  setwddir("实验报告")
  
  type <- gsub(pattern = "非靶向代谢",replacement = "全谱代谢",x = type)
  type <- gsub(pattern = "-EMDB",replacement = "",x = type)
  type <- gsub(pattern = "-PMDB",replacement = "",x = type)
  
  filepath <- dir(path = "/data/hstore4/database/exmould",
                  pattern = paste0(type, "-.*.xlsx"),
                  full.names = T)
  
  if (length(filepath) > 0) {
    base::print(getwd())
    file.copy(max(filepath), paste0(type, ".xlsx"), copy.mode = F)
    Sys.chmod(paths = paste0(type, ".xlsx"),mode = "0777",use_umask = F)
  } else {
    stop(paste0("不存在", type, "类型实验报告"))
  }
  
  setwd("../")
}



exfilecopyoebio.UntargetOut <- function(data) {
  stop("外来数据无实验报告")
}


exreportoebio <- function(data,...) {
  UseMethod("exreportoebio")
}


exreportoebio.default <- function(data,
                                  type = data$info$basic$项目类型,
                                  mode = "LM",
                                  report = T) {
  setwddir("实验报告")
  
  type <- gsub(pattern = "非靶向代谢",replacement = "全谱代谢",x = type)
  type <- gsub(pattern = "-EMDB",replacement = "",x = type)
  type <- gsub(pattern = "-PMDB",replacement = "",x = type)
  
  if (report) {
    base::print(getwd())
    oebio =paste0("/data/hstore4/database/mould/report/Oebio-实验/",type,"/","Report.py")
    system(paste("/opt/conda/bin/python3", oebio, '-m' , mode, sep=" "), intern = T)
  }
  else {
    stop(paste0("不生成","分析报告"))
  }
  
  setwd("../")
}


exreportoebio.UntargetOut <- function(data) {
  stop("外来数据无实验报告")
}


