#' @export
createproject <- function(project_id) {
  if (!grepl(pattern = "^[0-9A-Za-z-]*$", x = project_id)) {
    stop("项目编号不符合规则")
  }
  
  linkpath <- projectpath(path = project_id)
  
  if(file.exists(linkpath)){
    return()
  }
  
  # setwd(projectpath())
  path <- paste0("/data/hstore4/lumingos/metalib/",format(Sys.time(), "%Y%m"))
  path <- paste0(path,"/",project_id)
  createdir(path)
  createdir(paste0(path, "/raw/实验报告"))
  createdir(paste0(path, "/raw/搜库数据"))
  createdir(paste0(path, "/raw/质谱数据/GCMS"))
  createdir(paste0(path, "/raw/质谱数据/LCMS/neg"))
  createdir(paste0(path, "/raw/质谱数据/LCMS/pos"))
  createdir(paste0(path, "/raw/质谱数据/mzml/GCMS"))
  createdir(paste0(path, "/raw/质谱数据/mzml/LCMS/neg"))
  createdir(paste0(path, "/raw/质谱数据/mzml/LCMS/pos"))
  createdir(paste0(path, "/raw/质谱数据/发送数据"))
  file.copy(databasepath(path = "script/samplecheck/AutoSampleCheck.exe"),to = path,overwrite = T)
  system(command = paste0("ln -s '",fs::path_rel(path = path,start = projectpath()),"' '",linkpath,"'"))
  
  createanalysis(project_id = project_id,analysis_id = paste0(project_id,"-aa"),arrangeana = F)
  
}

#' @export
createanalysis <- function(project_id,
                           analysis_id,
                           arrangeana = T) {
  
  # if (!grepl(pattern = "^[0-9A-Za-z-]*$", x = project_id)) {
  #   stop("项目编号不符合规则")
  # }
  # if (!grepl(pattern = "^[0-9A-Za-z-]*$", x = analysis_id) | !grepl(pattern = project_id, x = analysis_id)) {
  #   stop("分析编号不符合规则")
  # }
  
  setwd(projectpath())
  
  if (!dir.exists(project_id)) {
    stop("项目目录不存在，请先生成项目目录")
  }
  
  setwd(project_id)
  
  path <- paste0("/data/hstore4/project/",format(Sys.time(), "%Y%m"))
  analysispath <- paste0(path,"/",project_id,"/",analysis_id)
  linkpath <- paste0(getwd(),"/",analysis_id)
  linkpath2 <- paste0(getwd(),"/raw")
  linkpath3 <- paste0(getwd(),"/raw/实验报告")
  createdir(analysispath)
  if(!file.exists(linkpath)){
    system(command = paste0("ln -s '",analysispath,"' '",linkpath,"'"))
  }
  if(!file.exists(paste0(analysispath,"/raw"))){
    system(command = paste0("ln -s '",linkpath2,"' '",analysispath,"/raw'"))
  }
  if(!file.exists(paste0(analysispath,"/实验报告"))){
    system(command = paste0("ln -s '",linkpath3,"' '",analysispath,"/实验报告'"))
  }
  
  GetAnalystInfo(analysis_id = analysis_id,savepath = analysis_id)
  setwd(analysis_id)
  if(file.exists("分析确认单.xlsx")){
    file.copy("分析确认单.xlsx",to = "实验报告",overwrite = T)
    if(arrangeana){
      ArrangeInfo()
    }
  }
}

#' @export
cleanproject <- function(){
  allproject <- list.files(projectpath(),full.names = T,recursive = F)
  for( i in 1:length(allproject)){
    if(fs::is_link(allproject[i])){
      if(!file.exists(allproject[i])){
        base::print(allproject[i])
        unlink(allproject[i])
        next 
      }
    }
    allproject2 <- list.files(allproject[i],full.names = T,recursive = F)
    if(length(allproject2) == 0){
      next
    }
    for ( j in 1:length(allproject2)) {
      if(fs::is_link(allproject2[j])){
        if(!file.exists(allproject2[j])){
          base::print(allproject2[j])
          unlink(allproject2[j])
        }
      }
    }
  } 
}
