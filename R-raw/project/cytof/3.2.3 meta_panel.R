#!/opt/conda/bin/Rscript

#' 构建metadata文件
#' 
#' @param flowsetObject flowSet对象的变量名
#' @param file_name md的file_name列，默认从flowSet对象中提取
#' @param sample_id md的sample_id列，默认同file_name
#' @param condition md的condition列，默认Ref
#' @param patient_id md的patient_id列，默认Patient1...
#' @param savemd 逻辑，是否保存flowSet对象的metadata文件
#' @param savepath 同原始存放fcs文件的路径，保存metadata文件
#' @param mdname 保存的metadata文件名，默认metadata.xlsx，可指定
create_md <- function(flowsetObject,
                      file_name = NULL,
                      sample_id = NULL,
                      condition = NULL,
                      patient_id = NULL,
                      savepath = "rawdata/",
                      ...){
  if(file.exists(paste0(savepath,"/metadata.xlsx"))){
    print("metadata.xlsx have already exists")
    return()
  }else{
    if(is.null(file_name)){
      file_name <- flowsetObject@phenoData@data$name
    }
    
    if(is.null(sample_id)){
      sample_id <- c()
      for(i in 1:length(file_name)){
        sample_id_tmp <- strsplit(file_name[i], "*\\.fcs$")[[1]]
        sample_id <- c(sample_id, sample_id_tmp)                                # sample_id
      }
    }
    
    if(is.null(condition)){
      condition <- rep("Ref", length(file_name))                                # condition
    }
    
    if(is.null(patient_id)){
      patient_id <- paste0("patient", 1:length(file_name))                      # patient_id
    }
    
    # 保存metadata
    md <- data.frame(file_name,
                     sample_id,
                     condition,
                     patient_id)
    
    savexlsx1(data = md,
              filename = paste0(savepath,"/metadata.xlsx"),
              sheet = "metadata")
    
    return(md)
  }
}


#' create_panel
#' 
#' 构建panel文件
#' 
#' @param flowsetObject flowSet对象的变量名
#' @param fcs_colname panel的fcs_colname列，默认从flowSet对象中提取
#' @param antigen  panel的antigen列，默认从flowSet对象中提取
#' @param marker_class panel的condition列，默认type
#' @param savepanel 逻辑，是否保存flowSet对象的panel文件
#' @param savepath 保存panel的路径
#' @param panelname 保存的panel文件名,可指定
create_panel <- function(flowsetObject,
                         fcs_colname = NULL,
                         antigen = NULL,
                         marker_class = NULL,
                         savepanel = T,
                         savepath = "rawdata/",
                         ...){
  if(file.exists(paste0(savepath ,"/panel.xlsx"))){
    print("panel.xlsx have already exists")
    return()
  }else{
    if(is.null(fcs_colname) & is.null(antigen) & is.null(marker_class)){
      print("panel.xlsx create from guessPanel")
      # guessPanel
      gsPaneldata <- guessPanel(flowsetObject[[1]]) 
      if(length(which((gsPaneldata$`use_channel`) == 0)) != 0){
        gsPaneldata1.1 <- gsPaneldata[-which((gsPaneldata$`use_channel`) == 0), ]
        gsPaneldata1.2 <- gsPaneldata1.1[is.na(match(gsPaneldata1.1$fcs_colname, gsPaneldata1.1$desc0)), ]
      }else{
        gsPaneldata1.2 <- gsPaneldata
      }
      fcs_colname <- unclass(gsPaneldata1.2$fcs_colname)                        # fcs_colname
      
      if(all(c("DNA1", "DNA2") %in% gsPaneldata1.2$desc)){
        anti  <- unlist(unclass(gsPaneldata1.2$desc))
      }else if(all(c("DNA1", "DNA2") %in% gsPaneldata1.2$desc0)){
        anti  <- unlist(unclass(gsPaneldata1.2$desc0))}else{
          anti  <- unlist(unclass(gsPaneldata1.2$antigen))}
      
      marker_class <- rep("type", length(fcs_colname))                          # marker_class
      antigen <- c()                                                            # antigen
      for(i in 1:length(anti)){
        tmp_anti <- unlist(str_split(anti[i], "/"))[1]
        antigen <- c(antigen, tmp_anti)}
      # 保存metadata
      panel <- data.frame(fcs_colname, antigen , marker_class)
      
      if(length(grep("191",panel$fcs_colname) &  grep("193",panel$fcs_colname)) != 0){
        panel <- panel[-grep('191', panel$fcs_colname),]
        panel <- panel[-grep('193', panel$fcs_colname),]
      }else{panel <- panel}
      
      if( length(which(panel$antigen %in% c("DNA1", "DNA2"))) == 0 ){
        panel <- panel
      }else{
        panel <- panel[-which(panel$antigen %in% c("DNA1", "DNA2")), ]
      }
      if( length(which(panel$antigen %in% "Live")) == 0 ){
        panel <- panel
      }else{
        panel <- panel[-which(panel$antigen %in% "Live"), ]
      }
      
      if(length(which(tolower(panel$fcs_colname) == tolower(panel$antigen))) != 0){
        panel <- panel[-which(tolower(panel$fcs_colname) == tolower(panel$antigen)),]
      }else{
        panel <- panel
      }
      
      savexlsx1(data = panel,
                filename = paste0(savepath,"/panel.xlsx"),
                sheet = "panel")
      
      return(panel)
    }
  }
}


meta_panel <- function(flowsetObject,
                       ...){
  suppressMessages(library("openxlsx"))
  suppressMessages(library("stringr"))
  suppressMessages(library("CATALYST"))
  
  print("读取rds文件")
  fs_object <- readdata(flowsetObject)
  print("开始处理metadata.xlsx文件")
  create_md(flowsetObject=fs_object, ...)
  print("开始处理panel.xlsx文件")
  create_panel(flowsetObject=fs_object, ...)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-fo","--flowsetObject",
                      default = "./rawdata/flowSet_object.rds", 
                      help = "读取fcs的结果文件：**.rds")
  parser$add_argument("-sp","--savepath",
                      default = "rawdata/", 
                      help = "保存metadata和panel文件的路径，默认rawdata")
  
  # metadata(md)特有参数
  parser$add_argument("-fn","--file_name", default = NULL, help = "md的file_name列，默认从flowSet对象中提取")
  parser$add_argument("-si","--sample_id", default = NULL, help = "md的sample_id列，默认同file_name")
  parser$add_argument("-cn","--condition", default = NULL, help = "md的condition列，默认Ref")
  parser$add_argument("-pi","--patient_id", default = NULL, help = "md的patient_id列，默认Patient1...")
  # panel特有参数
  parser$add_argument("-fc","--fcs_colname", default = NULL, help = "panel的fcs_colname列，默认从flowSet对象中提取")
  parser$add_argument("-at","--antigen", default = NULL, help = "panel的antigen列，默认从flowSet对象中提取")
  parser$add_argument("-mc","--marker_class", default = NULL, help = "panel的condition列，默认type")
  
  args <- parser$parse_args()
  
  metaargs <- do.call(what = meta_panel, args = args)
}
## 注意：每个目录和文件必须严格存在！rawdata是固定的下机数据存放目录 flowSet_object.rds是固定的文件名