#' 项目数据处理
#'
#' @param project_id 项目编号
#' @param analysis_id 分析编号
#' @param process 处理方式，现有pre（预处理）和ana（分析）
#'
#' @export
metanalystflow <- function(project_id,
                           analysis_id,
                           process = "pre",
                           mode = "LM") {
  setwd(projectpath())
  
  if (dir.exists(project_id)) {
    setwd(project_id)
  } else {
    stop(paste0("无", project_id, "目录"))
  }
  if (dir.exists(analysis_id)) {
    setwd(analysis_id)
  } else {
    stop(paste0("无", analysis_id, "目录"))
  }
  
  
  if (process == "pre") {
    commandR <- "source /etc/profile;source ~/.bashrc;metaboanalysis_2023 -t metabo_result_pre --report F --force T -ob F"
  } else if (process == "ana") {
    commandR <- paste0("source /etc/profile;source ~/.bashrc;metaboanalysis_2023 --mode ",mode)
  } else {
    stop("请选择正确分析步骤")
  }
  
  system(paste0(commandR, " >> log.txt 2>&1"),
         intern = F,
         ignore.stdout = F,
         ignore.stderr = F,
         wait = F)
}

#' @export
metafetch <- function(data = readregistration(),...) {
  # data <- readregistration()
  UseMethod("metafetch",data)
}

#' @export
metafetch.default <- function(data,
                              savepath = "raw",
                              ...) {
  if (!file.exists("raw.RData")) {
    data <- readregistration("内部分析单.xlsx")
    data <- autodataprocess(data,saminfo = "批次信息.xlsx")
    save(data, file = "raw.RData")
  } else {
    load(file = "raw.RData")
  }
  
  if(!file.exists(paste0(savepath,"/org.txt"))){
    species <- data[["info"]][["basic"]][["species"]]
    savetxt(data = species,filename = paste0(savepath,"/org.txt"))
  }
  if(!file.exists(paste0(savepath,"/project_title.txt"))){
    projecttitle <- data[["info"]][["basic"]][["项目报告"]]
    savetxt(data = projecttitle ,filename = paste0(savepath,"/project_title.txt"))
  }
  
  organizeobj(obj = data,...)
  
  returndata <- list(class = class(data),
                     pro_type = data$info$basic$项目类型,
                     compare = data$info$basic$比较分析,
                     qc = data$info$basic$QC质控)
  
  return(returndata)
}

#' @export
metafetch.UntargetBoth <- function(data,
                                   savepath = "raw",
                                   ...) {
  if (!file.exists("raw.RData")) {
    data <- readregistration("内部分析单.xlsx")
    lcdata <- readregistration("内部分析单-lc.xlsx")
    lcdata <- autodataprocess(lcdata,saminfo = "批次信息-lc.xlsx",adjustsavepath = "批次校正-lc")
    gcdata <- readregistration("内部分析单-gc.xlsx")
    gcdata <- autodataprocess(gcdata,saminfo = "批次信息-gc.xlsx",adjustsavepath = "批次校正-gc")
    
    save(lcdata, gcdata, file = "raw.RData")
    
    lcinfo <- lcdata$data$predata$information
    row.names(lcinfo) <- paste0(row.names(lcinfo),"~LCMS")
    colnames(lcinfo)[colnames(lcinfo)=="Fragmentation Score"] <- "Fragmentation Score/Total spectrum similarity"
    gcinfo <- gcdata$data$predata$information
    row.names(gcinfo) <- paste0(row.names(gcinfo),"~GCMS")
    colnames(gcinfo)[colnames(gcinfo)=="Total spectrum similarity"] <- "Fragmentation Score/Total spectrum similarity"
    lcinfo$`Ion mode` <- paste0("LCMS-",lcinfo$`Ion mode`)
    gcinfo$`Ion mode` <- "GCMS"
    for (namelist in colnames(gcinfo)[!(colnames(gcinfo) %in% colnames(lcinfo))]) {
      lcinfo[,namelist] <- NA
    }
    for (namelist in colnames(lcinfo)[!(colnames(lcinfo) %in% colnames(gcinfo))]) {
      gcinfo[,namelist] <- NA
    }
    
    information_data <- rbind(lcinfo,gcinfo)
    data$data$predata$information <- information_data
    
    lc_data <- lcdata$data$data$data
    row.names(lc_data) <- paste0(row.names(lc_data),"~LCMS")
    gc_data <- gcdata$data$data$data
    row.names(gc_data) <- paste0(row.names(gc_data),"~GCMS")
    lc_data <- lc_data[,colnames(lc_data) %in% colnames(gc_data)]
    gc_data <- gc_data[,colnames(gc_data) %in% colnames(lc_data)]
    all_data <- rbind(lc_data,gc_data)
    information_data <- information_data[row.names(information_data) %in% row.names(all_data),,drop = F]
    
    # LC的level1 优先保留 ,其次与GC合并的的保留GC
    df1 <- information_data[((information_data$`Ion mode` %in% c("LCMS-neg","LCMS-pos")) & information_data$level == "Level 1"),]
    df2 <- information_data[information_data$`Ion mode` == "GCMS",]
    df3 <- information_data[((information_data$`Ion mode` %in% c("LCMS-neg","LCMS-pos")) & information_data$level %in% c("Level 2","Level 3","Level 4")),]
    
    information_data2 <- rbind(df1,df2)
    information_data2 <- information_data2[!duplicated(information_data2$cid),,drop = F]
    information_data2 <- information_data2[!duplicated(tolower(information_data2$Metabolites)),,drop = F]
    
    information_data3 <- rbind(information_data2,df3)
    information_data3 <- information_data3[!duplicated(information_data3$cid),,drop = F]
    information_data3 <- information_data3[!duplicated(tolower(information_data3$Metabolites)),,drop = F]
    
    # 先按照打分排,再按照level排序
    information_data3 <- information_data3[order(information_data3$Score, decreasing = T),]
    information_data3 <- information_data3[order(information_data3$level, decreasing = F),]
    
    # information_data <- information_data[!duplicated(information_data$cid),,drop = F]
    # information_data <- information_data[!duplicated(tolower(information_data$Metabolites)),,drop = F]
    all_data <- all_data[row.names(all_data) %in% row.names(information_data3),,drop = F]
    data$data$data$data <- all_data
    lc_data <- lcdata$data$zeroprocess$data
    row.names(lc_data) <- paste0(row.names(lc_data),"~LCMS")
    gc_data <- gcdata$data$zeroprocess$data
    row.names(gc_data) <- paste0(row.names(gc_data),"~GCMS")
    lc_data <- lc_data[,colnames(lc_data) %in% colnames(gc_data)]
    gc_data <- gc_data[,colnames(gc_data) %in% colnames(lc_data)]
    all_data <- rbind(gc_data,lc_data)
    all_data <- all_data[row.names(all_data) %in% row.names(information_data3),,drop = F]
    data$data$zeroprocess$data <- all_data
    
    data$data$predata$information <- information_data3
    
    save(data, lcdata, gcdata, file = "raw.RData")
  } else {
    load(file = "raw.RData")
  }
  
  if(!file.exists(paste0(savepath,"/org.txt"))){
    species <- data[["info"]][["basic"]][["species"]]
    savetxt(data = species,filename = paste0(savepath,"/org.txt"))
  }
  if(!file.exists(paste0(savepath,"/project_title.txt"))){
    projecttitle <- data[["info"]][["basic"]][["项目报告"]]
    savetxt(data = projecttitle ,filename = paste0(savepath,"/project_title.txt"))
  }
  organizeobj(obj = data,...)
  
  returndata <- list(class = class(data),
                     pro_type = data$info$basic$项目类型,
                     compare = data$info$basic$比较分析,
                     qc = data$info$basic$QC质控)
  
  return(returndata)
}

#' @export
metafetch.UntargetBothEmdb <- metafetch.UntargetBoth

#' @export
metafetch.UntargetBothPmdb <- metafetch.UntargetBoth
