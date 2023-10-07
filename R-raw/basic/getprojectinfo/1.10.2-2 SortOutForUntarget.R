#!/opt/conda/bin/Rscript

#' GetRawDataInfo
#'
#' 获取raw，wiff数据路径
#'
#' @param sampleinfo 样本信息
#' @param analysisid 项目分析标号
#' @param listname 存储列名
#' @param datawd 数据路径
#' @param mzmlwd mzml路径
#'
#' @export
GetRawDataInfo <- function(sampleinfo,
                           listname,
                           datawd,
                           mzmlwd) {
  sample <- sampleinfo
  if (dir.exists(datawd)) {
    massfile <- c(
      list.files(path = datawd, pattern = "\\.D$"),
      list.files(path = datawd, pattern = "\\.wiff$"),
      list.files(path = datawd, pattern = "\\.raw$")
    )
    
    if (length(massfile) == 0) {
      warning(paste0("未发现 ", datawd, " 目录下质谱数据"), immediate. = T)
      
      mzmlfile <- list.files(path = mzmlwd, pattern = "\\.mzML$")
      
      if (length(mzmlfile) == 0) {
        warning(paste0("未发现 ", mzmlwd, " 目录下mzML数据"), immediate. = T)
      } else {
        mzmlname <- gsub(pattern = "\\.mzML", replacement = "", x = mzmlfile)
        for (i in order(nchar(sample$`实验名称`), decreasing = T)) {
          if (any(grepl(pattern = sample$`实验名称`[i], x = mzmlname))) {
            sample[i, paste0(listname, "-mzml")] <- paste0(mzmlwd, "/", mzmlfile[which(grepl(pattern = sample$`实验名称`[i], x = mzmlname))[1]])
            sample[i, paste0(listname, "-name")] <- mzmlname[which(grepl(pattern = sample$`实验名称`[i], x = mzmlname))[1]]
            mzmlfile <- mzmlfile[-which(grepl(pattern = sample$`实验名称`[i], x = mzmlname))[1]]
            mzmlname <- mzmlname[-which(grepl(pattern = sample$`实验名称`[i], x = mzmlname))[1]]
          } else {
            warning(paste0("在", mzmlwd, "中未发现", sample$`实验名称`[i], "样本"), immediate. = T)
          }
        }
      }
    } else {
      datamode <- gsub(pattern = ".*\\.", replacement = "", massfile[1])
      mzmlfile <- list.files(path = mzmlwd, pattern = "\\.mzML$")
      if (length(mzmlfile) == 0) {
        warning(paste0(datawd, " 下质谱数据未转格式，现进行转格式流程,转格式结果将保存到 ", mzmlwd), immediate. = T)
        
        msconvert(mzwd = datawd,
                  path = mzmlwd,
                  raw = datamode,
                  wait = T)
        
        mzmlfile <- list.files(path = mzmlwd, pattern = "\\.mzML$")
      }
      
      massname <- gsub(pattern = paste0("\\.", datamode), replacement = "", x = massfile)
      mzmlname <- gsub(pattern = "\\.mzML", replacement = "", x = mzmlfile)
      
      for (i in order(nchar(sample$`实验名称`), decreasing = T)) {
        if (any(grepl(pattern = sample$`实验名称`[i], x = mzmlname))) {
          sample[i, paste0(listname, "-mzml")] <- paste0(mzmlwd, "/", mzmlfile[which(grepl(pattern = sample$`实验名称`[i], x = mzmlname))[1]])
          sample[i, paste0(listname, "-name")] <- mzmlname[which(grepl(pattern = sample$`实验名称`[i], x = mzmlname))[1]]
          mzmlfile <- mzmlfile[-which(grepl(pattern = sample$`实验名称`[i], x = mzmlname))[1]]
          mzmlname <- mzmlname[-which(grepl(pattern = sample$`实验名称`[i], x = mzmlname))[1]]
        } else {
          warning(paste0("在", mzmlwd, "中未发现", sample$`实验名称`[i], "样本"), immediate. = T)
        }
        
        if (any(grepl(pattern = sample$`实验名称`[i], x = massname))) {
          sample[i, paste0(listname, "-raw")] <- paste0(datawd, "/", massfile[which(grepl(pattern = sample$`实验名称`[i], x = massname))[1]])
          sample[i, paste0(listname, "-name")] <- massname[which(grepl(pattern = sample$`实验名称`[i], x = massname))[1]]
          massfile <- massfile[-which(grepl(pattern = sample$`实验名称`[i], x = massname))[1]]
          massname <- massname[-which(grepl(pattern = sample$`实验名称`[i], x = massname))[1]]
        } else {
          warning(paste0("在", datawd, "中未发现", sample$`实验名称`[i], "样本"), immediate. = T)
        }
      }
    }
  } else {
    warning(paste0("未发现 ", datawd, " 目录"), immediate. = T)
  }
  
  return(sample)
}

#' ArrangeInfo
#'
#' 将内部分析单整理为分析确认单
#'
#' @param type 类型
#' @param xlsxname 读取数据表格，默认为内部分析单.xlsx
#'
#' @export
ArrangeInfo <- function(type = NA,
                        fromname = "分析确认单.xlsx",
                        xlsxname = "内部分析单.xlsx",
                        overwrite = F,
                        keggspecies = NA) {
  args <- getFunc_Paras()
  
  print("分析确认单整理中")
  
  #自动删除 分析确认单.xlsx里的 质控样本信息表
  if("质控样本信息表" %in% readxl::excel_sheets("分析确认单.xlsx")){
    wb1 <- openxlsx::loadWorkbook("分析确认单.xlsx")
    openxlsx::removeWorksheet(wb1,sheet="质控样本信息表")
    openxlsx::saveWorkbook(wb1, "分析确认单.xlsx",overwrite=TRUE)
  }
  
  wb <- openxlsx::createWorkbook()
  
  if (file.exists(fromname)) {
    
    sheetname <- getsheetname(fromname)
    
    if ("分析基本信息" %in% sheetname) {
      
      print("~~分析基本信息整理中")
      info <- readdata(filename = fromname, sheet = "分析基本信息")
      if(length(info[info$key == "映射物种", "value"]) == 0){
        info["映射物种",] <- c("映射物种",keggspecies)
      }
      if(is.na(keggspecies)){
        if(is.na(info[info$key == "映射物种", "value"])){
          keggspecies <- getkeggspecies(species = info[info$key == "样本物种", "value"])
        }else{
          keggspecies <- info[info$key == "映射物种", "value"]
        }
      }
      if(!is.na(keggspecies)){
        info <- info[info$key != "映射物种", ]
        info["映射物种",] <- c("映射物种",keggspecies)
        addkeggspecies(species = info[info$key == "样本物种", "value"],keggspecies = keggspecies)
      }else{
        stop(paste0("在",databasepath(path = "database/物种映射/map.xlsx"),"文件中添加物种信息"))
      }
      if (is.na(type)) {
        type <- if (length(info[info$key == "项目类别", "value"]) == 0) {
          info[info$key == "项目类型", "value"]
        } else {
          if(info[info$key == "项目类型", "value"]=="精准靶向代谢"){
            info[info$key == "项目类型", "value"]
          }else{
            info[info$key == "项目类别", "value"]
          }
        }
      } else {
        info[info$key == "项目类别","value"] <- type
      }
      if(info[info$key == "项目类型", "value"]=="精准靶向代谢"){
        type2 <- paste0(type,"-",info[info$key == "项目类别","value"])
      }else{
        type2 <- type
      }
      type2 <- gsub(pattern = "非靶向代谢",replacement = "全谱代谢",x = type2)
      type2 <- gsub(pattern = "-EMDB",replacement = "",x = type2)
      type2 <- gsub(pattern = "-PMDB",replacement = "",x = type2)
      if(is.na(grep("-b",info[info$key=="项目编号",2])[1])){
        info[info$key=="任务单号",2] <- paste0(info[info$key=="项目编号",2],"-b1")
      }else{
        info[info$key=="任务单号",2] <- info[info$key=="项目编号",2]
      }
      if(is.na(info[info[, 1] == "客户名称", 2]) | is.na(info[info[, 1] == "联系人", 2])){
        kh<-na.omit(c(info[info[, 1] == "客户名称", 2],info[info[, 1] == "联系人", 2]))[1]
        names(kh)<-"联系人"
        info["项目报告", ] <- c(
          "项目报告",
          paste0(
            info[info[, 1] == "项目编号", 2],
            "-", kh,
            "-", type2, "结题报告"
          )
        )
      }else if (info[info[, 1] == "客户名称", 2] == info[info[, 1] == "联系人", 2]){
        kh<-info[info[, 1] == "客户名称", 2]
        names(kh)<-"联系人"
        info["项目报告", ] <- c(
          "项目报告",
          paste0(
            info[info[, 1] == "项目编号", 2],
            "-", kh,
            "-", type2, "结题报告"
          )
        )
      }else{
        kh<-c(info[info[, 1] == "客户名称", 2],info[info[, 1] == "联系人", 2])
        names(kh)<-c("客户名称","联系人")
        info["项目报告", ] <- c(
          "项目报告",
          paste0(
            info[info[, 1] == "项目编号", 2],
            "-", info[info[, 1] == "客户名称", 2],
            "-", info[info[, 1] == "联系人", 2],
            "-", type2, "结题报告"
          )
        )
      }
      if(xlsxname == "内部分析单.xlsx"){
        init<-data.frame("目录"=c("项目类型","任务单号","项目编号",names(kh),"物种及样本","映射物种"),
                         "信息"=c(
                           type2,
                           info[info[, 1] == "任务单号", 2],
                           info[info[, 1] == "项目编号", 2],
                           kh,
                           ifelse(is.na(grep(info[info[, 1] == "样本物种", 2],info[info[, 1] == "样本类型", 2])[1]),paste0(info[info[, 1] == "样本物种", 2],info[info[, 1] == "样本类型", 2]),info[info[, 1] == "样本类型", 2]),
                           info["映射物种", 2]
                         ))
        
        savetxt(init,"init.txt",quote=F)
      }
      
      parameter <- getmysqldata(dbname = "parameter",
                                table = "lcparameter",
                                wherename = "type",
                                wheredata = type)
      
      if(dim(parameter)[1] == 0){
        
        stop(paste0("现不支持", type, "项目类型,请增加此项目类型"))
        
      }else{
        
        info["S3", ] <- c("S3", parameter[1, "S3"])
        info["索引", ] <- c("索引", "")
        info["ID", ] <- c("ID", parameter[1, "ID"])
        info["Metabolites", ] <- c("Metabolites", parameter[1, "Metabolites"])
        info["Compound ID", ] <- c("Compound ID", parameter[1, "Compound ID"])
        info["Ion mode", ] <- c("Ion mode", parameter[1, "Ion mode"])
        info["KEGG", ] <- c("KEGG", parameter[1, "KEGG"])
        info["预处理", ] <- c("预处理", "")
        info["missvalue", ] <- c("missvalue", parameter[1, "missvalue"])
        info["rsd", ] <- c("rsd", parameter[1, "rsd"])
        info["zeroprocess", ] <- c("zeroprocess", parameter[1, "zeroprocess"])
        info["samplenormalization", ] <- c("samplenormalization", parameter[1, "samplenormalization"])
        info["datatransformation", ] <- c("datatransformation", parameter[1, "datatransformation"])
        info["datascaling", ] <- c("datascaling", parameter[1, "datascaling"])
        info["分析参数", ] <- c("分析参数", "")
        info["log10L", ] <- c("log10L", parameter[1, "log10L"])
        info["PCAscaleC", ] <- c("PCAscaleC", parameter[1, "PCAscaleC"])
        
        if (!dir.exists("./raw/")) {
          if(dir.exists("../raw/")){
            system(command = "ln -s ../raw raw")
          }else{
            stop("raw目录不存在，请确认")
          }
        }
        
        print("~~搜库结果匹配中")
        info["原始文件", ] <- c("原始文件", "")
        info["LC原始文件", ] <- c("LC原始文件", "")

        info["negID", ] <- c("negID",
                             list.files(path = "./raw/搜库数据", pattern = "NEG-ID.csv$", full.names = T,recursive = T)[1])
        if(is.na(info["negID",2])){
          info["negID", ] <- c("negID",
                               list.files(path = "./raw/搜库数据", pattern = "NEG.((txt)|(xlsx))$", full.names = T,recursive = T)[1])
        }
        info["posID", ] <- c("posID",
                             list.files(path = "./raw/搜库数据", pattern = "POS-ID.csv$", full.names = T,recursive = T)[1])
        if(is.na(info["posID",2])){
          info["posID", ] <- c("posID",
                               list.files(path = "./raw/搜库数据", pattern = "POS.((txt)|(xlsx))$", full.names = T,recursive = T)[1])
        }
        info["negM", ] <- c("negM",
                            list.files(path = "./raw/搜库数据", pattern = "NEG-M.csv$", full.names = T,recursive = T)[1])
        info["posM", ] <- c("posM",
                            list.files(path = "./raw/搜库数据", pattern = "POS-M.csv$", full.names = T,recursive = T)[1])
        
        if (any(!is.na(info[c("negID", "posID", "negM", "posM"), 2]))) {
          info["原始文件", ] <- c("原始文件", "有")
          if (all(is.na(info[c("negID", "posID", "negM", "posM"), 2]))) {
            warning("LCMS匹配到部分搜库结果，请核对搜库结果命名是否准确", immediate. = T)
            print(info[c("negID", "posID", "negM", "posM"), ])
          }
        }
        
        # 判断是否EMDB和PMDB
        if(any(grepl("EMDB-",info[c("negID", "posID", "negM", "posM"),2])) & (type == "非靶向代谢-LCMS" | type == "非靶向代谢-双平台")){
          newargs <- args
          newargs$type <- paste0(type,"-EMDB")
          do.call(what = ArrangeInfo,args = newargs)
          return()
        }
        if(any(grepl("PMDB-",info[c("negID", "posID", "negM", "posM"),2])) & (type == "非靶向代谢-LCMS" | type == "非靶向代谢-双平台")){
          newargs <- args
          newargs$type <- paste0(type,"-PMDB")
          do.call(what = ArrangeInfo,args = newargs)
          return()
        }
        
        info["GC原始文件", ] <- c("GC原始文件",
                              list.files(path = "./raw/搜库数据", pattern = "^Area_.*", full.names = T,recursive = T)[1])
        info["其他原始文件", ] <- c("其他原始文件", "")
        info["数据矩阵", ] <- c("数据矩阵",
                            list.files(path = "./raw/搜库数据", pattern = "数据矩阵.xlsx$", full.names = T,recursive = T)[1])
        if (any(!is.na(info[c("数据矩阵"), 2]))) {
          info["原始文件", ] <- c("原始文件", "有")
        }
        
      }
      
      if ("样本信息" %in% sheetname & "样本分组信息" %in% sheetname) {
        print("~~样本信息及分组信息整理中")
        
        sample <- readdata(filename = fromname, sheet = "样本信息")
        if("批次" %in% colnames(sample)){
          sample <- sample[,colnames(sample) != "批次",drop = F]
        }
        
        
        if (dim(sample)[1] == 0) {
          stop("无样本信息")
        } else if (dim(sample)[1] < 10) {
          warning("~~~~样本数小于10，请确认", immediate. = T)
        }
        
        # addsheet1(data = sample[, 1:2],wb = wb,sheet = "样本信息")
        
        if("样本分组信息" %in% getsheetname(fromname)){
          samplegroup <- readdata(filename = fromname, sheet = "样本分组信息")
        }else{
          stop("~无样本分组信息")
        }
        
        if (dim(samplegroup)[1] == 0 ) {
          stop("~无样本分组信息")
        } else if (dim(samplegroup)[1] == 0) {
          warning("~~~~样本分组除QC外仅一组分组，请确认", immediate. = T)
        }
        
        samplegroup1 <- samplegroup[samplegroup$type == "常规", ]
        samplegroup2 <- samplegroup[samplegroup$type == "扩展", ]
        
        for (i in 1:dim(samplegroup1)[1]) {
          sample[sample$`样本分析名称` %in% unlist(strsplit(samplegroup1[i, 2], split = ",")), "分组"] <- samplegroup1[i, 1]
        }
        
        # 查找QC
        qcnameneg <- list.files(path = "./raw/质谱数据/LCMS/neg",pattern = "QC")
        qcnameneg <- basename(qcnameneg)
        qcnameneg <- gsub(pattern = "\\..*",replacement = "",x = qcnameneg)
        qcnameneg <- gsub(pattern = "-NEG",replacement = "",x = qcnameneg)
        qcnameneg <- gsub(pattern = "-neg",replacement = "",x = qcnameneg)
        qcnamepos <- list.files(path = "./raw/质谱数据/LCMS/pos",pattern = "QC")
        qcnamepos <- basename(qcnamepos)
        qcnamepos <- gsub(pattern = "\\..*",replacement = "",x = qcnamepos)
        qcnamepos <- gsub(pattern = "-POS",replacement = "",x =  qcnamepos)
        qcnamepos <- gsub(pattern = "-pos",replacement = "",x =  qcnamepos)
        qcnamegc <- list.files(path = "./raw/质谱数据/GCMS/",pattern = "QC",include.dirs = T)
        qcnamegc <- basename(qcnamegc)
        qcnamegc <- gsub(pattern = "\\..*",replacement = "",x = qcnamegc)
        qcname <- c()
        if(length(qcnameneg) > 0){
          qcname <- c(qcname,qcnameneg)
        }
        if(length(qcnamepos) > 0){
          qcname <- c(qcname,qcnamepos)
        }
        if(length(qcnamegc) > 0){
          qcname <- c(qcname,qcnamegc)
        }
        qcname <- qcname[!duplicated(qcname)]
        qcname <- gsub("[0-9]+.*-QC","QC",qcname)
        print("~~质控样本信息整理中")
        if (length(qcname) == 0) {
          warning("~~~~未发现质控样本信息", immediate. = T)
          info["QC", ] <- c("QC质控", "无")
          info["rsd", ] <- c("rsd", NA)
        } else {
          
          qcname <- qcname[order(qcname)]
          info["QC", ] <- c("QC质控", "有")
          qcdata <- data.frame("样本下机名称" = qcname)
          qcdata$`样本分析名称` <- qcdata$`样本下机名称`
          qcdata$`分组` <- "QC"
          
          samplegroup3 <- samplegroup[0, ]
          for (qcclass in unique(qcdata$`分组`)) {
            samplegroup3[qcclass, ] <- c(qcclass, paste(qcdata[qcdata$`分组` == qcclass, "样本分析名称"], collapse = ","), "常规")
          }
          samplegroup <- rbind(samplegroup3, samplegroup)
          
          sample <- rbind(qcdata,sample)
        }
        
        if(length(list.files(path = "./raw/",pattern = "^上机对照表.xlsx$")) > 0 ){
          experiment <- readdata(filename = "./raw/上机对照表.xlsx",sheet = 1)
          sample <- merge(x = sample,y = experiment,by = "样本分析名称",sort = T,all.x = T)
          sample[is.na(sample[, "实验名称"]),"实验名称"] <- sample[is.na(sample[, "实验名称"]),"样本分析名称"]
        }else{
          sample[, "实验名称"] <- sample[, "样本分析名称"]
        }
        
        if(!("批次" %in% colnames(sample))){sample[, "批次"] <- 1}
        sample[is.na(sample[, "批次"]), "批次"] <- 1
        
        sample <- sample[, c("样本下机名称", "实验名称","样本分析名称", "批次", "分组")]
        if(file.exists("./实验报告/拟靶向脂质.xlsx")){
          weiv<-readdata("./实验报告/拟靶向脂质.xlsx",sheet="称重")
          names(weiv)<-c("样本分析名称","称重/体积")
          sample<-dplyr::left_join(sample,weiv,by="样本分析名称")
          sample$`称重/体积`[is.na(sample$`称重/体积`)]<-round(mean(sample$`称重/体积`,na.rm = T),2)
        }else{
          if(!("称重/体积" %in% colnames(sample))){sample[, "称重/体积"] <- NA}
        }
        if(!("处理" %in% colnames(sample))){sample[, "处理"] <- NA}
        if ("差异比较信息" %in% sheetname) {
          print("~~差异比较信息整理中")
          comparegroup <- readdata(filename = fromname, sheet = "差异比较信息")
          
          
          if (dim(comparegroup)[1] == 0) {
            warning("~~~~未发现差异比较信息", immediate. = T)
            info["比较分析", ] <- c("比较分析", "无")
          } else {
            info["比较分析", ] <- c("比较分析", "有")
            
            for (i in 1:dim(comparegroup)[1]) {
              if (is.na(comparegroup[i, "comment"])) {
                comparegroup[i, "compare"] <- paste(comparegroup[i, "case"], comparegroup[i, "control"], sep = "/")
              } else if (comparegroup[i, "comment"] == "") {
                comparegroup[i, "compare"] <- paste(comparegroup[i, "case"], comparegroup[i, "control"], sep = "/")
              } else {
                if (all(unlist(strsplit(comparegroup[i, "comment"], split = "/")) %in% samplegroup$group)) {
                  comparegroup[i, "compare"] <- paste0(unique(c(
                    comparegroup[i, "case"],
                    comparegroup[i, "control"],
                    unlist(strsplit(comparegroup[i, "comment"], split = "/"))
                  )),
                  collapse = "/"
                  )
                } else {
                  comparegroup[i, "compare"] <- paste(comparegroup[i, "case"], comparegroup[i, "control"], sep = "/")
                }
              }
            }
            
            comparegroup <- comparegroup[!duplicated(comparegroup[, c("paired", "compare")]), ]
            
            for (i in 1:dim(comparegroup)[1]) {
              comparegroup[i, "num"] <- length(unlist(strsplit(comparegroup[i, "compare"], split = "/")))
            }
            
            j <- 1
            
            for (i in 1:dim(comparegroup)[1]) {
              if (any(unlist(strsplit(comparegroup[i, "compare"], split = "/")) %in% samplegroup2$group)) {
                comparegroup[i, "class"] <- paste0("ExtendGroup", j)
                sample[, paste0("ExtendGroup", j)] <- "Unassigned"
                
                
                if (any(duplicated(unlist(strsplit(samplegroup[samplegroup$group %in% unlist(strsplit(comparegroup[i, "compare"], split = "/")), "samples"],
                                                   split = ","
                ))))) {
                  stop(paste0(comparegroup[i, "compare"], "比较组中有重复样本"))
                } else {
                  for (k in unlist(strsplit(comparegroup[i, "compare"], split = "/"))) {
                    sample[sample$`样本分析名称` %in% unlist(strsplit(samplegroup[samplegroup$group %in% k, "samples"], split = ",")), paste0("ExtendGroup", j)] <- k
                  }
                }
                
                j <- j + 1
              } else {
                comparegroup[i, "class"] <- "Group"
              }
            }
            
            comparegroup[, "log10L"] <- parameter[1, "log10L"]
            comparegroup[, "PCAscaleC"] <- parameter[1, "PCAscaleC"]
            comparegroup[, "PLSscaleC"] <- parameter[1, "PLSscaleC"]
            comparegroup[, "OPLSscaleC"] <- parameter[1, "OPLSscaleC"]
            comparegroup[, "PLSpermI"] <- parameter[1, "PLSpermI"]
            comparegroup[, "OPLSpermI"] <- parameter[1, "OPLSpermI"]
            comparegroup[comparegroup$num > 2, "PLSpermI"] <- 200
            comparegroup[comparegroup$num > 2, "OPLSpermI"] <- 0
            comparegroup[, "UnivariateAnalysis"] <- "ttest"
            comparegroup[, "p.adjust.method"] <- "BH"
            
            comparegroup[, "VIP"] <- 1
            comparegroup[, "FC"] <- NA
            comparegroup[, "Pvalue"] <- 0.05
            comparegroup[, "adjPvalue"] <- NA
            comparegroup[, "errorFC"] <- 2
            comparegroup[, "order"] <- "VIP"
            
            addsheet1(data = comparegroup,wb = wb,sheet = "差异比较信息")
          }
        } else {
          warning("~~~~未发现差异比较信息", immediate. = T)
          info["比较分析", ] <- c("比较分析", "无")
        }
        
        if (dim(samplegroup)[1] <= 300) {
          samplegroup[, "shape"] <- rep(c(21, 22, 24, 23, 25), 60)[1:dim(samplegroup)[1]]
          # shape<-c("圆","方","上三","菱","下三")
          fill <- c("green4", "blue3", "firebrick",
                            "gold", "darkviolet", "darkorange",
                            "skyblue3", "olivedrab3", "dodgerblue3",
                            "aquamarine2", "deeppink3", "slateblue3",
                            "brown2", "palegreen2", "chocolate2",
                            "antiquewhite3", "steelblue1", "violetred1",
                            "burlywood3", "pink1", "slategray2",
                            "orangered1", "cyan3", "yellow4",
                            "red", "plum", "greenyellow",
                            "mediumpurple2", "tan1", "magenta")
                            
          if (dim(samplegroup)[1] <= 30) {
            fill <- fill[1:dim(samplegroup)[1]]
          } else {
            set.seed(111)
            colorrange <- colors()[c(1:150, 361:650)]
            colorrange <- colorrange[!(colorrange %in% fill)]
            fill <- c(fill, sample(colorrange, dim(samplegroup)[1] - length(fill)))
          }
          
          samplegroup[, "fill"] <- fill
          samplegroup[, "colour"] <- rep("dimgrey", 300)[1:dim(samplegroup)[1]]
        } else {
          stop("样本分组最多支持300组，现超过限制，请减少样本分组")
        }
        
        
        sample[, "GCMS-raw"] <- NA
        sample[, "GCMS-mzml"] <- NA
        sample[, "GCMS-name"] <- NA
        sample[, "LCMS-neg-raw"] <- NA
        sample[, "LCMS-neg-mzml"] <- NA
        sample[, "LCMS-neg-name"] <- NA
        sample[, "LCMS-pos-raw"] <- NA
        sample[, "LCMS-pos-mzml"] <- NA
        sample[, "LCMS-pos-name"] <- NA
        
        
        if (dir.exists("./raw/质谱数据")) {
          massfile <- c(
            list.files(path = "./raw/质谱数据/GCMS", pattern = "\\.D$", recursive = F, full.names = T),
            list.files(path = "./raw/质谱数据", pattern = "\\.wiff$", recursive = T, full.names = T),
            list.files(path = "./raw/质谱数据", pattern = "\\.raw$", recursive = T, full.names = T),
            list.files(path = "./raw/质谱数据", pattern = "\\.mzML$", recursive = T, full.names = T)
          )
          
          if (length(massfile) == 0) {
            warning("是否为外来数据，未发现质谱数据", immediate. = T)
          } else {
            sample <- GetRawDataInfo(sampleinfo = sample,
                                     listname = "GCMS",
                                     datawd = "./raw/质谱数据/GCMS",
                                     mzmlwd = "./raw/质谱数据/mzml/GCMS")
            sample <- GetRawDataInfo(sampleinfo = sample,
                                     listname = "LCMS-neg",
                                     datawd = "./raw/质谱数据/LCMS/neg",
                                     mzmlwd = "./raw/质谱数据/mzml/LCMS/neg")
            sample <- GetRawDataInfo(sampleinfo = sample,
                                     listname = "LCMS-pos",
                                     datawd = "./raw/质谱数据/LCMS/pos",
                                     mzmlwd = "./raw/质谱数据/mzml/LCMS/pos")
            
            if (any(c(
              any(!apply(sample[, grepl(pattern = "GCMS", x = colnames(sample))], 2, function(data) {
                any(is.na(data))
              })),
              any(!apply(sample[, grepl(pattern = "LCMS-neg", x = colnames(sample))], 2, function(data) {
                any(is.na(data))
              })),
              any(!apply(sample[, grepl(pattern = "LCMS-neg", x = colnames(sample))], 2, function(data) {
                any(is.na(data))
              }))
            ))) {
              
              
            } else {
              warning("是否为外来数据，未发现与样本关联质谱数据", immediate. = T)
            }
          }
        } else {
          warning("是否为外来数据，未发现质谱数据", immediate. = T)
        }
        
        print(type)
        if (grepl(pattern = "-双平台",x = type)) {
          print("双平台项目,分为lc与gc分别生成登记单")
          lcargs <- args
          lcargs$type <- "全谱代谢-LCMS"
          lcargs$xlsxname <- "内部分析单-lc.xlsx"
          do.call(what = ArrangeInfo,args = lcargs)
          gcargs <- args
          gcargs$type <- "全谱代谢-GCMS"
          gcargs$xlsxname <- "内部分析单-gc.xlsx"
          do.call(what = ArrangeInfo,args = gcargs)
          
          #20221222 当GCLC QC数不一致,分析单里删除多的QC
          gcinfo <- readxlsx(filename="内部分析单-gc.xlsx",sheet="样本基本信息")
          lcinfo <- readxlsx(filename="内部分析单-lc.xlsx",sheet="样本基本信息")
          
          gcinfo <- gcinfo[!(is.na(gcinfo$`GCMS-name`) & gcinfo$分组 == "QC"),]
          gcinfoqc <- gcinfo[gcinfo$分组 == "QC",]
          gcinfo <- gcinfo[gcinfo$分组 != "QC",]
          lcinfo <- lcinfo[!(is.na(lcinfo$`LCMS-pos-name`) & lcinfo$分组 == "QC"),]
          lcinfoqc <- lcinfo[lcinfo$分组 == "QC",]
          lcinfo <- lcinfo[lcinfo$分组 != "QC",]
          if(nrow(lcinfoqc)!=0){
            if(dim(lcinfoqc)[1] >= dim(gcinfoqc)[1]){
              infoqc <- lcinfoqc
              infoqc[1:dim(gcinfoqc)[1],"GCMS-raw"] <- gcinfoqc$`GCMS-raw`
              infoqc[1:dim(gcinfoqc)[1],"GCMS-mzml"] <- gcinfoqc$`GCMS-mzml`
              infoqc[1:dim(gcinfoqc)[1],"GCMS-name"] <- gcinfoqc$`GCMS-name`
            }else{
              infoqc <- gcinfoqc
              infoqc[1:dim(lcinfoqc)[1],"LCMS-neg-raw"] <- lcinfoqc$`LCMS-neg-raw`
              infoqc[1:dim(lcinfoqc)[1],"LCMS-neg-mzml"] <- lcinfoqc$`LCMS-neg-mzml`
              infoqc[1:dim(lcinfoqc)[1],"LCMS-neg-name"] <- lcinfoqc$`LCMS-neg-name`
              infoqc[1:dim(lcinfoqc)[1],"LCMS-pos-raw"] <- lcinfoqc$`LCMS-pos-raw`
              infoqc[1:dim(lcinfoqc)[1],"LCMS-pos-mzml"] <- lcinfoqc$`LCMS-pos-mzml`
              infoqc[1:dim(lcinfoqc)[1],"LCMS-pos-name"] <- lcinfoqc$`LCMS-pos-name`
            }
            infoqc$样本分析名称 <- paste0("QC",1:dim(infoqc)[1])
          }else{
            infoqc<-lcinfoqc
          }
          
          sample <- sample[sample$分组 != "QC",]
          sample <- rbind(sample,infoqc)
          sample <- sample[!(is.na(sample$`GCMS-name`) & sample$分组 == "QC"),]
          sample <- sample[!(is.na(sample$`LCMS-pos-name`) & sample$分组 == "QC"),]
          samplegroup[samplegroup$group == "QC","samples"] <- paste(sample[sample[,"分组"]=="QC","样本分析名称"],collapse = ",")
          
          gcinfo <- rbind(gcinfo,infoqc)
          gcinfo <- gcinfo[!(is.na(gcinfo$`GCMS-name`) & gcinfo$分组 == "QC"),]
          lcinfo <- rbind(lcinfo,infoqc)
          lcinfo <- lcinfo[!(is.na(lcinfo$`LCMS-pos-name`) & lcinfo$分组 == "QC"),]
          
          savexlsx1(data = gcinfo,filename = "内部分析单-gc.xlsx",sheet = "样本基本信息")
          savexlsx1(data = samplegroup,filename = "内部分析单-gc.xlsx",sheet = "样本分组信息")
          savexlsx1(data = lcinfo,filename = "内部分析单-lc.xlsx",sheet = "样本基本信息")
          savexlsx1(data = samplegroup,filename = "内部分析单-lc.xlsx",sheet = "样本分组信息")
          
        }
        
        sample[,"newname1"] <- stringi::stri_match_last_regex(str = sample$样本分析名称,pattern = "(.*[^0-9])([0-9]+)")[,2]
        sample[is.na(sample[,"newname1"]),"newname1"] <- sample$样本分析名称[is.na(sample[,"newname1"])]
        sample[,"newname2"] <- as.numeric(stringi::stri_match_last_regex(str = sample$样本分析名称,pattern = "(.*[^0-9])([0-9]+)")[,3])
        sample[is.na(sample[,"newname2"]),"newname2"] <- sample$样本分析名称[is.na(sample[,"newname2"])]
        sample <- sample[order(sample[,c("newname2")]),]
        sample <- sample[order(sample[,c("newname1")]),]
        sample <- sample[order(sample[,c("分组")]),]
        sample <- sample[,colnames(sample) != "newname2"]
        sample <- sample[,colnames(sample) != "newname1"]
        
        addsheet1(data = sample,wb = wb,sheet = "样本基本信息")
        addsheet1(data = samplegroup,wb = wb,sheet = "样本分组信息")
      } else {
        stop("分析确认单.xlsx中未找到样本信息或样本分组信息")
      }
      
      if (info[c("原始文件"), 2] != "有") {
        warning("未找到搜库结果或质谱数据，请确认是否存在文件或命名不正确")
      }
      
      addsheet1(data = info,wb = wb,sheet = "分析基本信息")
    } else {
      stop("分析确认单.xlsx中未找到分析基本信息")
    }
    
    try({
      savewb(wb = wb,filename = xlsxname,overwrite = overwrite)
    })
  } else {
    stop("未找到分析确认单.xlsx")
  }
  
  print("complete!")
}

#' @export
ArrangeInfoinpath <- function(savepath = "./",
                              ...){
  runinpath(path = savepath,
            moudle = ArrangeInfo,
            ...)
}

#' @export
getkeggspecies <- function(species = "人"){
  
  data <- readdata(databasepath(path = "database/物种映射/map.xlsx"))
  species <- gsub(pattern = " ",replacement = "",x = species)
  if(species %in% data[,2]){
    return(data[data[,2] == species,1][1])
  }else{
    return(NA) 
  }
}

#' @export
addkeggspecies <- function(species = "人",keggspecies = "hsa"){
  
  data <- readdata(databasepath(path = "database/物种映射/map.xlsx"))
  if(!is.na(keggspecies)){
    data[dim(data)[1]+1,] <- c(keggspecies,species) 
    data <- data[!duplicated(data[,2]),]
    savexlsx(data = data,filename = databasepath(path = "database/物种映射/map.xlsx",exist = F),sheet = "映射关系")
  }
  
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-fn","--fromname",default = "分析确认单.xlsx",help = "分析单文件名")
  parser$add_argument("-xn","--xlsxname",default = "内部分析单.xlsx",help = "保存文件名")
  parser$add_argument("-tp","--type",default = "",help = "项目类型")
  parser$add_argument("-ks","--keggspecies",default = "ko", help = "物种")
  parser$add_argument("-ow","--overwrite",default = F,help = "是否覆盖原始分析单",action='store_true')
  
  args <- parser$parse_args()
  
  if(args$type == ""){args$type <- NULL}
  
  result <- do.call(what = ArrangeInfo,args = args) 
  
}
