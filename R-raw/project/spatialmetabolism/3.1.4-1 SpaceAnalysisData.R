#!/opt/conda/bin/Rscript

#' 空代数据比较分析
#'
#' @param maxnode 分析样本点上限
#' @param onetoone 逻辑，是否结果一对一
#' @param qualitativefrom 定性结果文件路径
#' @param datafrom 数据路径
#' @param savepath 存储路径
#' @param moderange 正负离子模式
#'
#' @export
SpaceAnalysisData <- function(maxnode = 1000,
                              onetoone = F,
                              infopath = "项目登记单.xlsx",
                              qualitativefrom = "./sample/qualitative/Qualitative.xlsx",
                              datafrom = "./sample/area/data/",
                              savepath = "./sample/analysis/",
                              moderange = c("neg", "pos")) {
  print("~~~~~新分析流程-20230320~~~~~")
  
  wd <- getwd()
  
  set.seed(123)
  
  if (!file.exists(infopath)) {
    stop("不存在项目登记单.xlsx")
  }
  
  # 项目信息获取
  info <- readdata(filename = infopath, sheet = "项目登记单")
  species <- info[info[, 1] == "映射简写", 2]
  # 分组信息获取
  group <- readdata(filename = infopath, sheet = "分组信息")
  # 比较组信息获取
  compare <- readdata(filename = infopath, sheet = "比较组信息")
  compare <- compare[!is.na(compare$比较组),]
  
  if (dim(compare)[1] == 0) {
    return(F)
  }
  
  for (mode in moderange) {
    compare1 <- compare[compare$`模式` == mode | compare$`模式` == "both", ]
    if (dim(compare1)[1] == 0) {
      print(paste0(mode, "模式下未发现比较组"))
      next
    }
    
    for (testmode in c("pixel_level", "sample_level")) {
      
      if(testmode == "pixel_level"){
        compare2 <- compare1[compare1$`比较模式` == testmode | compare1$`比较模式` == "all" |compare1$`比较模式` == "both", ]
      }else{
        compare2 <- compare1[compare1$`比较模式` == testmode | compare1$`比较模式` == "mean" |compare1$`比较模式` == "both", ]
      }
      
      compare2 <- compare2[!duplicated(compare2$`比较组`), ]
      
      if (dim(compare2)[1] == 0) {
        print(paste0(mode, "模式下未发现", testmode, "类型比较分析"))
        next
      }
      
      for (i in 1:dim(compare2)[1]) {
        groupname <- compare2[i, "比较组"]
        print(paste0("###################进行", mode, "模式下", testmode, "类型的", groupname, "组分析##################"))
        groupname <- unlist(strsplit(groupname, split = "/"))
        
        classname <- NULL
        mzdata <- NULL
        
        for (k in seq_len(length(groupname))) {
          testname <- groupname[k]
          testgroup <- group[group$`分组` == testname, ]
          testgroup <- testgroup[testgroup$`模式` == mode | testgroup$`模式` == "both", ]
          if (dim(testgroup)[1] == 0) {
            stop(paste0(testname, "分组未在分组信息中找到相应信息，请完善分组信息"))
          }
          testgroup <- paste(testgroup$样品, testgroup$选区, sep = "-")
          
          print(paste0("数据读取:",datafrom, testgroup[1], "-", mode, "-", testmode, ".txt"))
          testdata <- read.table(
            file = paste0(datafrom, testgroup[1], "-", mode, "-", testmode, ".txt"),
            sep = "\t", header = T, check.names = F
          )
          if (dim(testdata)[2] > maxnode + 1) {
            testdata <- testdata[, c(1, sample(x = 2:dim(testdata)[2], size = maxnode)), drop = F]
          }
          if (length(testgroup) > 1) {
            for (j in 2:length(testgroup)) {
              print(paste0("数据读取:",datafrom, testgroup[j], "-", mode, "-", testmode, ".txt"))
              testdata1 <- read.table(
                file = paste0(datafrom, testgroup[j], "-", mode, "-", testmode, ".txt"),
                sep = "\t", header = T, check.names = F
              )
              testdata1 <- testdata1[, c(T, !(colnames(testdata1)[-1] %in% colnames(testdata))), drop = F]
              if (dim(testdata1)[2] > maxnode + 1) {
                testdata1 <- testdata1[, c(1, sample(x = 2:dim(testdata1)[2], size = maxnode)), drop = F]
              }
              testdata <- cbind(testdata, testdata1[, -1, drop = F])
            }
          }
          
          if (!is.null(mzdata)) {
            testdata <- testdata[, c(T, !(colnames(testdata)[-1] %in% colnames(mzdata))), drop = F]
          }
          
          if (dim(testdata)[2] > maxnode + 1) {
            testdata <- testdata[, c(1, sample(x = 2:dim(testdata)[2], size = maxnode)), drop = F]
          }
          # 用0代替由高斯处理后的负值
          testdata[testdata < 0] <- 0
          
          testset <- data.frame(
            samplename = names(testdata)[-1],
            class = testname, stringsAsFactors = F
          )
          
          
          classname <- rbind(classname, testset)
          if (is.null(mzdata)) {
            mzdata <- testdata
          } else {
            mzdata <- cbind(mzdata, testdata[, -1, drop = F])
          }
        }
        
        metadata <- readdata(filename = qualitativefrom, sheet = mode)
        metadata[, "Ion mode"] <- mode
        mzdata <- cbind(metadata, mzdata[, -1, drop = F])
        mzdata[, "mz"] <- format(mzdata[, "mz"], digits = 5, nsmall = 5, trim = T)
        mzdata <- mzdata[!is.na(mzdata$Metabolites), ]
        
        colnames(classname) <- c("sample", "group")
        
        data <- meta::ArrangeInfoFromData(data = mzdata,
                                          index = data.frame(
                                            "索引" = c("ID", "Metabolites", "Compound ID", "Ion mode", "KEGG"),
                                            "列名" = c("mz", "Metabolites", "HMDB", "Ion mode", "KEGG"),
                                            stringsAsFactors = F
                                          ),
                                          group = classname,
                                          compare = data.frame("比较组" = compare2[i, "比较组"], "配对" = F, "分组" = "group", stringsAsFactors = F),
                                          type = "空间代谢组",
                                          rsd = NA,
                                          missvalue = NA,
                                          species = species,
                                          saveregistration = F)
        
        if(onetoone | !(paste0(mode, "-all") %in% getsheetname(qualitativefrom))){
          onetoone2 <- T
          class(data) <- "SpaceM2"
        }else{
          class(data) <- "SpaceM"
          
          information2 <- readdata(filename = qualitativefrom, sheet = paste0(mode, "-all"))
          information2[, "Ion mode"] <- mode
          information2[, "mz"] <- format(information2[, "mz"], digits = 5, nsmall = 5, trim = T)
          data$data$predata$information2 <- information2
        }
        
        # 比较组运算
        # data <- meta::datastatistics(data) # 各比较组多元统计参数填写
        # data <- meta::automulstatistics(data) # 各比较组多元统计运算
        data <- lmbio::run_mulstatisticsana(data,all = F)
        unlink(x = "oecloud",recursive = T)
        unlink(x = ".snakemake",recursive = T)
        data <- meta::getvip(data) # 获取vip数据，VIP.xlsx或脚本内数据
        data <- meta::countdif(data)
        
        filepath <- paste0(savepath, paste(groupname, collapse = "_"), "-", testmode, "/", mode)
        setwddir(filename = filepath)
        
        data <- meta::createlist_old(data,move = F)#7.1由createlist_old改为createlist
        save(data, file = "raw.RData")
        
        setwd(wd)
      }
    }
  }
  
  # 正负离子联合结果生成
  
  difpath <- dir(path = savepath, full.names = T, recursive = F)
  
  for (p in difpath) {
    path_dir <- dir(path = p, full.names = T, recursive = F)
    negdir <- paste0(p, "/neg")
    posdir <- paste0(p, "/pos")
    
    if (dir.exists(negdir) & dir.exists(posdir)) {
      print(paste0("开始合并 \'", p, "\' 下neg和pos的 \'差异代谢物.xlsx\' 的代谢物列表"))
      path_dir_ <- dir(path = path_dir,
                       pattern = "差异代谢物$",
                       full.names = T,
                       recursive = T,
                       include.dirs = T)
      difpath_ <- dir(path = path_dir_,
                      pattern = "差异代谢物.xlsx",
                      full.names = T)
      
      # neg
      negdir <- grep(pattern = "neg", x = difpath_, value = T)
      
      if (length(negdir) == 0) {
        next
      }
      
      # pos
      posdir <- grep(pattern = "pos", x = difpath_, value = T)
      
      if (length(posdir) == 0) {
        next
      }
      
      sheetname <- getsheetname(filename = negdir)
      
      for(i in 1:length(sheetname)){
        difdata1 <- readdata(filename = negdir, sheet = sheetname[i])
        n <- max(which(grepl(pattern = "average", x = names(difdata1))))
        difdata1 <- difdata1[, 1:n]
        
        
        difdata2 <- readdata(filename = posdir, sheet = sheetname[i])
        n <- max(which(grepl(pattern = "average", x = names(difdata2))))
        difdata2 <- difdata2[, 1:n]
        
        difdata <- rbind(difdata1, difdata2)
        
        setwddir(paste0(p, "/all"))
        
        savexlsx1(data = difdata,
                  filename = "差异代谢物.xlsx",
                  sheet = sheetname[i])
        
        if(i == length(sheetname)){
          print("~~~~~~合并完成，并保存~~~~~~")
          meta::keggrich(name = "差异代谢物.xlsx",
                         needgroup = i,
                         meta = "Metabolites",
                         logfc = "log2(FC)",
                         kegg = "KEGG",
                         species = species)
        }
        setwd(wd)
      }
    }
  }
  return(T)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-q","--qualitativefrom",default = "./sample/qualitative/Qualitative.xlsx",
                      help = "定性结果路径,默认./sample/qualitative/Qualitative.xlsx")
  parser$add_argument("-i","--infopath",default = "项目登记单.xlsx",
                      help = "项目登记单路径,默认项目登记单.xlsx")
  parser$add_argument("-d","--datafrom",default = "./sample/area/data/",
                      help = "数据路径,默认./sample/area/data/")
  parser$add_argument("-s","--savepath",default = "./sample/analysis/",
                      help = "存储路径,默认./sample/analysis/")
  parser$add_argument("-m","--maxnode",default = 1000, type= "double",help = "长宽分辨率比")
  parser$add_argument("-mr","--moderange",default = c("neg","pos"),help = "正负离子模式",nargs = "+",
                      choices = c("neg","pos"))
  
  args <- parser$parse_args()
  
  writeinfo()
  
  createdir(filename = args$savepath,linkdir = T)
  
  result <- do.call(what = SpaceAnalysisData,args = args)
  
  writeinfo(endtime = T)
  
}
