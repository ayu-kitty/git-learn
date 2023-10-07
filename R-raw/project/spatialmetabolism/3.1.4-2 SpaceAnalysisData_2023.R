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
SpaceAnalysisData_2023 <- function(infopath = "项目登记单.xlsx",
                                   qualitativefrom = "sample/qualitative/Qualitative.xlsx",
                                   datafrom = "sample/area/data",
                                   savepath = "sample/analysis",
                                   moderange = c("neg", "pos"),
                                   maxnode = 100,
                                   addana = T) {
  
  wd <- getwd()
  
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
  
  sheetname <- getsheetname(qualitativefrom)
  if("neg" %in% sheetname){
    qualitativeneg <- readdata(filename = qualitativefrom, sheet = "neg")
  }else{qualitativeneg <- NULL}
  if("pos" %in% sheetname){
    qualitativepos <- readdata(filename = qualitativefrom, sheet = "pos")
  }else{qualitativepos <- NULL}
  qualitativedata <- rbind(qualitativeneg,qualitativepos)
  
  movebackground(db = "background",
                 dbfrom = "/data/hstore4/database/kegg/kegg/",
                 keggmapfrom = "/data/hstore4/database/kegg/kegg/",
                 org = species,
                 inputfile = qualitativedata)
  
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
        print(paste0("~",mode, "模式下未发现", testmode, "类型比较分析"))
        next
      }
      
      compare3 <- compare2[,"比较组",drop = F]
      
      groupname <- group[group$模式 == mode | group$模式 == "both","分组"]
      groupname <- unique(groupname)
      
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
        
        print(paste0("数据读取:",datafrom,"/",testgroup[1], "-", mode, "-", testmode, ".txt"))
        testdata <- read.table(file = paste0(datafrom,"/",testgroup[1], "-", mode, "-", testmode, ".txt"),
                               sep = "\t", header = T, check.names = F)
        
        # if (dim(testdata)[2] > maxnode + 1) {
        #   testdata <- testdata[, c(1, sample(x = 2:dim(testdata)[2], size = maxnode)), drop = F]
        # }
        
        if (length(testgroup) > 1) {
          for (j in 2:length(testgroup)) {
            print(paste0("数据读取:",datafrom,"/",testgroup[j], "-", mode, "-", testmode, ".txt"))
            testdata1 <- read.table(file = paste0(datafrom,"/",testgroup[j], "-", mode, "-", testmode, ".txt"),
                                    sep = "\t", header = T, check.names = F)
            
            testdata1 <- testdata1[, c(T, !(colnames(testdata1)[-1] %in% colnames(testdata))), drop = F]
            # if (dim(testdata1)[2] > maxnode + 1) {
            #   testdata1 <- testdata1[, c(1, sample(x = 2:dim(testdata1)[2], size = maxnode)), drop = F]
            # }
            testdata <- cbind(testdata, testdata1[, -1, drop = F])
          }
        }
        
        # testdata  <<- testdata 
        if(testmode == "pixel_level"){
          if(addana){
            
            if(dim(testdata)[2]-1 <= maxnode){
              maxnode2 <- dim(testdata)[2]-1
              num <- 5
            }else{
              maxnode2 <- maxnode
              num <- ceiling(((dim(testdata)[2]-1)/maxnode))*5
              if(num > 100){
                num <- 100
              }
            }
            
            if(maxnode2 > 10){
              testdata2 <- testdata[,1,drop = F]
              # testdata3 <<- testdata
              for(i in 1:maxnode2){
                testdata2[,paste0(testname,"-",mode,"-",i)] <- apply(testdata[,sample(2:(dim(testdata)[2]-1),size = num),drop = F],MARGIN = 1,FUN = mean)
              }
              testdata <- testdata2
            }
            
          }else{
            if (!is.null(mzdata)) {testdata <- testdata[, c(T, !(colnames(testdata)[-1] %in% colnames(mzdata))), drop = F]}
            if (dim(testdata)[2] > maxnode + 1) {testdata <- testdata[, c(1, sample(x = 2:dim(testdata)[2], size = maxnode)), drop = F]}
          }
        }

        # 用0代替由高斯处理后的负值
        testdata[testdata < 0] <- 0
        
        testset <- data.frame(samplename = names(testdata)[-1],
                              class = testname,
                              stringsAsFactors = F)
        
        classname <- rbind(classname, testset)
        if (is.null(mzdata)) {
          mzdata <- testdata
        } else {
          mzdata <- cbind(mzdata, testdata[, -1, drop = F])
        }
      }
      
      mzdata <- mzdata[,!duplicated(colnames(mzdata)),drop = F]
      mzdata[, "mz"] <- format(mzdata[, "mz"], digits = 5, nsmall = 5, trim = T)
      metadata <- readdata(filename = qualitativefrom, sheet = mode)
      metadata[, "Ion mode"] <- mode
      mzdata <- merge(metadata,mzdata,by = "mz",sort = F)
      mzdata <- mzdata[!is.na(mzdata$Metabolites), ]
      mzdata <- mzdata[,c(2:1,3:dim(mzdata)[2])]
      
      class <- classname
      samplename <- unique(class[,1])
      class2 <- data.frame("sample" = samplename,row.names = samplename)
      class2 <- class2[,0,drop = F]
      groupname <- unique(class[,2])
      for ( i in 1:length(groupname)) {
        class2[class[class[,2] == groupname[i],1],groupname[i]] <- groupname[i]
      }
      classname <- class2
      
      organizedata(rawdatafile = mzdata,
                   rawclassfile = classname,
                   rawcomparefile = compare3,
                   datafile = paste0(savepath,"/raw/",mode,"/",testmode,"/rawdata/datafile.txt"),
                   infofile = paste0(savepath,"/raw/",mode,"/",testmode,"/rawdata/infofile.txt"),
                   classfile = paste0(savepath,"/raw/",mode,"/",testmode,"/rawdata/classfile.yaml"),
                   classtypefile = paste0(savepath,"/raw/",mode,"/",testmode,"/rawdata/classtype.xlsx"),
                   comparefile = paste0(savepath,"/raw/",mode,"/",testmode,"/rawdata/compare.yaml"))
      
      system(paste0("source /etc/profile;source ~/.bashrc;",
                    packagepath(path = "python/lmbio/advancedanalysis/DiffAnalysis.py"),
                    " --omic M",
                    " --org ",species,
                    " --tempfilepath ",paste0(savepath,"/raw/",mode,"/",testmode),
                    " --projectpath ",paste0(savepath,"/result/",mode,"/",testmode),
                    " --oecloud -ff 0 -vf 1 -pf 0.05 --report F"))
      
    }
  }
  
  # 正负离子联合结果生成
  for (testmode in c("pixel_level", "sample_level")) {
    
    negfile <- list.files(paste0(savepath,"/raw/","neg","/",testmode,"/diff_filter"))
    posfile <- list.files(paste0(savepath,"/raw/","pos","/",testmode,"/diff_filter"))
    file <- intersect(negfile,posfile)
    base::print(file)
    for ( i in 1:length(file)) {
      negdata <- get_diffana_data_2(data = paste0(savepath,"/raw/","neg","/",testmode,"/diff_filter/",file[i]),savedata = F,datafrom = "filterdata",adddata = F)
      posdata <- get_diffana_data_2(data = paste0(savepath,"/raw/","pos","/",testmode,"/diff_filter/",file[i]),savedata = F,datafrom = "filterdata",adddata = F)
      negdata <- negdata[,colnames(negdata) %in% colnames(posdata),drop = F]
      posdata <- posdata[,colnames(posdata) %in% colnames(negdata),drop = F]
      bothdata <- rbind(negdata,posdata)
      group <- gsub(pattern = ".rds",replacement = "",x = file[i])
      group <- gsub(pattern = "diff_filter-",replacement = "",x = group)
      savexls(data = bothdata,
              filename = paste0(savepath,"/result/","both","/",testmode,"/通路富集/",group,".xls"))
      
      # base::print(paste0("source /etc/profile;source ~/.bashrc;",
      #                    packagepath(path = "python/lmbio/protein/functionenrich.py"),
      #                    " --omic M",
      #                    " --org ",species,
      #                    ifelse("ChEBI" %in% colnames(bothdata)," --enrichtype all"," --enrichtype K"),
      #                    " --diffile ",paste0(savepath,"/result/","both","/",testmode,"/通路富集/",group,".xls"),
      #                    " --savepath ",paste0(savepath,"/result/","both","/",testmode,"/通路富集")))
      
      system(paste0("source /etc/profile;source ~/.bashrc;",
                    packagepath(path = "python/lmbio/protein/functionenrich.py"),
                    " --omic M",
                    " --org ",species,
                    ifelse("ChEBI" %in% colnames(bothdata)," --enrichtype all"," --enrichtype K"),
                    " --diffile ",paste0(savepath,"/result/","both","/",testmode,"/通路富集/",group,".xls"),
                    " --savepath ",paste0(savepath,"/result/","both","/",testmode,"/通路富集")))
      
    }
    
  }
  
  return(T)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-q","--qualitativefrom",default = "sample/qualitative/Qualitative.xlsx",
                      help = "定性结果路径,默认sample/qualitative/Qualitative.xlsx")
  parser$add_argument("-i","--infopath",default = "项目登记单.xlsx",
                      help = "项目登记单路径,默认项目登记单.xlsx")
  parser$add_argument("-d","--datafrom",default = "sample/area/data",
                      help = "数据路径,默认sample/area/data")
  parser$add_argument("-s","--savepath",default = "sample/analysis",
                      help = "存储路径,默认sample/analysis")
  parser$add_argument("-m","--maxnode",default = 100, type= "double",help = "长宽分辨率比")
  parser$add_argument("-mr","--moderange",default = c("neg","pos"),help = "正负离子模式",nargs = "+",
                      choices = c("neg","pos"))
  
  args <- parser$parse_args()
  
  writeinfo()
  
  createdir(filename = args$savepath,linkdir = T)
  
  result <- do.call(what = SpaceAnalysisData_2023,args = args)
  
  writeinfo(endtime = T)
  
}
