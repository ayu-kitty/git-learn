#!/opt/conda/bin/Rscript

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser1 <- ArgumentParser(Parsername = "parser1",add_help = F)
  group1 <- parser1$add_argument_group('基本参数')
  group1$add_argument("-m","--massrange",default = NULL,type= "double",help = "mz范围,默认NULL",
                      nargs = 2,dest = "mass.range")
  group1$add_argument("-re","--resolution",default = 5, type= "double",help = "分辨率,默认5")
  group1$add_argument("-u","--units",default = "ppm",help = "分辨率单位,默认ppm")
  group1$add_argument("-mr","--moderange",default = c("neg","pos"),help = "正负离子模式",nargs = "+",
                      choices = c("neg","pos"))
  group1$add_argument("-a","--asp",default = 1, type= "double",help = "长宽分辨率比")
  group1$add_argument("-t","--type",default = "空间代谢组",help = "空代产品线类型:空间代谢组 or 中药空间代谢组 or 空间代谢组-waters")    
  group1$add_argument("-tr","--tolerance",default = 10, type= "double",help = "峰对齐参数,默认5")
  group1$add_argument("-w","--waters",default = F,help = "是否使用waters仪器",action ='store_true')

  parser2 <- parser1$newparser(Parsername = "parser2",add_help = F)
  group2 <- parser2$add_argument_group('定性参数')
  group2$add_argument("-p","--ppm",default = 5, type= "double",help = "定性ppm参数")
  group2$add_argument("-pn","--polymernum",default = 3, type= "integer",help = "胶聚合物最小数量")
  group2$add_argument("-pr","--polymerratio",default = 0.5, type= "double",help = "聚合物最大表达范围")
  group2$add_argument("-pm","--polymerminmz",default = 300, type= "double",help = "聚合物最小mz")
  group2$add_argument("-st","--species_tissue",default = "小鼠",help = "物种_组织",required = T)
  group2$add_argument("-od","--outdatabse",default = NULL,help = "外部数据库，如是中药数据库输入tcm")
  
  parser3 <- parser1$newparser(Parsername = "parser3",add_help = F)
  group3 <- parser3$add_argument_group('标准化参数')
  group3$add_argument("-nm","--normethod",default = "none", help = "是否不进行内标归一化处理,不进行填写none",
                      choices = c("none","reference","tic"))
  group3$add_argument("-nrn","--norreferencemzneg",default = 89.02319, type= "double", help = "reference归一化模式下,选择的负离子",nargs = "+")
  group3$add_argument("-nrp","--norreferencemzpos",default = 317.1140, type= "double", help = "reference归一化模式下,选择的正离子",nargs = "+")
  group3$add_argument("-sd","--stripdeal",default = F, help = "是否不进行条纹处理",action='store_true',
                      dest = "stripdeal")
  group1$add_argument("-nb","--number",default = 100000, type= "integer",help = "分段矫正参数")

  
  parser4 <- parser1$newparser(Parsername = "parser4",add_help = F)
  group4 <- parser4$add_argument_group('批次校正参数')
  group4$add_argument("-sb","--sampletobatch",default = NULL,help = "根据样本提供批次,如1:A1,B1,C1 2:A2,B2,C2",
                      nargs = "+")
  group4$add_argument("-cm","--combatmod",default = NULL,help = "根据样本提供感兴趣分组,如A:A1,A2 B:B1,B2 C:C1,C2",
                      nargs = "+")
  
  parser5 <- parser1$newparser(Parsername = "parser5",add_help = F)
  group5 <- parser5$add_argument_group('对齐参数')
  group5$add_argument("-nrf","--negrefmz",default = c(89.0244,124.0074,133.0142,171.1391,255.233,327.233,554.262,714.5079,834.5283), type= "double",help = "负离子的质量轴校正",nargs = "+")
  group5$add_argument("-prf","--posrefmz",default = c(104.1075,147.1128,149.0233,156.042,239.1642,301.141,556.2766,578.2585,606.2942,734.5694,772.5238,798.5401,826.5705,844.5334), type= "double",help = "正离子的质量轴校正",nargs = "+")
  
  parser <- parser1$newparser(parents = list("[parser1,parser2,parser3,parser4,parser5]"),Parsername = "parser_all")
  parser2_1 <- parser1$newparser(parents = list("[parser1,parser2]"),Parsername = "parser2_1")
  parser3_1 <- parser1$newparser(parents = list("[parser1,parser3]"),Parsername = "parser3_1")
  parser4_1 <- parser1$newparser(parents = list("[parser1,parser4]"),Parsername = "parser4_1")
  parser5_1 <- parser1$newparser(parents = list("[parser1,parser5]"),Parsername = "parser5_1")
  args <- parser$parse_args()
  args1 <- parser1$parse_known_args()[[1]]
  args2_1 <- parser2_1$parse_known_args()[[1]]
  args3_1 <- parser3_1$parse_known_args()[[1]]
  args4_1 <- parser4_1$parse_known_args()[[1]]
  args5_1 <- parser4_1$parse_known_args()[[1]]
  
  writeinfo()
  
  if(!is.null(args4_1$sampletobatch)){
    txt <- args$sampletobatch
    txt <- strsplit(txt,split = ":")
    txt <- Reduce(rbind,txt)
    txt <- as.data.frame(txt)
    txt <- apply(txt, 1, splitdata)
    txt <- Reduce(rbind, txt)
    colnames(txt) <- c("batch","samplename")
    args4_1$sampletobatch <- txt
    print("~查看批次组别")
    print(txt)
    if(any(duplicated(txt$samplename))){
      stop("样本有重复分组")
    }
  }
  
  if(!is.null(args4_1$combatmod)){
    txt <- args4_1$combatmod
    txt <- strsplit(txt,split = ":")
    txt <- Reduce(rbind,txt)
    txt <- as.data.frame(txt)
    txt <- apply(txt, 1, splitdata)
    txt <- Reduce(rbind, txt)
    colnames(txt) <- c("modlist","samplename")
    args4_1$combatmod <- txt
    print("~查看感兴趣分组")
    print(txt)
    if(any(duplicated(txt$samplename))){
      stop("样本有重复分组")
    }
  }
  
  # 峰对齐
  args5_1$number <- args3_1$number
  createdir(filename = "./sample/peak",linkdir = T)
  result <- do.call(what = imzmlReferenceMZ,args = args5_1)
  gc(reset = TRUE)
  
  args1$number <- NULL
  args1$negrefmz <- NULL
  args1$posrefmz <- NULL

  # 归一化
  cleanenv()
  createdir(filename = "./sample/imzml-pre",linkdir = T)
  result <- do.call(what = imzmlNormalization,args = args3_1)
  gc(reset = TRUE)
  
  args1$waters <- NULL
  
  # 去除背景
  cleanenv()
  result <- do.call(what = imzmlrmbgmz,args = args1)
  gc(reset = TRUE)
  
  # 定性
  cleanenv()
  result <- do.call(what = imzmlpreQualitative,args = args2_1)
  gc(reset = TRUE)
  
  args2_1$waters <- NULL
  args3_1$waters <- NULL
  args4_1$waters <- NULL
  
  # 预处理
  cleanenv()
  createdir(filename = "./sample/pre",linkdir = T)
  createdir(filename = "./sample/final",linkdir = T)
  args1$tolerance <- NULL
  result <- do.call(what = MulPreDeal,args = args1)
  gc(reset = TRUE)
  
  # 批次校正
  cleanenv()
  createdir(filename = "./sample/adjustdata",linkdir = T)
  result <- do.call(what = Combatspace,args = args4_1)
  gc(reset = TRUE)
  
  # 聚类分析
  cleanenv()
  createdir(filename = "./sample/cluster",linkdir = T)
  result <- MulanaCluster()
  gc(reset = TRUE)
  
  # 数据提取
  if(length(list.files(path = "./sample/final/",pattern = "^qc_data")) > 0){
    samplefrom <- "./sample/adjustdata/"
  }else{
    samplefrom <- "./sample/final/"
  }
  
  createdir(filename = "./sample/area",linkdir = T)
  GetAllData(asp = args$asp,samplefrom = samplefrom)
  GetClusterData(asp = args$asp,samplefrom = samplefrom)
  GetAreaData(asp = args$asp,samplefrom = samplefrom)
  GetMulData(asp = args$asp,samplefrom = samplefrom)
  gc(reset = TRUE)
  
  # 比较分析
  createdir(filename = "./sample/analysis",linkdir = T)
  SpaceAnalysisData()
  gc(reset = TRUE)
  
  # 生成项目报告
  movetoreport(asp = args$asp)
  makespacereport(type = args$type)
  
  writeinfo(endtime = T)
}
