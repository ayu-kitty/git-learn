################################################################################
################################################################################
# / Rscript Pipeline Name: 蛋白质组学WGCNA4_分析流程
# / Author: jingxinxing
# / Date: 2021-04(20210513)
# / Version: 4.2.0
# /  New:
#      1.绘图代码优化，根据数据自动调整图形的大小
#      2.为后面的大升级做准备
#      3.主函数RUNWGCNA参数优化
################################################################################
################################################################################

#' WGCNA整套分析
#'
#' @param path 运行路径
#' @param wgcna_project_type 运行类型，G基因、P蛋白、M代谢
#' @param MiceImputeMethod 空值填充方法
#' @param kegg_species_code 物种
#' @param wgcna_nthreads 运算核心数
#' @param MinModuleSize 最小模块数
#' @param type 运算方式
#' @param NAFillingMethod 缺失值填充方式
#' @param BioFunEnrichAnalysis 是否进行富集分析
#' @param projectname 项目名
#' @param exprmat_file 表达矩阵名
#' @param traits_file 性状矩阵名
#' @param corType 相关性计算方式
#' @param removesample 移除样本
#'
#' @export
RUNWGCNA <- function(path = "./",
                     wgcna_project_type = "P",
                     projectname = paste0("WGCNA_",wgcna_project_type),
                     exprmat_file = "ExprData_raw.xlsx",
                     traits_file = "Sample_Traits.xlsx",
                     MiceImputeMethod = 5,
                     kegg_species_code = "hsa",
                     wgcna_nthreads = 5,
                     MinModuleSize = 5, ## 用于模块检测的最小模块数
                     type = "unsigned",
                     NAFillingMethod = "min",
                     BioFunEnrichAnalysis = F,
                     corType = "pearson",
                     removesample = NULL) {

  ## ########################
  ## 测试参数
  # path = "./"
  # wgcna_project_type = "P"
  # projectname = paste0("WGCNA_",wgcna_project_type)
  # exprmat_file = "ExprData_raw.xlsx"
  # traits_file = "Sample_Traits.xlsx"
  # MiceImputeMethod = NULL
  # kegg_species_code = "hsa"
  # wgcna_nthreads = 5
  # type = "unsigned"
  # NAFillingMethod = "min"
  # BioFunEnrichAnalysis = F
  # removesample = NULL
  # MinModuleSize = 5
  # corType = "pearson"
  ############################

  ## 2.WGCNA分析开始日期时间
  pt <- proc.time()

  wgcna_project_type_cn  <- switch(wgcna_project_type,
                                   "G" = "基因组、转录组",
                                   "P" = "蛋白质组质",
                                   "M" = "代谢物质组学")
  wgcna_project_type_en  <- switch(wgcna_project_type,
                                   "G" = "Gene",
                                   "P" = "Protein",
                                   "M" = "Metabolite")
  wgcna_project_type_name  <- switch(wgcna_project_type,
                                     "G" = "GeneName",
                                     "P" = "Protein_ID",
                                     "M" = "Metabolites")

  print(paste0(wgcna_project_type_cn,"WGCNA分析及绘图日志开始时间：",Sys.time()))


  print("************************************************************************")
  print("分析开始：Start.........................................................")
  print("************************************************************************")

  ## 3.设置WGCNA分析工作主目录
  setwd(path)
  wd <- getwd()
  print(paste0("当前工作路径为：",wd))

  ## 4.批量导入R包
  ### 4.1 导入分析用的R包 ###
  analibpackages <- c("WGCNA","mice","openxlsx","meta","dplyr","stringr","ggplot2","pheatmap","igraph")
  sapply(analibpackages, library, character.only = T)

  ## 6.WGCNA分析R环境变量及线程设置
  ### 6.1 options()函数允许用户设置和检查影响R计算和显示其结果的方式的各种全局选项
  options(stringsAsFactors = FALSE)
  op <- par(no.readonly = TRUE) # 将R默认图形环境变量保存为op，为后面的环境变量恢复提供基础

  ### 6.2 打开多线程
  print("----------------------------------------------------------------------")
  enableWGCNAThreads(nThreads = wgcna_nthreads) # WGCNA分析进程的线程设置
  print(paste0("当前WGCNA分析进程的线程数为：", WGCNAnThreads(), "个"))
  print("----------------------------------------------------------------------")

  ## 7.创建WGCNA分析目录系统
  createdir(projectname,force = T)
  createdir(paste0(projectname,"/1.Data_Processing/1.1Data"))

  # 一、数据处理
  # WGCNA分析的输入数据是基因、蛋白或代谢物表达矩阵数据和样品性状矩阵数据，经过对数据进行Z-score标准化、批次效应校正、缺失值处理后得到WGCNA分析所需的Cleandate。WGCNA后续分析都是基于Cleandata进行的。

  ## 1.数据处理
  ### 1.1 读取原始的蛋白表达矩阵数据
  print(paste0("原始表达矩阵文件名为： ",exprmat_file))
  trusted_prot_data <- readdata(filename = exprmat_file,sheet = 1,rowNames = T)
  ### 1.4 读取sample_information.xlsx文件中的样品和比较组信息
  print(paste0("性状文件名为： ",traits_file))
  traitData <- readdata(filename = traits_file,sheet = 1,rowNames = T)

  # 复制一份蛋白原始表达矩阵数据到WGCNA_P/1.Data_Processing/1.1Data/目录内
  file.copy(from = exprmat_file, to = paste0("./",projectname,"/1.Data_Processing/1.1Data/"), copy.mode = TRUE)
  file.copy(from = traits_file, to = paste0("./",projectname,"/1.Data_Processing/1.1Data/"), copy.mode = TRUE)

  ### 1.5 获取Clean可信蛋白-样品数据矩阵
  sam_name <- row.names(traitData)
  sam_ma <- match(sam_name, colnames(trusted_prot_data)) # 获取样品在数据集中的位置索引
  trusted_prot_clean_data <- trusted_prot_data[,unique(sam_ma)] # 根据样品名称位置索引提取干净的蛋白样品表达矩阵数据

  ### 1.6 Clean可信蛋白-样品数据矩阵的行、列转置：t()函数
  trusted_prot_clean_data_t <- t(trusted_prot_clean_data)

  ### 1.7 表达矩阵数据的标准化（Normalization）处理Z-score法（其他标准化方法可以后续填加进来）：scale()函数；
  trusted_prot_clean_data_t <- scale(trusted_prot_clean_data_t, center = T, scale = T)
  trusted_prot_clean_data_t <- as.data.frame(trusted_prot_clean_data_t)

  ### 2.1 蛋白表达矩阵数据格式和内容检测：主要是蛋白和样品的缺失值检测和处理
  gsp <- goodSamplesGenes(trusted_prot_clean_data_t, verbose = 3) # good samples proteins, gsp
  # Flagging genes and samples with too many missing values.....step 1

  ### 2.2 表达矩阵的缺失数据的剔除，主要从样品和蛋白两个维度进行，对于蛋白就是从样品和蛋白两个维度进行阈值的筛选，超过缺失值筛选阈值的样品和蛋白就要剔除掉！
  if (!gsp$allOK) {
    # Optionally, print the protein and sample names that were removed:
    if (sum(!gsp$goodGenes)>0)
      printFlush(paste("Removing proteins:", paste(names(trusted_prot_clean_data_t)[!gsp$goodGenes], collapse = ", ")));
    if (sum(!gsp$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(trusted_prot_clean_data_t)[!gsp$goodSamples], collapse = ", ")));
    # Remove the offending proteins and samples from the data:
    trusted_prot_clean_data_t <- trusted_prot_clean_data_t[gsp$goodSamples, gsp$goodGenes]
  }

  ### 2.0 表达矩阵数据缺失值处理NA或NaN，R包Mice
  if(any(is.na(trusted_prot_clean_data_t))){

    # 后续添加缺失值填充流程
    if(is.null(NAFillingMethod)){
      print("默认使用最小值进行替换")
      trusted_prot_clean_data_t[is.na(trusted_prot_clean_data_t)] <- min(trusted_prot_clean_data_t[!is.na(trusted_prot_clean_data_t)])
    }else if(NAFillingMethod == "min"){
      trusted_prot_clean_data_t[is.na(trusted_prot_clean_data_t)] <- min(trusted_prot_clean_data_t[!is.na(trusted_prot_clean_data_t)])
    }else{
      print("请填写正确缺失值处理方式")
    }

  }else{
    print("表达数据中无空值，跳过缺失值填充阶段")
  }

  ## 1.0 保存预处理数据
  setwd(paste0(wd,"/",projectname,"/1.Data_Processing/1.1Data/"))
  savexlsx5(data = trusted_prot_clean_data_t,
            filename = "ExprData_nafill.xlsx",
            sheet = "ExprData_nafill")

  ### 2.4 保存表达矩阵的样品聚类矩阵数据（目录：WGCNA_P/1.Data_Processing/2.Sample_filtered）
  setwdfile(paste0(wd,"/",projectname,"/1.Data_Processing/1.2Sample_filtered/"),force = F)

  dist_matrix <- dist(trusted_prot_clean_data_t)
  sampletree <- hclust(dist_matrix, method = "average")
  dist_matrix <- as.matrix(dist_matrix)
  savexlsx5(data = dist_matrix,
            filename = "Distance_matrix.xlsx",
            sheet = "Distance_matrix")

  ### 2.6 绘制样品距离树距离点图（点图型的，可能有的客户需要，所以暂时注释掉）
  # height <- sort(sampletree$height, decreasing = TRUE)
  # plot(height, cex = 1.5)
  # h <- height[2] # 提取排序第二的样品Height值
  # abline(h = h, col = "red")

  ### 2.7 绘制划线样品距离树
  #### PDF
  pdf(file = "SampleTree.pdf", width = dim(dist_matrix)[1]/5, height = 12)
  height <- sort(sampletree$height, decreasing = TRUE)
  h <- height[2] # 提取排序第二的样品Height值
  plot(sampletree,
       main = " ",
       sub="", xlab="", col = "Blue", lwd = 2, cex = 1) # 绘制样品距离树图
  # abline(h = h+1, col = "red") # 划线
  dev.off()
  #### PNG
  png("SampleTree.png", width = 100*dim(dist_matrix)[1]/5, height = 1200)
  height <- sort(sampletree$height, decreasing = TRUE)
  h <- height[2] # 提取排序第二的样品Height值
  plot(sampletree,
       main = " ",
       sub="", xlab="", col = "Blue", lwd = 4, cex = 2)
  # abline(h = h+1, col = "red") # 画线
  dev.off()

  ### 2.8 样品平均表达条形图展示(3.0版新增分析结果)
  meanExpressionByArray <- apply(trusted_prot_clean_data_t, 1, mean, na.rm=T)
  NumberMissingByArray <- apply( is.na(data.frame(trusted_prot_clean_data_t)),1, sum)

  #### PDF
  pdf("Mean expression across samples.pdf", width = dim(dist_matrix)[1]/5, height = 12)
  barplot(meanExpressionByArray,
          xlab = "Sample", ylab = "Mean expression",
          # main ="Mean expression across samples",
          col = "red",
          names.arg = names(meanExpressionByArray), cex.names = 0.5, cex.axis = 1.1,
          axes = T,
          axisnames = T,
          axis.lty = 2,
          beside = T,
          width = 1)
  dev.off()

  #### PNG
  png("Mean expression across samples.png", width = 100*dim(dist_matrix)[1]/5, height = 1200)
  barplot(meanExpressionByArray,
          xlab = "Sample", ylab = "Mean expression",
          # main ="Mean expression across samples",
          col = "red",
          names.arg = names(meanExpressionByArray), cex.names = 0.5, cex.axis = 2,
          axes = T,
          axisnames = T,
          axis.lty = 2,
          beside = T,
          width = 1)
  dev.off()

  ### 2.9 手动处理离群样品：G22，可以考虑传参的方式去除，比如我想去除A样品，这边就把A传进来，然后剔除它；更优的办法是自动化去除离群样品
  if(is.null(removesample)){
    proExprData <- trusted_prot_clean_data_t
  }else{
    proExprData <- trusted_prot_clean_data_t[!(rownames(trusted_prot_clean_data_t) %in% removasample),]
  }

  print(paste0("proExprData数据样品数和", wgcna_project_type_cn,"个数分别为："))
  print(dim(proExprData)) # 查看proExprData数据大小

  ### 3.0 WGCNA分析蛋白表达矩阵Clean数据保存（去除离群样品后的）
  setwd(paste0(wd,"/",projectname,"/1.Data_Processing/1.1Data/"))
  savexlsx5(data = proExprData,
            filename = "ExprData_clean.xlsx",
            sheet = "ExprData_clean")

  ### 3.1 读取样品性状矩阵数据，绘制样品层次聚类树和样品性状热图。当前R工作路径在：./WGCNA_P/1.Data_Processing/1.1Data/中
  sampleTree2 <- hclust(dist(proExprData), method = "average") # Re-cluster samples

  ### 3.2 去除掉样品性状矩阵中的离群样品
  if(is.null(removesample)){
    newtraitData <- traitData
  }else{
    newtraitData <- traitData[!(rownames(traitData) %in% removasample),]
  }

  print("对比离群样品处理前和处理后的性状矩阵大小：")
  print(dim(traitData)) # 原始的形状数据大小
  print(dim(newtraitData)) # 原始的形状数据大小
  print(paste0("性状个数为：",dim(traitData)[2],"个！"))
  ### 3.3 导出去掉离群样品的性状数据
  savexlsx5(data = newtraitData,
            filename = "Sample_Traits_clean.xlsx",
            sheet = "Sample_Traits_clean")

  ### 3.4 将样品性状值转为颜色标记，此时的工作路径仍然是在：./WGCNA_P/1.Data_Processing/1.1Data/中
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  newtraitData2 <- as.data.frame(lapply(newtraitData, as.numeric), stringsAsFactors = F)
  rownames(newtraitData2) <- rownames(newtraitData)
  newtraitData2[is.na(newtraitData2)] <- 0
  traitColors <- numbers2colors(newtraitData2, signed = FALSE, naColor = "grey")

  ### 3.5 绘制样品表达矩阵和样品性状相关的树状图和热图Plot the sample dendrogram and the heatmap colors underneath.
  setwd(paste0(wd,"/",projectname,"/1.Data_Processing/1.2Sample_filtered/"))

  #### PDF
  pdf("Sample_dendrogram_and_trait_heatmap.pdf", width = dim(proExprData)[1]/5, height = dim(newtraitData)[2]/1.2)
  plotDendroAndColors(sampleTree2,
                      traitColors,
                      groupLabels = names(newtraitData),
                      main = " ",
                      autoColorHeight = TRUE,
                      addGuide = TRUE,
                      guideAll = TRUE)
  dev.off()
  #### PNG
  png("Sample_dendrogram_and_trait_heatmap.png", width = 100*dim(proExprData)[1]/5, height = 100*dim(newtraitData)[2]/1.2)
  plotDendroAndColors(sampleTree2,
                      traitColors,
                      groupLabels = names(newtraitData),
                      main = " ",
                      autoColorHeight = TRUE,
                      addGuide = TRUE,
                      guideAll = TRUE)
  dev.off()

  # 二、网络构建
  # ***
  # 蛋白共表达网络是无尺度（scale-free）的权重蛋白网络。无尺度特性，又称作无标度特性，是指网络的度分布满足幂律分布。幂律分布这一特性，正说明无尺度网络的度分布是呈集散分布：大部分的节点只有比较少的连接，而少数节点有大量的连接，这样的少数的多连接节点被称为共表达网络中的Hub节点。
  # 为了尽量满足无尺度网络分布前提条件，需要选择邻接矩阵权重参数power的取值。设定power值从1-30，分别计算其得到的网络对应的相关系数和网络的平均连接度。其中，相关系数越高（最大为1），网络越接近无尺度网络分布，但同时还需要保证一定的蛋白连接度，所以这个power值在相关系数足够大的同时保证基因的连接度较大。
  # 本次分析结果中所用的power值为：10（R^2 >= 8.5）

  ## 1.构建共表达网络Contrustion Co-expression Network
  ### 1.1 软阈值确定β_Power_selected
  powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
  ### 1.2 设置网络类型为：默认为type = unsigned
  # type <- "unsigned" #
  sft <- pickSoftThreshold(proExprData,
                           powerVector = powers,
                           networkType = type,
                           verbose=5)

  power <- sft$powerEstimate
  print(paste("=====================软阈值为：", power, "======================="))

  # 本次分析结果中所用的power值为：（R^2 >= 8.5）`r print(power)`（R^2 >= 8.5）

  ## 2.1 使用经验阈值
  nSamples <- nrow(proExprData)
  if (is.na(power)) {
    power = ifelse(nSamples < 20, ifelse(type == type, 9, 18),
                   ifelse(nSamples < 30, ifelse(type == type, 8, 16),
                          ifelse(nSamples < 40, ifelse(type == type, 7, 14),
                                 ifelse(type == type, 6, 12))
                   )
    )
  }

  ### 1.3 保存sft数据
  setwdfile(paste0(wd,"/",projectname,"/2.Network_Construction"),force = F)
  savexlsx1(data = sft$fitIndices,
            filename = "Soft_Threshold_identification.xlsx",
            sheet = "fitIndices")

  ### 1.4 自定义画布参数
  # op <- par(no.readonly = TRUE)

  #### 1.7.1 PDF
  pdf("Soft_Threshold_identification.pdf", width = 10, height = 5)
  par(mfrow = c(1,2))
  cex1 = 0.9
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",
       type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,
       cex=cex1,
       col="blue")
  abline(h=0.85, col="red") # RsquaredCut = 0.85, desired minimum scale free topology fitting index R^2.
  ## 软阈值Soft threshold与平均连通性Mean connectivity
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",
       ylab="Mean Connectivity",
       type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5],
       labels=powers,
       cex=cex1,
       col="blue")
  dev.off()

  #### 1.7.2 PNG
  png("Soft_Threshold_identification.png", width = 1000, height = 500)
  par(mfrow = c(1,2))
  cex1 = 0.9
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",
       type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,
       cex=cex1,
       col="blue")
  abline(h=0.85, col="red")
  ## 软阈值Soft threshold与平均连通性Mean connectivity
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",
       ylab="Mean Connectivity",
       type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5],
       labels=powers,
       cex=cex1,
       col="blue")
  dev.off()

  ## 2.2 无尺度拓扑的网络构建与评价
  #   评估结果如下图所示。近似直线关系（高R^2值，图中的右图）示出近似无标度拓扑。在大多数应用中，我们发现当选择高power来定义邻接矩阵时，无标度拓扑至少近似满足。
  ### 1.9检验选定的β值下记忆网络是否逼近 scale free
  if (type == "unsigned") {
    ADJ1_cor <- abs(WGCNA::cor(proExprData, use = "p" ))^power # 计算临接矩阵Adjacency
  } else if (type == "signed") {
    ADJ1_cor <- (0.5*(1 + WGCNA::cor(proExprData, use = "p" )))^power
  } else if (type == "signed hybrid") {
    ADJ1_cor <- WGCNA::cor(proExprData, use = "p" )^power
  } else if (type == "distance") {
    ADJ1_cor <- (1-(dist_matrix/max(dist_matrix))^2)^power
  } else {
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("主函数RUNWGCNA的参数type为空，或者选项不在给定的范围内，请确认后重新进行分析！")
  }

  ### 2.0计算连接度加权网络Calculates connectivity of a weighted network
  if (nrow(ADJ1_cor) < 5000) { # 蛋白少（<5000）的时候使用的代码
    k <- as.vector(apply(ADJ1_cor, 2, sum, na.rm = TRUE))
    head(k)
  } else if (nrow(ADJ1_cor) >= 5000) {
    k <- softConnectivity(datExpr = proExprData, power = power) # 蛋白多的时候使用的代码
  }

  #### 2.1.1 PDF
  pdf("Connectivity_and_Scale_free_topology_visualization.pdf", width = 10, height = 5)
  par(mfrow=c(1,2))
  hist(k, col = "grey", main = "Histogram of Connectivity with\n Power: 30") # 绘制Histogram of k
  scaleFreePlot(k, col = "blue", main = "Check Scale free topology (Power: 30)\n") # 绘制Scale free topology network if Power
  dev.off()

  #### 2.1.2 PNG
  png("Connectivity_and_Scale_free_topology_visualization.png", width = 1000, height = 500)
  par(mfrow=c(1,2))
  hist(k, col = "grey", main = "Histogram of Connectivity with\n Power: 30") # 绘制Histogram of k
  scaleFreePlot(k, col = "blue", main="Check Scale free topology (Power: 30)\n") # 绘制Scale free topology network if Power
  dev.off()

  ### 2.2经验power (无满足条件的power时选用)
  # 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
  # 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。

  ### 画布恢复默认变量
  par(op)

  ### 2.3构建共表达网络
  ### 一步法网络构建One-step network construction and module detection
  ### power: 上一步计算的软阈值
  ### maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)，这里是蛋白；
  ### 4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可以处理3万个
  ### 计算资源允许的情况下最好放在一个block里面。
  ### corType: pearson or bicornumericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
  ### saveTOMs：最耗费时间的计算，存储起来，供后续使用
  ### mergeCutHeight: 合并模块的阈值，越大模块越少
  nProteins <- ncol(proExprData)
  print(paste0(wgcna_project_type_cn, "总数为：",nProteins))
  nSamples <- nrow(proExprData)
  print(paste0("样品总数为：",nSamples))
  print("========================================================================")
  corFnc <- ifelse(corType == "pearson", WGCNA::cor, bicor)
  ### 2.4对二元变量，如样本性状信息计算相关性时，或基因(这里是蛋白）表达严重依赖于疾病状态时，需设置下面参数
  maxPOutliers <- ifelse(corType == "pearson",1,0.05)
  ### 2.5关联样品性状的二元变量时，设置
  robustY <- ifelse(corType == "pearson",T,F)

  ### 2.6Automatic network construction and module detection
  # This function performs automatic network construction and module detection on large expression datasets in a block-wise manner.
  if (!is.null(MinModuleSize)) { # 函数参数MinModuleSize不为空，即用户设置了参数的选项
    net <- blockwiseModules(proExprData, # Expression data
                            power = power,
                            maxBlockSize = nProteins,
                            TOMType = type,
                            minModuleSize = MinModuleSize,
                            reassignThreshold = 1e-6,
                            mergeCutHeight = 0.15,
                            numericLabels = TRUE,
                            pamRespectsDendro = FALSE,
                            saveTOMs = TRUE,
                            corType = corType,
                            maxPOutliers = maxPOutliers,
                            loadTOMs = TRUE,
                            saveTOMFileBase = "ExpressMatrix_.tom", verbose = 3)
  } else { # MinModuleSize参数的选项为空
    mf10.1 <- function(x) { # 自定义函数根据表达矩阵预估模块数量
      y <- round(x*(1/5))
      return(y)
    }
    MinModuleSize <- ncol(proExprData) %>% mf10.1
    net <- blockwiseModules(proExprData, # Expression data
                            power = power, # soft-thresholding
                            maxBlockSize = nProteins,
                            TOMType = type,
                            minModuleSize = MinModuleSize,
                            reassignThreshold = 1e-6,
                            mergeCutHeight = 0.15,
                            numericLabels = TRUE,
                            pamRespectsDendro = FALSE,
                            saveTOMs = TRUE,
                            corType = corType,
                            maxPOutliers = maxPOutliers,
                            loadTOMs = TRUE,
                            saveTOMFileBase = "ExpressMatrix_.tom", verbose = 3,
                            indent = 4)
  }

  ### 2.7根据模块中蛋白数目的多少，降序排列，依次编号为1-最大模块数、0-grey表示未分入任何模块的蛋白
  print("========================================================================")
  print(paste0(wgcna_project_type_cn,"共表达网络中的识别的模块如下："))
  print(table(net$colors))
  print("MEsgrey(ME0)算在内，总共鉴定到的模块数量为：")
  print(ncol(net$MEs))
  print(names(net$MEs))

  # 三、模块识别及可视化分析 #

  ## 6.模块识别和可视化Module_Identification_and_Visualization
  ### 此时的工作目录切换到./WGCNA_P/3.Module_Identification_and_Visualization
  setwdfile(paste0(wd,"/",projectname,"/3.Module_Identification_and_Visualization/"),force = F)

  ### 层级聚类树展示各个模块，灰色的为未分类到模块的蛋白
  moduleLabels <- net$colors # Convert labels to colors for plotting

  ### 保存moduleLabels为本地文件
  moduleLabels_number <- as.data.frame(moduleLabels)
  moduleLabels_number <- data.frame(rownames(moduleLabels_number), moduleLabels_number[,1])
  colnames(moduleLabels_number) <- c(wgcna_project_type_name,"moduleLabels")

  ### 数字代码到颜色的转换，获得moduleColors对象
  moduleColors <- labels2colors(moduleLabels)
  moduleLabels_Colors <- as.data.frame(moduleColors)

  ### 合并为一个新的数据框：ModuleColorsLabels
  ModuleColorsLabels <- cbind(moduleLabels_number, moduleLabels_Colors)

  ### 保存ModuleColorsLabels为xlsx文件
  savexlsx1(data = ModuleColorsLabels,
            filename = "ModuleColorsLabels.xlsx",
            sheet = "ModuleColorsLabels")

  ### 合并模块Merge Dynamic
  MEDissThres <- 0.5
  merge_modules <- mergeCloseModules(proExprData, moduleColors, cutHeight = MEDissThres, verbose = 3)
  mergedColors <- merge_modules$colors # 合并后的颜色
  mergedMEs <- merge_modules$newMEs # 新模块的特征向量蛋白

  print("模块合并之后的各模块的名称为：")
  print(names(mergedMEs))

  #### PDF
  pdf("Module_Cluster_Dendrogram_Merged_dynamic.pdf", width = 13, height = 6)
  plotDendroAndColors(net$dendrograms[[1]],
                      cbind(moduleColors, mergedColors),
                      c(paste("Dynamic Tree Cut(", ncol(net$MEs),")"), paste("Merged dynamic(", ncol(mergedMEs), ")")),
                      dendroLabels = FALSE,
                      hang = 0.03,
                      addGuide = TRUE,
                      guideHang = 0.05)
  # cex.lab = 1.2, cex.rowText = 1.2, cex.colorLabels = 1.2, cex.dendroLabels = 0.8)
  dev.off()

  #### PNG
  png("Module_Cluster_Dendrogram_Merged_dynamic.png", width = 1300, height = 700)
  plotDendroAndColors(net$dendrograms[[1]],
                      cbind(moduleColors, mergedColors),
                      c(paste("Dynamic Tree Cut(", ncol(net$MEs),")"), paste("Merged dynamic(", ncol(mergedMEs), ")")),
                      dendroLabels = FALSE,
                      hang = 0.03,
                      addGuide = TRUE,
                      guideHang = 0.05)
  # cex.lab = 1.2, cex.rowText = 1.2, cex.colorLabels = 1.2, cex.dendroLabels = 1.2)
  dev.off()

  # 四、模块与样品、性状的关联分析

  ## 7.模块与额外信息进行关联分析Modules_with_External_Information_Associations_Analysis
  ### 7.1模块数据提取和保存Modules_data_Extract_and_Save
  MEs <- net$MEs
  MEs_col <- MEs # MEs和MEs_col是一样的，只是列名用颜色名称替换了对应的数字编码
  colnames(MEs_col) <- paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
  MEs_col <- orderMEs(MEs_col)

  ### 保存数据MEs为本地文件
  savexlsx5(data = MEs_col,filename = "net_MEs.xlsx",sheet = "net_MEs")

  ### 保存合并之后的模块为本地文件
  savexlsx5(data = mergedMEs,filename = "Module_Eigengenes.xlsx",sheet = "Module_Eigengenes")

  ### 导出Cytoscape网络数据
  ### 提取指定模块的蛋白
  probes <- colnames(proExprData) ## 我们例子里面的probe就是蛋白名

  ### 如果采用分步计算，或设置的blocksize>=总基因数，直接load计算好的TOM结果，否则需要再计算一遍，比较耗费时间
  TOM <- TOMsimilarityFromExpr(proExprData,
                               power = power,
                               corType = corType,
                               networkType = type)

  print("============================================")
  print("TOM矩阵-Topological overlap matrix的行列数：")
  print(dim(TOM))
  print("============================================")

  dimnames(TOM) <- list(probes, probes)

  ### 先导出总的Cytoscape数据
  ### Export the network into edge and node list files Cytoscape can read threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在cytoscape中再调整
  cyt <- exportNetworkToCytoscape(TOM,
                                  weighted = TRUE,
                                  threshold = 0,
                                  nodeNames = probes,
                                  nodeAttr = moduleColors,
                                  includeColNames = TRUE)

  ### 此时的工作目录切换到./WGCNA_P/3.Module_Identification_and_Visualization
  setwdfile(paste0(wd,"/",projectname,"/4.Modules_with_External_Information_Associations_Analysis/4.1Modules_Samples_Associations_Analysis"),force = F)

  ### Module Colors: MEblue	MEturquoise MEyellow MEbrown MEred MEgreen ###
  ###
  ################################################################################
  #@           循环处理各个鉴定到的模块基因、蛋白质或代谢物的表达情况           #@
  ################################################################################

  moduleLabels_Colors_unique <- unique(moduleLabels_Colors$moduleColors)
  no_need_module <- which(moduleLabels_Colors_unique == "grey")
  moduleLabels_Colors_unique_need <- moduleLabels_Colors_unique[-no_need_module]

  ### 自定义函数绘制模块特有蛋白的表达条形图：ModuleEigenproteinBarplot
  ModuleEigenproteinBarplot <- function(module = "blue") {
    # 1.模块颜色传参确定
    mecolor <- paste0("ME",module)
    # 2.提取相同颜色的模块特有蛋白数据
    color_MEs_data <- MEs_col %>% as.data.frame() %>% dplyr::select(mecolor)
    # 3.绘制图形
    ggplot(data = color_MEs_data, mapping = aes(x = rownames(color_MEs_data), y = color_MEs_data[,1])) +
      geom_bar(stat = "identity", colour = module, fill = module) +
      theme_bw() +
      xlab("") +
      ylab(paste0("Eigen", wgcna_project_type_en, " Expression")) +
      geom_text(label = signif(color_MEs_data[,1], digits = 1),
                size = (200/nrow(color_MEs_data)+1), check_overlap = T, col = "black",
                position = "identity", hjust = "middle", vjust = -1, fontface = "bold") +
      theme(axis.title.y = element_text(size=14), axis.text.x = element_text(angle = 90, size = nrow(color_MEs_data)/25, face = "bold"))
    ggsave(filename = paste0(module,"_ModuleEigen", wgcna_project_type_en, "_barplot.pdf"),
           width = nrow(color_MEs_data)/8, height = 7.5, units = "in",limitsize = FALSE)
    ggsave(filename = paste0(module,"_ModuleEigen", wgcna_project_type_en, "_barplot.png"),
           width = nrow(color_MEs_data)/8, height = 7.5, units = "in",limitsize = FALSE)
  }

  ### WGCNA分析clean表达矩阵数据
  pro_data_clean_new <- t(proExprData)
  pro_data_clean_new <- cbind(rownames(pro_data_clean_new), pro_data_clean_new)
  colnames(pro_data_clean_new)[1] <- wgcna_project_type_name

  ### 循环处理各模块数据
  for (m in moduleLabels_Colors_unique_need) { # 遍历创建颜色模块目录
    print(paste0("当前处理",m,"模块！"))
    print("=====================================================================")
    # 工作目录设置，设置到颜色模块中
    setwdfile(path =paste0(wd,"/",projectname,"/4.Modules_with_External_Information_Associations_Analysis/4.1Modules_Samples_Associations_Analysis/",m),force = T)

    Module_Color <- m
    inModule_Color <- (moduleColors == Module_Color)
    modProbes_Color <- probes[inModule_Color]
    Color_TOM <- TOM[modProbes_Color, modProbes_Color]
    mcindex <- which(moduleColors == m)
    color_moduleColors <- moduleColors[mcindex]

    colorModule_cyt <- exportNetworkToCytoscape(Color_TOM,
                                                weighted = TRUE,
                                                threshold = 0,
                                                nodeNames = modProbes_Color,
                                                nodeAttr = color_moduleColors)

    ### 第一部分：绘制模块特有蛋白的表达条形图 ###
    ModuleEigenproteinBarplot(module = m)

    ### 数据保存
    # savexlsx1(data = colorModule_cyt$edgeData,filename = "linkdata.xlsx",sheet = "edgeData")
    # savexlsx1(data = colorModule_cyt$nodeData,filename = "linkdata.xlsx",sheet = "nodeData")

    print("*********************************************************************")
    print(paste0("BARPLOT is OK! | ",m,"模块的特征",wgcna_project_type_cn,"表达条形图绘制完成！"))
    print("*********************************************************************")

    ### 第二部分：绘制模块的蛋白表达箱线图 ###
    color_module_proteins_label <- colorModule_cyt$nodeData
    colnames(color_module_proteins_label) <- c("nodeName", "altName", "nodeAttr")
    rownames(color_module_proteins_label) <- color_module_proteins_label$nodeName
    color_module_proteins_label$altName <- NULL # 去掉altName列

    ### 获取属于颜色模块中的蛋白表达矩阵数据：color_Module_data = Module ProteinSampleExpression Data
    color_Module_data <- merge.data.frame(x = color_module_proteins_label, y = pro_data_clean_new,
                                          by.x = "nodeName", by.y = wgcna_project_type_name, all.x = TRUE)
    ### 保存颜色模块的蛋白表达矩阵数据
    savexlsx1(data = color_Module_data,filename = paste0(m,"_Module_data.xlsx"),sheet = paste0(m,"_Module_data"))

    color_boxplot_data <- pro_data_clean_new[pro_data_clean_new[,wgcna_project_type_name] %in% color_Module_data$nodeName,]
    color_boxplot_data <- as.data.frame(color_boxplot_data)
    color_boxplot_data <- reshape2::melt(data = color_boxplot_data,
                                         id.vars = wgcna_project_type_name,
                                         variable.name = "Sample_name",
                                         value.name = "Expression")
    color_boxplot_data$Expression <- as.numeric(color_boxplot_data$Expression)

    savexlsx1(data = color_boxplot_data,filename = paste0(m,"_module_barplot_data.xlsx"),sheet = paste0(m,"_module_barplot_data"))

    ### 不分组绘图
    pp <- ggplot(data = color_boxplot_data,
                 aes(x = Sample_name, y = Expression)) +
      geom_boxplot(col = "black", fill = m) +
      theme_bw() +
      xlab("") +
      ylab(paste0(wgcna_project_type_en, " Expression")) +
      theme(axis.title.y = element_text(size=14),
            plot.margin = unit(rep(0.6,4),'lines'),
            aspect.ratio = .6,
            axis.text.x = element_text(angle = 90,
                                       size = nSamples/25,
                                       face = "bold"))

    ggsave(plot = pp,filename = paste0(m,"_Module_boxplot.pdf"), width = nSamples/8, height = 7.5, units = "in",limitsize = FALSE)
    ggsave(plot = pp,filename = paste0(m,"_Module_boxplot.png"), width = nSamples/8, height = 7.5, units = "in",limitsize = FALSE)

    ### 分组绘图：bygroup
    # Sample_Group <- as.factor(color_boxplot_data$Sample_group)
    # ggplot(data = color_boxplot_data,
    #        aes(x = Sample_name, y = Expression, colour = Sample_Group)) +
    #   geom_boxplot() +
    #   theme_bw() +
    #   xlab("") +
    #   ylab(paste0(wgcna_project_type_en, " Expression")) +
    #   theme(axis.title.y = element_text(size=14),
    #         axis.text.x = element_text(angle = 90,
    #                                    size = nSamples/25,
    #                                    face = "bold"))
    # ggsave(filename = paste0(m,"_Module_boxplot_group.pdf"), width = nSamples/8, height = 7.5, units = "in",limitsize = FALSE)
    # ggsave(filename = paste0(m,"_Module_boxplot_group.png"), width = nSamples/8, height = 7.5, units = "in",limitsize = FALSE)

    print("*********************************************************************")
    print(paste0("BOXPLOT is OK! | ",m,"模块的", wgcna_project_type_cn,"表达箱线图绘制完成！"))
    print("*********************************************************************")

    ### 第三部分：绘制模块的蛋白表达聚类热图 ###
    color_Module_heatmap_data <- color_Module_data %>% as.data.frame() %>% dplyr::select(-nodeAttr)

    rownames(color_Module_heatmap_data) <- color_Module_heatmap_data$nodeName
    color_Module_heatmap_data$nodeName <- NULL
    color_Module_heatmap_data_new <- as.data.frame(lapply(color_Module_heatmap_data, as.numeric))
    rownames(color_Module_heatmap_data_new) <- rownames(color_Module_heatmap_data)
    ### PDF ###
    auto_heatmap(data = color_Module_heatmap_data_new,
                 type = c("png","pdf"),
                 name = paste0(m,"_Module_heatmap"),
                 cluster_rows = T,
                 cluster_cols = T,
                 display_numbers = F,
                 number_format = "%.3f",
                 show_colnames = T,
                 show_rownames = T,
                 angle_col = "45")

    print("*********************************************************************")
    print(paste0("HEATMAP is OK! | ",m,"模块的", wgcna_project_type_cn,"表达聚类热图绘制完成！"))
    print("*********************************************************************")
    print(paste0(m,"模块的",wgcna_project_type_cn,"表达情况分析完成！"))
    print("=====================================================================")
  }

  ################################################################################
  ################################################################################

  ### Iterative garbage collection.
  collectGarbage()

  ## 4.2 模块与性状关联分析

  ### 4.2.1 模块与性状的相关性及显著性热图

  ### 7.3模块与性状关联分析Modules_Traits_Associations_Analysis
  ### 切换工作目录到./WGCNA_P/4.Modules_with_External_Information_Associations_Analysis/4.2Modules_Traits_Associations_Analysis
  setwdfile(paste0(wd,"/",projectname,"/4.Modules_with_External_Information_Associations_Analysis/4.2Modules_Traits_Associations_Analysis/"),force = F)

  if (corType == "pearson") {
    modTraitCor <- cor(MEs_col, newtraitData, use = "p")
    modTraitP <- corPvalueStudent(modTraitCor, nSamples)
  } else {
    modTraitCorP <- bicorAndPvalue(MEs_col, newtraitData, robustY=robustY) # robustY=robustY
    modTraitCor <- modTraitCorP$bicor
    modTraitP <- modTraitCorP$p
  }

  ### signif表示保留几位小数
  textMatrix <- paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
  dim(textMatrix) <- dim(modTraitCor)
  # sizeGrWindow(10, 5)

  ### 保存数据为本地文件
  savexlsx5(data = modTraitCor,filename = "Module-Trait_Cor_data.xlsx",sheet = "Module-Trait_Cor_data")
  savexlsx5(data = modTraitP,filename =  "Module-Trait_P-value_data.xlsx",sheet = "Module-Trait_P-value_data")

  ### 绘制性状与模块的相关性与显著性热图
  par(mfrow = c(1,1))
  ### PDF ###
  pdf(file = "Module-Trait_relationship.pdf", width = dim(modTraitCor)[2], height = dim(modTraitCor)[1]*1.5)
  WGCNA::labeledHeatmap(Matrix = modTraitCor,
                        xLabels = colnames(newtraitData),
                        yLabels = colnames(MEs_col),
                        cex.lab = 1,
                        cex.lab.x = 1.2,
                        ySymbols = colnames(MEs_col),
                        colorLabels = FALSE,
                        colors = blueWhiteRed(50),
                        textMatrix = textMatrix,
                        setStdMargins = TRUE,
                        cex.text = 0.8,
                        zlim = c(-1,1),
                        main = paste("Module-trait relationships heatmap"))
  dev.off()
  ### PNG ###
  png(file = "Module-Trait_relationship.png", width = 100*dim(modTraitCor)[2], height = 100*dim(modTraitCor)[1]*1.5)
  WGCNA::labeledHeatmap(Matrix = modTraitCor,
                        xLabels = colnames(newtraitData),
                        yLabels = colnames(MEs_col),
                        cex.lab = 1,
                        cex.lab.x = 1.2,
                        ySymbols = colnames(MEs_col),
                        colorLabels = FALSE,
                        colors = blueWhiteRed(50),
                        textMatrix = textMatrix,
                        setStdMargins = TRUE,
                        cex.text = 1.5,
                        zlim = c(-1,1),
                        main = paste("Module-trait relationships heatmap"))
  dev.off()

  par(op) # 恢复默认参数
  #
  ### 计算模块与蛋白的相关性矩阵
  if (corType=="pearson") {
    geneModuleMembership <- as.data.frame(cor(proExprData, MEs_col, use = "p"))
    MMPvalue <- as.data.frame(corPvalueStudent(
      as.matrix(geneModuleMembership), nSamples))
  } else {
    geneModuleMembershipA <- bicorAndPvalue(proExprData, MEs_col, robustY=robustY)
    geneModuleMembership <- geneModuleMembershipA$bicor
    MMPvalue <- geneModuleMembershipA$p
  }

  ### 计算性状与蛋白的相关性矩阵，只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。
  if (corType=="pearson") {
    geneTraitCor <- as.data.frame(cor(proExprData, newtraitData, use = "p"))
    geneTraitP <- as.data.frame(corPvalueStudent(as.matrix(geneTraitCor), nSamples))
  } else {
    geneTraitCorA <- bicorAndPvalue(proExprData, newtraitData, robustY=robustY)
    geneTraitCor <- as.data.frame(geneTraitCorA$bicor)
    geneTraitP <- as.data.frame(geneTraitCorA$p)
  }

  ### 4.2.2 筛选目标模块和性状

  ### 最后把模块蛋白显著性矩阵和模块Membership相关性矩阵联合起来,指定感兴趣模块进行分析
  ### 卡阈值进行筛选与模块有显著相关性的样本性状，并进行分析和绘图
  modNames <- substring(colnames(MEs_col), 3)

  # sizeGrWindow(10, 10)
  par(mfrow = c(1,1))

  ### 自定义函数ModuleTraitCorScatterPlot，绘制目标模块与性状的相关性点图
  ModuleTraitCorScatterPlot <- function(module = "yellow", pheno = "R1") {
    # 1.传入参数模块和性状
    # module <- "yellow"
    # pheno <- "R1"
    # 2.获取的目标数据列
    module_column <- match(module, modNames)
    pheno_column <- match(pheno, colnames(traitData))
    # 3.获取模块内的基因/蛋白/代谢物
    moduleGenes <- moduleColors == module
    # 4.绘制模块与性状相关性点图
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                       abs(geneTraitCor[moduleGenes, pheno_column]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste(wgcna_project_type_en, " significance for", pheno),
                       main = paste("Module membership vs. ", wgcna_project_type_en, " significance\n"),
                       cex.main = 1.2, cex.lab = 1.3, cex.axis = 1.2, col = module, abline = TRUE,
                       abline.color = "red", abline.lty = 1, cex = 2, pch = 19)
  }

  ## 获取全部模块和全部样品性状数据
  pvalue_target_pheno_list <- colnames(modTraitP)
  ## 循环遍历所有目标性状进行相关性分析 ##
  for (dd in pvalue_target_pheno_list) {
    setwdfile(paste0(wd,"/",projectname,"/4.Modules_with_External_Information_Associations_Analysis/4.2Modules_Traits_Associations_Analysis/",dd),force = T)

    samp_target_pheno_list <- row.names(modTraitP)[modTraitP[,dd] < 0.05]
    samp_target_pheno_list <- c(samp_target_pheno_list,row.names(modTraitCor)[modTraitCor[,dd] < -0.5 | modTraitCor[,dd] > 0.5])
    samp_target_pheno_list <- unique(samp_target_pheno_list)
    if(length(samp_target_pheno_list) == 0){
      print("注意：目标性状没有产生！")
      print("使用所有的目标性状进行分析绘图！")
      samp_target_pheno_list <- row.names(modTraitP)
    }

    for (mm in samp_target_pheno_list) {
      ### PDF ###
      pdf(file = paste0("Module_membership(",gsub("ME", "", mm),")_vs_significance(",dd,").pdf"), width = 10, height = 10)
      ModuleTraitCorScatterPlot(module = gsub("ME", "", mm), pheno = dd)
      dev.off()
      ### PNG ###
      png(file = paste0("Module_membership(",gsub("ME", "", mm),")_vs_significance(",dd,").png"), width = 1000, height = 1000)
      ModuleTraitCorScatterPlot(module = gsub("ME", "", mm), pheno = dd)
      dev.off()
    }
  }

  # 五、关键蛋白共表达网络和TOM矩阵的可视化分析
  ## 5.1 关键基因/蛋白/代谢物共表达网络分析
  ## 8.模块Hub关键蛋白的共表达网络和TOM矩阵的可视化Hub_Proteins_Co-Expression_Network_and_TOM_Visualization
  setwdfile(paste0(wd,"/",projectname,"/5.Hub_Network_&_TOM_Visualization/5.1Modules_Hub_Network/"),force = F)

  for (M in moduleLabels_Colors_unique_need) { # 遍历创建颜色模块目录
    print(paste0("当前处理",M,"模块！"))
    print("=====================================================================")

    Module_Color <- M
    inModule_Color <- (moduleColors == Module_Color)
    modProbes_Color <- probes[inModule_Color]
    Color_TOM <- TOM[modProbes_Color, modProbes_Color]
    mcindex <- which(moduleColors == M)
    color_moduleColors <- moduleColors[mcindex]
    colorModule_cyt <- exportNetworkToCytoscape(Color_TOM,
                                                weighted = TRUE,
                                                threshold = 0,
                                                nodeNames = modProbes_Color,
                                                nodeAttr = color_moduleColors)

    color_module_proteins_label <- colorModule_cyt$nodeData
    colnames(color_module_proteins_label) <- c("nodeName", "altName", "nodeAttr")
    rownames(color_module_proteins_label) <- color_module_proteins_label$nodeName
    color_module_proteins_label$altName <- NULL # 去掉altName列
    #
    color_Module_data <- merge.data.frame(x = color_module_proteins_label, y = pro_data_clean_new,
                                          by.x = "nodeName", by.y = wgcna_project_type_name, all.x = TRUE)
    color_Module_datExpr <- t(color_Module_data)
    color_Module_intramodularConnet <- intramodularConnectivity.fromExpr(color_Module_datExpr,
                                                                         colors = color_Module_data$nodeAttr,
                                                                         corFnc = "cor",
                                                                         corOptions = "use = 'p'",
                                                                         # weights = NULL,
                                                                         distFnc = "dist",
                                                                         distOptions = "method = 'euclidean'",
                                                                         networkType = type,
                                                                         power = power,
                                                                         scaleByMax = FALSE,
                                                                         ignoreColors = if (is.numeric(colors)) 0 else "grey",
                                                                         getWholeNetworkConnectivity = TRUE)

    ### 根据模块内连接度由高往低进行排序，取top50蛋白作为Hub Protein进行Hub蛋白进行网络图可视化
    IntramConnet_color <- color_Module_intramodularConnet$kTotal
    color_Module_Hub_protein <- data.frame(color_Module_data, IntramConnet_color)

    colnames(color_Module_Hub_protein)[1:2] <- c(wgcna_project_type_name,"Module_Color")
    color_Module_Hub_protein$IntramConnet_color <- as.numeric(color_Module_Hub_protein$IntramConnet_color)
    color_Module_Hub_protein <- color_Module_Hub_protein[order(-color_Module_Hub_protein$IntramConnet_color),] # 由高到低排序

    if(dim(color_Module_Hub_protein)[1] > 50){
      color_Module_Hub_protein_top50 <- color_Module_Hub_protein[1:50,]
    }else{
      color_Module_Hub_protein_top50 <- color_Module_Hub_protein
    }
    if(wgcna_project_type=="P")
   { ### 导出全部的模块蛋白连接度数据
    savexlsx1(data = color_Module_Hub_protein,
              filename = paste0(M,"_Module_Hub_proteins.xlsx"),
              sheet = paste0(M,"_Module_Hub_proteins"))
    ### 导出top50的蛋白连接数据
    savexlsx1(data = color_Module_Hub_protein_top50,
              filename = paste0(M,"_Module_Hub_proteins_top50.xlsx"),
              sheet = paste0(M,"_Module_Hub_proteins_top50"))
   }else if(wgcna_project_type=="M")
   {
     ### 导出全部的模块代谢物连接度数据
    savexlsx1(data = color_Module_Hub_protein,
              filename = paste0(M,"_Module_Hub_Metabolites.xlsx"),
              sheet = paste0(M,"_Module_Hub_Metabolites"))
    ### 导出top50的代谢物连接数据
    savexlsx1(data = color_Module_Hub_protein_top50,
              filename = paste0(M,"_Module_Hub_Metabolites_top50.xlsx"),
              sheet = paste0(M,"_Module_Hub_Metabolites_top50"))
   }
    print("*********************************************************************")
    print(paste0("HUB数据保存完毕！ | ",M,"模块的关键",wgcna_project_type_cn,"全部和top50的数据保存完成！"))
    print("*********************************************************************")

    ### 绘制关键基因/蛋白/代谢物共表达网络图
    colorModule_cyt_df <- colorModule_cyt$edgeData
    colorModule_cyt_df[,4:6] <- NULL
    color_Module_Hub_protein_top50_ID <- color_Module_Hub_protein_top50[,1:2]
    color_Hub_Network <- merge.data.frame(color_Module_Hub_protein_top50_ID, colorModule_cyt_df, by.x = wgcna_project_type_name, by.y = "fromNode", all.y = TRUE)
    color_Hub_Network_new <- color_Hub_Network[complete.cases(color_Hub_Network),]
    color_Hub_Network_new <- color_Hub_Network_new[order(-color_Hub_Network_new$weight),]
    color_Hub_Network_plot <- color_Hub_Network_new %>%
      dplyr::select(wgcna_project_type_name, toNode, weight, Module_Color) # 数据框列的位置的重新排列
    if(dim(color_Hub_Network_plot)[1] > 50){
      color_Hub_Network_plot <- color_Hub_Network_plot[1:50,] # 取前50个蛋白
    }
    color_Hub_Net_plot_ig <- graph_from_data_frame(color_Hub_Network_plot, directed = F) # 参数directed为真，则有箭头，为假则无箭头方向
    ### PDF ###
    pdf(file = paste0(M,"_Hub_Network.pdf"), width = 10, height = 10)
    plot(color_Hub_Net_plot_ig,
         layout = layout_with_kk,
         # layout = layout_in_circle,
         vertex.color = M,
         vertex.size = 10,
         vertex.label.cex = 0.9,
         vertex.label.dist = 1.5,
         vertex.label.color = "black",
         main = paste0("Module: ", M," of Hub ",wgcna_project_type_en))
    dev.off()
    ### PNG ###
    png(file = paste0(M,"_Hub_Network.png"), width = 1000, height = 1000)
    plot(color_Hub_Net_plot_ig,
         layout = layout_with_kk,
         # layout = layout_in_circle,
         vertex.color = M,
         vertex.size = 10,
         vertex.label.cex = 0.9,
         vertex.label.dist = 1.5,
         vertex.label.color = "black",
         main = paste0("Module: ", M," of Hub ", wgcna_project_type_en))
    dev.off()
    #
    print(paste0("共表达网络图绘制完成！ | ", M,"模块的关键",wgcna_project_type_cn,"共表达网络分析和绘图完成！"))
    print("=====================================================================")
  }

  ## 5.2 TOM矩阵的可视化分析
  ### 5.2.1 所有基因/蛋白/代谢物的TOM矩阵可视化
  setwdfile(paste0(wd,"/",projectname,"/5.Hub_Network_&_TOM_Visualization/5.2TOM_Visualization/"),force = F)

  TOM <- TOMsimilarityFromExpr(proExprData,
                               power=power,
                               corType=corType,
                               networkType=type)
  # load(net$TOMFiles[1], verbose = T) # 导入计算好的TOM
  TOM <- as.matrix(TOM) # 数据类型转换
  probes <- colnames(proExprData)
  dimnames(TOM) <- list(probes, probes)
  dissTOM <- 1-TOM
  plotTOM <- dissTOM^7
  diag(plotTOM) <- NA

  ### 可视化dissTOM
  if (dim(dissTOM)[1] >= 500) { # 基因、蛋白或代谢物的数量大于等于1000个，执行此代码
    print("TOM加权相关性矩阵太大，无法全部绘制，所以只提供随机400个蛋白的")
    ### 5.2.2 随机选取400个蛋白的small-TOM可视化

    nSelect <- 500
    # For reproducibility, we set the random seed
    set.seed(100) # 设定随机种子
    select <- sample(nProteins, size = nSelect)
    selectTOM <- dissTOM[select, select]
    selectTree <- hclust(as.dist(selectTOM), method = "average") # 构建蛋白分层聚类树
    selectColors <- moduleColors[select]
    ## Open a graphical window，设置画布大小
    sizeGrWindow(10,10)
    # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing the color palette; setting the diagonal to NA also improves the clarity of the plot
    plotDiss <- selectTOM^7
    diag(plotDiss) <- NA

    #### 绘制随机选择的TOM矩阵（400）
    # ## PDF Save
    pdf("Network heatmap (400).pdf", width = 10, height = 10)
    TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap (400)")
    dev.off()
    # ## PNG Save
    png("Network heatmap (400).png", width = 1200, height = 1200)
    TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap (400)")
    dev.off()

  } else { # 基因、蛋白或代谢物的数量小于1000个，执行此代码
    ### PDF ###
    pdf(file = "Network_Heatmap_All.pdf", width = 10, height = 10)
    TOMplot(dissTOM,
            net$dendrograms,
            moduleColors,
            main = "Network heatmap (All)",
            terrainColors = T,
            setLayout = T)
    dev.off()
    ### PNG ###
    png(file = "Network_Heatmap_All.png", width = 1000, height = 1000)
    TOMplot(dissTOM,
            net$dendrograms,
            moduleColors,
            main = "Network heatmap (All)",
            terrainColors = T,
            setLayout = T)
    dev.off()

  }

  ################################################################################
  setwd(wd)
  print(paste0("当前工作路径为：",getwd()))
  ################################################################################

  ## 10.脚本耗时
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  print("脚本运行时间为：")

  print(proc.time() - pt)
  print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

  print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
  print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
}

