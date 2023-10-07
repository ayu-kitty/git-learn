#! /opt/conda/bin/Rscript

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  # Step1: Set up Seurat object for WGCNA
  # Step2: Construct metacells
  # Step3: Co-expression network analysis
  #    3.1 Set up the expression matrix
  #    3.2 Select soft-power threshold
  #    3.3 Construct co-expression network
  #    3.4 Optional: inspect the topoligcal overlap matrix (TOM)
  # Step4: Module Eigengenes and Connectivity
  #    4.1 Compute harmonized module eigengenes
  #    4.2 Compute module connectivity
  #    4.3 Compute hub gene signature scores
  # Step5: Module Network
  # Step6: diffGroup in Module
  # Step7: ModuleTraitCorrelation
  # Step8: Module Enrichment
  
  #================================================================================
  # package loading
  #================================================================================
  suppressWarnings({
    # single-cell analysis package
    suppressPackageStartupMessages(library(Seurat))
    # plotting and data science packages
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(cowplot))
    suppressPackageStartupMessages(library(patchwork))
    # co-expression network analysis packages:
    suppressPackageStartupMessages(library(WGCNA))
    suppressPackageStartupMessages(library(hdWGCNA))
    suppressPackageStartupMessages(library(harmony))
    suppressPackageStartupMessages(library(Cardinal))
    
    # gene enrichment packages
    suppressPackageStartupMessages(library(enrichR))
    suppressPackageStartupMessages(library(GeneOverlap))
    # others
    suppressPackageStartupMessages(library(openxlsx))
    suppressPackageStartupMessages(library(scales))
    suppressPackageStartupMessages(library(igraph))
    suppressPackageStartupMessages(library(ggforestplot))
    suppressPackageStartupMessages(library(ggrepel))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(reshape2))
    suppressPackageStartupMessages(library(R.utils))
    suppressPackageStartupMessages(library(paletteer))
    suppressPackageStartupMessages(library(gridExtra))
    suppressPackageStartupMessages(library(pheatmap))
    
    suppressPackageStartupMessages(library(optparse))
    
  })
  # ============= command line parameters setting=============================
  option_list = list(
    # make_option(c("-p","--script_path"), type = "character", default = NULL,
    # help = "hdWGCNA该代码本身所在的路径，可以是绝对路径，也可以是相对路径。"),
    make_option(c("-i","--input"), type = "character", default = NULL,
                help = "表达量矩阵."),
    make_option(c("-g","--group"), type = "character", default = NULL,
                help = "Group文件, 第一列为样本列, 其它列应至少有Group列和Sample列, 空代和空转中Group代表实际样本的临床分组/或者细胞来源的分区或者Celltype, 而Sample代表来源的样本" ),
    make_option(c("-a","--annotation"), type = "character", default = NULL,
                help = "注释文件, 如果输入矩阵第一列为 mz,或者gene ID等, 需要根据注释文件将mz或者gene ID转换为富集分析用的 C00001 号,或者gene Name, 便于提供富集分析用列表" ),
    make_option(c("-l","--orders"), type = "character", default = NULL,
                help = "输入既可以是一个只包含一列的文件,也可以是[]内的字符: [ N-epi,OLK-epi,OSCC-nepi,OSCC-ca2 ].指定Group列元素的顺序,绘图则也会按照这个顺序展示. 比如: c('N-epi','OLK-epi','OSCC-nepi','OSCC-ca2')" ),
    make_option(c("-s","--species"), type = "character", default = NULL,
                help = "物种, 比如为: hsa, mmu ....." ),
    make_option(c("-n","--hubgeneN"), type = "double", default = 25,
                help = "Top N Hubgene,默认为25" ),
    make_option(c("-z","--minModuleSize"), type = "double", default = 20,
                help = "构建Module 和 TOM的时候, 允许Module的最小基因数目. 默认为20. " ),
    make_option(c("-d","--deepSplit"), type = "double", default = 4,
                help = "0~4之间的数字,数字越大,模块个数越多,模块内的基因数目越少.默认为4." ),
    make_option(c("-t","--networkType"), type = "character", default = "signed",
                help = "signed代表构建Network的时候会考虑方向. unsigned/signed hybrid/signed" ),
    make_option(c("-c","--cormethod"), type = "character", default = "pearson",
                help = "相关性算法，可选方法： pearson / spearman / kendall. " ),
    make_option(c("-r","--traiscol"), type = "character", default = "all",
                help = "计算相关性时候所用到的性状的列名。默认为all，即分组文件中的所有列。否则输入列名，以逗号分隔，比如 [ p1,p2,p3,p4,p5]。" ),
    make_option(c("-k","--kvalue"), type = "double", default = 20,
                help = "nearest-neighbors parameter。默认为20。注：当数据集比较小的时候，默认的kvalue、min_cells、max_shared这三个参数，需要调整数值更小，否则会报错。" ),  
    make_option(c("-b","--min_cells"), type = "double", default = 25,
                help = "Metapixel中最小的cells数目，需要大于kvalue。默认为25。" ),  
    make_option(c("-x","--max_shared"), type = "double", default = 10,
                help = " maximum number of shared cells between two metacells。默认为10。" ),
    make_option(c("-o","--outputdir"), type = "character", default = "hdWGCNA_report",
                help = " 生成文件夹结果的命名。" ),
    make_option(c("-e","--wgcnaname"), type = "character", default = "WGCNA",
                help = "WGCNA对象的命名, 比如如果分析的为肿瘤不同恶性程度的Group数据, 可以命名为Malignancy, 会影响部分绘图的命名。" )
  );
  opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript hdWGCNA.R -i fpkm.xls -g group.xls -a anno.xls -s hsa -l N-epi,OLK-epi,OSCC-nepi,OSCC-ca2");
  opt = parse_args(opt_parser);
  
  # ---------------------------------------------------------
  #  ---- source functions
  
  theme_set(theme_cowplot())  # using the cowplot theme for ggplot
  set.seed(12345)  # set random seed for reproducibility
  enableWGCNAThreads(nThreads = 8)  # optionally enable multithreading
  
  #script_path <- opt$script_path
  input <- opt$input  #"DLM20223345_neg_NepiOLK....txt"
  group <- opt$group  #"group.xls"
  annotation <- opt$annotation # "数据矩阵.xlsx"
  orders <- opt$orders # c("N-epi","OLK-epi","OSCC-nepi","OSCC-ca2")
  species <- opt$species  # "hsa"
  hub_n_genes <- opt$hubgeneN  #25
  minModuleSize <- opt$minModuleSize # 20 # module的最小基因数
  deepSplit <- opt$deepSplit # 4, deepSplit: 0~4之间,值越大,module数目越多
  networkType <- opt$networkType # "signed"  # unsigned/signed hybrid/signed
  wgcna_name <- opt$wgcnaname # "Malignancy"
  cormethod <- opt$cormethod # 相关性算法 ： pearson, spearman, kandall
  traiscol <- opt$traiscol
  kvalue <- opt$kvalue
  min_cells <- opt$min_cells
  max_shared <- opt$max_shared
  dirprefix <- opt$outputdir
  
  #source(paste0(script_path,"/hdWGCNA_function.R"))
  suppressPackageStartupMessages(library(lmbio))
  ####################################################################################################
  ###-----------------------------------  step 0 : 数据读取 ---------------------------------------###
  # 不建议读原始的sscc的rds, 有NA（稀疏矩阵）
  print(">>>>>>>>>>>>>>>>>>>>>>>>>> reading matrix and group file.")
  wd <- getwd()
  
  # 读取矩阵文件，有些index是以数字开始(比如mz作为index)，如果是，前面会加上index-
  df_matrix <- read.csv(input,sep="\t",quote="",header=T,row.names=1,check.names=F,stringsAsFactors=FALSE, colClasses=c("character"))
  if(length((as.numeric(rownames(df_matrix)))) > 0){
    prefix_index <- "index-"
    rownames(df_matrix) <- paste0(prefix_index,rownames(df_matrix))
  }else{
    prefix_index <- "NONE"
  }
  
  df_group <- read.csv(group,sep="\t",quote="",header=T,row.names=1,check.names=F)
  
  if(is.null(orders) || orders == "None"){
    orders <- unique(df_group$Group)
  }else if(!is.null(orders) && grepl(",", orders) ){
    orders <- strsplit(orders, ",")[[1]]
    df_group$Group <- factor(df_group$Group, levels = orders)
    df_group <- df_group[order(df_group$Group), ]
  }else if(!is.null(orders) && !grepl(",", orders)){
    orders <- readLines(orders)
  }
  
  select_orders4net <- orders  #c("N-epi","OLK-epi","OSCC-nepi","OSCC-ca2")  ## 用于创建seurat_obj所用到的分组
  
  df_matrix <- df_matrix[,rownames(df_group)] ## 对matrix重新按照Group信息排序
  
  if(!identical(colnames(df_matrix),rownames(df_group))){
    stop("样本和分组中样本不对应,或者分组顺序不对应!")
  }
  
  print(">>>>>>>>>>>>>>>>>>>>>>>>>> create Seurat object file.")
  spectra_data <- df_matrix
  seurat_obj <- CreateSeuratObject(counts = spectra_data) #,meta.data = NULL ,assay = "Meta")
  
  table(is.na(seurat_obj@assays$RNA@data@x))  # 为什么SSCC的rds里面有空值 , 但umap的 rds 没有
  x <- seurat_obj@assays$RNA@counts
  sum(is.na(x))
  
  y <- seurat_obj@assays$RNA@scale.data
  sum(is.na(y))
  
  print(">>>>>>>>>>>>>>>>>>>>>>>>>> Step1: create raw dimension plot files.")
  seurat_obj <- FindVariableFeatures(object = seurat_obj)
  seurat_obj <- ScaleData(object = seurat_obj,features = VariableFeatures(object = seurat_obj))
  #seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj) )
  
  seurat_obj <- RunPCA(object = seurat_obj,dims=1:20)
  
  # 将所有的trais都加入到seurat_obj
  for (each_Trais in colnames(df_group)){
    print(each_Trais)
    seurat_obj[[`each_Trais`]] <- df_group[[`each_Trais`]]
  }
  # 运行 t-SNE
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:20,check_duplicates = FALSE)  
  ## 但最好找出dup的features及个数,并判定是否正常
  
  # 运行UMAP
  seurat_obj <- RunUMAP(object = seurat_obj,dims=1:20)
  
  # 运行harmony，原始用 'Group',有时候有报错。实际也应该根据样本'Sample'来去批次。
  seurat_obj <-  RunHarmony(object = seurat_obj,'Sample')
  
  ## 绘制 降维聚类图
  mycol_use <- mycol[1:length(unique(seurat_obj$Group))]
  
  ####################################################################################################
  ###-----------------------------------  step1.SetupForWGCNA  ---------------------------------------###
  # variable: use the genes stored in the Seurat object’s VariableFeatures.
  # fraction: use genes that are expressed in a certain fraction of cells for in 
  #           the whole dataset or in each group of    cells, specified by group.by.
  # custom: use genes that are specified in a custom list.
  
  print(">>>>>>>>>>>>>>>>>>>>>>>>>> step2: Setup WGCNA parameter.")
  dir.create(dirprefix)
  dirname1 <- paste0(dirprefix,"/1.WGCNA_setup")
  dir.create(dirname1)
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "variable", # the gene selection approach: fraction/variable(Seurat object’s VariableFeatures)/custom (gene_list)
    fraction = 0.01, #fraction of cells that a gene needs to be expressed in order to be included.
    #For example, fraction = 0.05 means that 5% of cells must express a gene (count > 0) for it to be included.
    wgcna_name = wgcna_name # the name of the hdWGCNA experiment
  )
  
  ################################################################################################
  ###--------------------------------  step2.Construct metacells  ---------------------------------###
  # Briefly, metacells are aggregates of small groups of similar cells originating from the same biological sample of origin. 
  # The k-Nearest Neighbors (KNN) algorithm is used to identify groups of similar cells to aggregate, and then the average or 
  # summed expression of these cells is computed, thus yielding a metacell gene expression matrix. 
  # We were originally motivated to use metacells in place of the original single cells because correlation network approaches 
  # such as WGCNA are sensitive to data sparsity.
  # The group.by parameter determines which groups metacells will be constructed in. 
  # We only want to construct metacells from cells that came from the same biological sample of origin, so it is critical to pass 
  # that information to hdWGCNA via the group.by parameter. Additionally, we usually construct metacells for each cell type separately. 
  # Thus, in this example, we are grouping by Sample and cell_type to achieve the desired result. 
  # 根据不同生物样本来进行metacell构建,同时demo中既基于不同样本,也基于不同的celltype构建了metacell,因此grouping by了Sample和cell_type
  # K值的选择: The number of cells to be aggregated k should be tuned based on the size of the input dataset, in general a lower 
  # number for k can be used for small datasets. We generally use k values between 20 and 75. The dataset used for this tutorial 
  # has 40,039 cells, ranging from 890 to 8,188 in each biological sample, and here we used k=25. The amount of allowable overlap 
  # between metacells can be tuned using the max_shared argument.
  
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("Group", "Sample"), # 也可以是c("cluster") .specify the columns in seurat_obj@meta.data to group by
    ident.group = 'Group', # set the Idents of the metacell seurat object
    reduction = 'harmony', # select the dimensionality reduction to perform KNN on,A dimensionality reduction stored in the Seurat object. Default = 'pca'
    k = kvalue, #20, # nearest-neighbors parameter
    min_cells = min_cells, #25, # min_cells must be greater then k value
    max_shared = max_shared, #10, # maximum number of shared cells between two metacells
    assay = 'RNA', # Assay to extract data for aggregation. Default = 'RNA'
    slot = "counts",
    mode = "average",
    target_metacells = 1000,
    max_iter = 5000,
    verbose = FALSE  # logical indicating whether to print additional information
  )
  
  # normalize metacell expression matrix:
  seurat_obj <- NormalizeMetacells(seurat_obj)
  metacell_obj <- GetMetacellObject(seurat_obj)
  
  seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
  seurat_obj <- RunPCAMetacells(seurat_obj,dims=1:20,features=VariableFeatures(seurat_obj))
  seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='Sample') ##?? seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='cluster')
  seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:20)
  
  RunTSNEMetacells <- function(seurat_obj, ...){
    seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- Seurat::RunTSNE(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, ...)
    seurat_obj
  }
  seurat_obj <- RunTSNEMetacells(seurat_obj, reduction='harmony', dims=1:20)
  
  plot_dimension(seurat_obj,DimPlotMetacells,'Metapixel','pca','PCA','Group','Sample',dirname1,mycol,12,6)
  #plot_dimension(seurat_obj,DimPlotMetacells,'Metapixel','harmony','Harmony','Group','Sample',"2.WGCNA_setup",mycol,12,6)
  plot_dimension(seurat_obj,DimPlotMetacells,'Metapixel','umap','UMAP','Group','Sample',dirname1,mycol,12,6)
  plot_dimension(seurat_obj,DimPlotMetacells,"Metapixel",'tsne','tSNE','Group','Sample',dirname1,mycol,12,6)
  
  
  table(is.na(seurat_obj@assays$RNA@scale.data))  # SSCC的rds里面有空值 , 但umap的 rds 没有
  x <- seurat_obj@assays$RNA@counts
  sum(is.na(x))
  
  y <- seurat_obj@assays$RNA@scale.data
  sum(is.na(y))
  
  ################################################################################################
  ###--------------------------------  step3.select soft power  ---------------------------------###
  
  dirname2 <- paste0(dirprefix,"/2.bestPower") #dirname3
  dir.create(dirname2)
  all_cluster <- unique(seurat_obj@meta.data$Group)
  
  ### ---------------   3.1 Set up the expression matrix
  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = select_orders4net, # the name of the group of interest in the group.by column ,可以选1个,也可以选多个list
    group.by='Group', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    use_metacells = TRUE, # A logical determining if we use the metacells (TRUE) or the full expression matrix (FALSE)
    assay = 'RNA', # using RNA assay
    slot = 'data' # using normalized data
  )
  
  ### ------------------3.2 Select soft-power threshold
  # Test different soft powers:
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),  # 这样power就是1:10和12:30(按照2的步长)
    use_metacells = TRUE,
    networkType = networkType, # you can also use "unsigned" or "signed hybrid"
    corFnc = "bicor", # Correlation function for the gene-gene correlation adjacency matrix.
    setDatExpr = FALSE, # logical flag indicating whether to run setDatExpr.
    group.by = NULL,
    group_name = NULL
  )
  
  # ---- plot the results:
  # 注: 源码中绘制soft power的时候,没有考虑slope负方向(一般可以乘以slope列數值的負方向，仅关注正值便可),
  #     因此可以手动重新计算[SFT.R.sq]列,加上方向
  
  # best_power选择的原则: 
  #     (1)一般认为大于0.8的话, 会符合无尺度分布, 所以如果有大于0.8的就选大于0.8的poewer最小值  
  #        也就是最早满足大于0.8; 
  #     (2)如果没有大于0.8的,则选择值最大的值。
  sign_x <- -sign(seurat_obj@misc[[wgcna_name]]$wgcna_powerTable$slope)
  seurat_obj@misc[[wgcna_name]]$wgcna_powerTable$new_SFT.R.sq <- seurat_obj@misc[[wgcna_name]]$wgcna_powerTable$SFT.R.sq * sign_x
  
  if(any(seurat_obj@misc[[wgcna_name]]$wgcna_powerTable$new_SFT.R.sq >= 0.8)){
    sele <- seurat_obj@misc[[wgcna_name]]$wgcna_powerTable[seurat_obj@misc[[wgcna_name]]$wgcna_powerTable$new_SFT.R.sq >= 0.8,]
    ord <- order(sele[,"new_SFT.R.sq"],decreasing = F)[1]
    best_power <- sele$Power[ord]
  }else{
    ord <- order(seurat_obj@misc[[wgcna_name]]$wgcna_powerTable$new_SFT.R.sq,decreasing = T)[1]
    best_power <- seurat_obj@misc[[wgcna_name]]$wgcna_powerTable$Power[ord]
  }
  
  plot_list <- PlotSoftPowers(seurat_obj,selected_power = best_power,plot_connectivity = TRUE)
  
  pdf(paste0(dirname1,"/SelectPower.pdf"))
  p <- wrap_plots(plot_list, ncol=2)
  print(p)
  dev.off()
  
  power_table <- GetPowerTable(seurat_obj)
  write_excel(power_table,"SelectPower",paste0(dirname1,"/SelectPower"))
  
  ### ------------------step3.3: best power
  seurat_obj <- ConstructNetwork(
    seurat_obj, 
    soft_power=best_power,  #### 可变,目前是手动输入
    setDatExpr=FALSE,
    tom_name = paste0("bestPower",best_power),
    tom_outdir = dirname2,
    overwrite_tom = TRUE,# name of the topoligical overlap matrix written to disk
    min_power = 3, # the smallest soft power to be selected if soft_power=NULL
    consensus = FALSE,  # flag indicating whether or not to perform Consensus network analysis
    wgcna_name = wgcna_name,
    blocks = NULL,
    maxBlockSize = 30000,
    randomSeed = 12345,
    corType = "pearson",
    consensusQuantile = 0.3,
    networkType = networkType,
    TOMType = "signed",
    TOMDenom = "min",
    scaleTOMs = TRUE,
    scaleQuantile = 0.8,
    sampleForScaling = TRUE,
    sampleForScalingFactor = 1000,
    useDiskCache = TRUE,
    chunkSize = NULL,
    deepSplit = deepSplit, # 0~4之间, 值越大, module越小数目越多. [ The higher the value (or if TRUE), the more and smaller clusters will be produced.] For method "hybrid", can be either logical or integer in the range 0 to 4. For method "tree", must be logical. In both cases, provides a rough control over sensitivity to cluster splitting. For the "hybrid" method, a finer control can be achieved via maxCoreScatter and minGap below.
    pamStage = FALSE,
    detectCutHeight = 0.99,  ### default = 0.995
    minModuleSize = minModuleSize,   ## Minimum Module size.
    mergeCutHeight = 0.2,  ### default = 0.2 模块合并时候的聚类树高度，值越大，合并的模块越多，最终生成的模块越少(测试了0.05~0.2, 结果一样,无变化,暂定default)
    saveConsensusTOMs = TRUE
  )
  
  ###  -----------------------  修改WGCNA的颜色, 自定义 ------------------------------------------
  
  # 更改 seurat_obj@misc[[wgcna_name]]$wgcna_net$colors
  raw_col <- seurat_obj@misc[[wgcna_name]]$wgcna_net$colors
  colour_raw2new_df <- colour_raw2new(raw_col,mycol_module)
  new_col <- change_colour(raw_col,colour_raw2new_df)
  names(new_col) <- names(raw_col)
  seurat_obj@misc[[wgcna_name]]$wgcna_net$colors <- new_col
  # 更改 seurat_obj@misc[[wgcna_name]]$wgcna_net$unmergedColors
  raw_unmerged_col <- seurat_obj@misc[[wgcna_name]]$wgcna_net$unmergedColors
  colour_raw2new_df <- colour_raw2new(raw_unmerged_col,mycol_module)
  new_unmerged_col <- change_colour(raw_unmerged_col,colour_raw2new_df)
  names(new_unmerged_col) <- names(raw_unmerged_col)
  seurat_obj@misc[[wgcna_name]]$wgcna_net$unmergedColors <- new_unmerged_col
  # 更改 seurat_obj@misc[[wgcna_name]]$wgcna_net$multiMEs[[wgcna_name]]$data 中的列名 MExxx
  raw_colnames <- colnames(seurat_obj@misc[[wgcna_name]]$wgcna_net$multiMEs[[wgcna_name]]$data)
  raw_colnames_x <- gsub("^ME","",raw_colnames)
  new_colnames <- change_colour(raw_colnames_x,colour_raw2new_df)
  new_colnames <- paste0("ME",new_colnames)
  colnames(seurat_obj@misc[[wgcna_name]]$wgcna_net$multiMEs[[wgcna_name]]$data) <- new_colnames
  # 更改 seurat_obj@misc[[wgcna_name]]$wgcna_modules$module 
  new_col1 <- change_colour(seurat_obj@misc[[wgcna_name]]$wgcna_modules$module,colour_raw2new_df)
  seurat_obj@misc[[wgcna_name]]$wgcna_modules$module <- new_col1
  seurat_obj@misc[[wgcna_name]]$wgcna_modules$module <- as.factor(seurat_obj@misc[[wgcna_name]]$wgcna_modules$module)
  # 更改 seurat_obj@misc[[wgcna_name]]$wgcna_modules$color 
  new_col2 <- change_colour(seurat_obj@misc[[wgcna_name]]$wgcna_modules$color,colour_raw2new_df)
  seurat_obj@misc[[wgcna_name]]$wgcna_modules$color <- new_col2
  
  net <- GetNetworkData(seurat_obj, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)
  
  
  pdf(paste0(dirname2,"/Dendrogram.pdf"))
  p <- PlotDendrogram(
    seurat_obj, 
    main='hdWGCNA Dendrogram',
    groupLabels =  , #"Module colors",
    wgcna_name = NULL,
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )  # 更多参数可以参考 WGCNA::plotDendroAndColors()
  print(p)
  dev.off()
  
  TOM <- GetTOM(seurat_obj)
  TOM_w <- TOM
  id <- rownames(TOM_w)
  TOM_w<-cbind(id,TOM_w)
  write_excel(TOM_w,"TOM",paste0(dirname2,"/bestPower",best_power,"_TOM"))
  
  dissTOM = 1-TOM# Transform dissTOM with a power to make moderately strong
  plotTOM = dissTOM^7
  
  pdf(paste0(dirname2,"/TOM-NetHeatmap.pdf"))
  TOMplot(
    plotTOM,
    net$dendrograms,
    Colors = net$colors,
    terrainColors = FALSE,
    main = "TOM - Network heatmap"
  )
  dev.off()
  
  save.image(file=paste0("","step2.RData"))
  
  ################################################################################################
  ###---------------------------  step4. Module Eigengenes and Connectivity  -----------------------------###
  # 会从四个角度计算： hMEs, MEs, scores, or average
  # 1. avergae：即模块特征基因的平均表达量
  # 2. MEs: 模块特征基因(MEs)是一种常用的度量来总结整个共表达模块的基因表达谱。
  #    简而言之，通过对包含每个模块的基因表达矩阵的子集执行主成分分析(PCA)来计算模块特征基因。每个PCA矩阵的第一个PC是MEs。
  # 3. hMEs：通过harmony进行批次矫正后的 MEs
  # 4. scores: Gene scoring analysis 基因评分分析是单细胞转录组学中一种流行的方法，用于计算一组基因的总体特征得分。
  #    hdWGCNA的ModuleExprScore()函数，用于计算每个模块给定数量的基因的基因分数，使用Seurat或UCell算法。
  #    基因评分是通过计算模块特征基因来总结模块表达的另一种方法。
  
  dirname3 <- paste0(dirprefix,"/3.ModuleEigengenes")
  dir.create(dirname3)
  dir.create(paste0(dirname3,"/individual"))
  
  ###------------  4.1 计算模块特征基因连接度，以及harmonized特征基因连接度 (对应后面算法的 MEs 以及 hMEs) 
  # 数据存储在 seurat_obj@misc[[`wgcna_name`]]$MEs 以及 seurat_obj@misc[[`wgcna_name`]]$hMEs   
  # hdWGCNA的ModuleEigengenes()函数可以计算single中的module eigengenes。
  # 同时可以通过Harmony进行MEs的 batch correction即批次矫正，以获取 harmonized module eigengenes (hMEs)
  # Compute harmonized module eigengenes. 
  # need to run ScaleData first or else harmony throws an error:  ##### 没试过harmony , 最好不要scale. 否则有na
  # seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
  seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
  
  # compute all MEs in the full single-cell dataset
  seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars="Sample"  #groups to harmonize by
  )
  
  ###-------- 4.2 Compute module connectivity 计算特征基因连接度 kMEs (对应后面算法的 kMEs)-------------------------
  # In co-expression network analysis, we often want to focus on the “hub genes”, those which are highly connected within each module. 
  # Therefore we wish to determine the eigengene-based connectivity, also known as kME, of each gene. 
  # hdWGCNA includes the ModuleConnectivity to compute the kME values in the full single-cell dataset, 
  # rather than the metacell dataset. This function essentially computes pairwise correlations between genes and module eigengenes. 
  # kME can be computed for all cells in the dataset, but we recommend computing kME in the cell type or group that was previously used to run ConstructNetwork.
  # 1. kME (eigengene connectivity) 特征基因连接度. 用于筛选Hubgene
  #    Hub genes are those that show most connections in the network as indicated by their high KME (eigengene connectivity) value
  
  # compute eigengene-based connectivity (kME):
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = 'Group', 
    group_name = select_orders4net
  )
  
  # rename the modules
  seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name = paste0(wgcna_name,"-M")
  )
  ### 4.3(1)计算Module中Hubgene平均表达量 (对应后面算法的 average )
  # 数据存储在 seurat_obj@misc[[`wgcna_name`]]$avg_modules   
  seurat_obj <- AvgModuleExpr(
    seurat_obj,
    n_genes = hub_n_genes
  )
  
  ### 4.3(2) Compute hub gene signature scores 计算模块内基因评分
  # 数据存储在 seurat_obj@misc[[`wgcna_name`]]$module_scores   
  
  
  # 报错处理
  # compute gene scoring for the top 25 hub genes by kME for each module
  # with Seurat method
  # 注意: 在这一步骤(seurat的utilities.R源码), 会使用sample函数对matrix进行抽样
  # 增加判断个数的输出
  
  
  n_genes <- hub_n_genes
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != 'grey']
  
  if(n_genes == "all"){
    gene_list <- lapply(mods, function(cur_mod){
      subset(modules, module == cur_mod) %>% .$gene_name
    })
  } else{
    gene_list <- lapply(mods, function(cur_mod){
      cur <- subset(modules, module == cur_mod)
      cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
        top_n(n_genes) %>% .$gene_name
    })
  }
  names(gene_list) <- mods
  
  # 计算bins后的个数, 并按照 70% 比例.
  # 在 seurat的AddModuleScore()函数中,会进行subsample操作: 先对所有的features进行bins操作, bins后对每个bins结果中抽样ctrl个features用作control样本.
  # AddModuleScore()函数的意义: 感兴趣的Hubgene基因抽出来，每一个细胞算一个这些基因表达的平均值，
  # 背景基因的平均值在于找每个基因的所在的bin，在该bin内随机抽取相应的ctrl个基因作为背景，
  # 最后所有的目标基因算一个平均值，所有的背景基因算一个平均值，两者相减就是该gene set 的score值
  
  nbin <-  24 ## 设定bins的个数
  print(length(seurat_obj@misc[[`wgcna_name`]]$datExpr[1,]))
  bins <- table(cut_number(x = runif(length(seurat_obj@misc[[`wgcna_name`]]$datExpr[1,])) , 
                           n = nbin, 
                           labels = FALSE, 
                           right = FALSE))
  x <- unique(bins)[1]
  
  # Computes a module score for each co-expression module using Seurat AddModuleScore or UCell.
  seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = n_genes,  # the number of genes to use for each module, ranked by kME. Setting n_genes = 'all' uses all of the genes in a module
    ctrl = x*0.7,  # Number of control features selected from the same bin per analyzed featurefloor(length(x)*0.7),  
    # ctrl应该小于length(x), 因为涉及到subsample, 这儿设定了 0.7比例
    # github上关于选择多少的建议: https://github.com/satijalab/seurat/issues/4558
    method='Seurat',
    nbin = 24,
    seed = 1,
    assay = NULL,
    name = 'Group',
    search = FALSE
  )
  
  # 排序
  order_x <- sort(colnames(seurat_obj@misc[[`wgcna_name`]]$MEs))
  
  seurat_obj@misc[[`wgcna_name`]]$MEs <- seurat_obj@misc[[`wgcna_name`]]$MEs[,order_x]
  seurat_obj@misc[[`wgcna_name`]]$hMEs <- seurat_obj@misc[[`wgcna_name`]]$hMEs[,order_x]
  
  order_y <- sort(colnames(seurat_obj@misc[[`wgcna_name`]]$module_scores))
  seurat_obj@misc[[`wgcna_name`]]$module_scores <- seurat_obj@misc[[`wgcna_name`]]$module_scores[,order_y]
  seurat_obj@misc[[`wgcna_name`]]$avg_modules <- seurat_obj@misc[[`wgcna_name`]]$avg_modules[,order_y]
  
  
  # 获取 MEs以及hMEs方法
  # harmonized module eigengenes:
  hMEs <- GetMEs(seurat_obj)
  
  # module eigengenes:
  MEs <- GetMEs(seurat_obj, harmonized=FALSE)
  
  # module scores
  module_score <- GetModuleScores(seurat_obj)
  
  # module average expression
  average_module <- GetAvgModuleExpr(seurat_obj)
  
  ###  绘图以及输出表格 
  # 1. 绘制所有features的kMEs图
  # We can visualize the genes in each module ranked by kME using the PlotKMEs function.
  pdf(paste0(dirname3,"/3.1_kME_allFeatures.pdf"),width=8,height=6)
  PlotKMEs(seurat_obj, 
           ncol=5,
           plot_widths = c(5, 2), # the relative width between the kME rank plot and the hub gene text
           n_hubs = hub_n_genes,
           text_size = 1.0,
           wgcna_name = NULL
  )
  dev.off()
  
  #####  make a featureplot of hMEs for each module
  reduction_type <- 'umap'
  # 模块相关性图
  ## average expression相关
  all_plot_1(seurat_obj,"average",FALSE,TRUE,reduction_type,dirname3,"/3.2_averageExpr_","3.2_Module_Corplot_averageExpr","averageExpr",TRUE)
  
  ## hMEs相关
  hMEs_yc <- "n"
  if(hMEs_yc == "y"){
    all_plot_1(seurat_obj,"hMEs",FALSE,TRUE,reduction_type,dirname3,"/3.2_hMEs_","3.2_Module_Corplot_hMEs","MEs",TRUE)
  }
  ## MEs相关
  all_plot_1(seurat_obj,"MEs",FALSE,TRUE,reduction_type,dirname3,"/3.2_MEs_","3.2_Module_Corplot_MEs","MEs",TRUE)
  ## Module scores相关
  error1 <- try({
    all_plot_1(seurat_obj,"scores",FALSE,TRUE,reduction_type,dirname3,"/3.3_ModuleScores_","3.2_Module_Corplot_scores","ModuleScores",TRUE)
  },silent=F)
  if ("try-error" %in% class(error1) ) {
    graphics.off()
    all_plot_1(seurat_obj,"scores",FALSE,TRUE,reduction_type,dirname3,"/3.3_ModuleScores_","3.2_Module_Corplot_scores","ModuleScores",FALSE)
  }
  
  ### ------------------------------------- 输出表格 ---------------------------------
  # 1.导出各模块所有的 gene
  all_genesinModule <- seurat_obj@misc[[`wgcna_name`]]$wgcna_modules
  write_excel(all_genesinModule,"wgcna_modules",paste0(dirname3,"/3.1_allFeatures_inModules"))
  
  # 获取 top gene 连接度 
  hub_df <- GetHubGenes(seurat_obj, n_hubs = hub_n_genes)
  write_excel(hub_df,"hubgene",paste0(dirname3,"/3.1_top",hub_n_genes,"_HubFeatures_inModules"))
  
  # 获取所有的gene
  #hub_all <- GetHubGenes(seurat_obj, n_hubs = length(seurat_obj@assays$RNA@counts[,1]))
  #write_excel(hub_all,"hubgene",paste0(dirname4,"/kME_allFeatures"))
  
  ### ------------------------------- 绘制 趋势图--------------------------------------
  
  # get hMEs from seurat object
  mods <- colnames(MEs); mods <- mods[mods != 'grey']
  
  ## 绘制小提琴图以及点图
  ## --------   average expression
  # add hMEs to Seurat meta-data:
  average_module_x <- average_module
  colnames(average_module_x) <- paste0("AVR_",colnames(average_module_x))
  average_module_x <- ScaleData(average_module_x)
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, average_module_x)
  #plot_VLNandDot(seurat_obj,wgcna_name,mods,modx,output1,output2,output3)
  plot_VLNandDot(seurat_obj,wgcna_name,hub_df,mods,"average","AVR_","3.3_Module2Cluster_DotPlot_averageExpr.pdf","_violin_averageExpr.pdf","3.3_violinAll_averageExpr.pdf",dirname3,mycol)
  
  ## --------   hMEs 
  # add hMEs to Seurat meta-data:
  hMEs_x <- hMEs
  colnames(hMEs_x) <- paste0("hMEs_",colnames(hMEs_x))
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, hMEs_x)
  #plot_VLNandDot(seurat_obj,wgcna_name,mods,modx,output1,output2,output3)
  plot_VLNandDot(seurat_obj,wgcna_name,hub_df,mods,"hMEs","hMEs_","3.3_Module2Cluster_DotPlot_hMEs.pdf","_violin_hMEs.pdf","3.3_violinAll_hMEs.pdf",dirname3,mycol)
  
  ## --------   MEs 
  # add hMEs to Seurat meta-data:
  MEs_x <- MEs
  colnames(MEs_x) <- paste0("MEs_",colnames(MEs_x))
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs_x)
  #plot_VLNandDot(seurat_obj,wgcna_name,mods,modx,output1,output2,output3)
  plot_VLNandDot(seurat_obj,wgcna_name,hub_df,mods,"MEs","MEs_","3.3_Module2Cluster_DotPlot_MEs.pdf","_violin_MEs.pdf","3.3_violinAll_MEs.pdf",dirname3,mycol)
  
  ## --------   module gene score
  # add hMEs to Seurat meta-data:
  module_score_x <- module_score
  colnames(module_score_x) <- paste0("MdScore_",colnames(module_score_x))
  module_score_x <- ScaleData(module_score_x)
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, module_score_x)
  #plot_VLNandDot(seurat_obj,wgcna_name,mods,modx,output1,output2,output3)
  plot_VLNandDot(seurat_obj,wgcna_name,hub_df,mods,"scores","MdScore_","3.3_Module2Cluster_DotPlot_module_score.pdf","_violin_module_score.pdf","3.3_violinAll_module_score.pdf",dirname3,mycol)
  
  #########  ---------------  绘制基于表达的 折线图和拟合线图 ---------------------
  
  ## plot Module Hubgene linePlot
  all_color <- paletteer_d( "ggsci::default_igv",n=51)
  all_color <- alpha(all_color,0.5)
  
  
  p_line <- NULL
  p_fitline <- NULL
  # All_hubgene <- NULL
  # M_group <- NULL
  
  plot_line_yn <- "n"
  
  module_type <- unique(hub_df$module)
  
  for (each in 1:length(module_type)){
    each_name <- module_type[each]
    each_module_hub_df <- hub_df[which(hub_df$module == each_name),]
    each_hub_df <- spectra_data[each_module_hub_df$gene_name,]
    index <- rownames(each_hub_df)
    each_hub_df_w <- cbind(index,each_hub_df)
    write_excel(each_hub_df_w,"module_hub",paste0(dirname3,"/individual/",each_name,"_Hubgene_table"))
    
    data_merged <-merge(df_group[,'Group',drop=F],t(each_hub_df),by.x='row.names',by.y='row.names')
    melt_data <- reshape2::melt(data_merged,id=c("Row.names","Group"))
    melt_data$value <- as.numeric(melt_data$value)
    row_mean <- aggregate(melt_data$value, by=list(melt_data$Group,melt_data$variable),mean)
    colnames(row_mean)<- c("Group","Feature","Value")
    fa <- select_orders4net
    row_mean$Group <- factor(row_mean$Group,levels=fa)
    length_period <- length(unique(row_mean$Group))
    raw_level <- levels(row_mean$Group)
    levels(row_mean$Group) <- c(1:length_period)
    
    # 不提供
    if(plot_line_yn == "y") {
      pdf(paste0(dirname3,"/individual/",each_name,"_hubLines_Expr.pdf"),width=10,height=6)
      p1 <- plotline(row_mean,"with","with",each_name,all_color)
      p11 <- plotline(row_mean,"without","without",each_name,all_color)
      print(p1)
      dev.off()
      
      if(each == 1){
        p_line[[1]] <- p11
      }else{
        p_line[[each]] <- p11
      }
    }
    pdf(paste0(dirname3,"/individual/",each_name,"_hubFitline_Expr.pdf")) #,width=8,height=6)
    p2 <- plot_smooth(row_mean,"#B09C8566","with","with",each_name,all_color)
    p22 <- plot_smooth(row_mean,"#B09C8566","without","without",each_name,all_color)
    
    print(p2)
    dev.off()
    if(each == 1){
      p_fitline[[1]] <- p22
    }else{
      p_fitline[[each]] <- p22
    }
  }
  
  ncol_n <- 3  ## 每行排几个
  
  if(plot_line_yn == "y") {
    pdf(paste0(dirname3,"/3.3_hubLines_all_Expr.pdf"),width=ncol_n*2,height=(2*(ceiling(2*length(module_type)/ncol_n-1))))
    grid.arrange(grobs =p_line, ncol = ncol_n)
    dev.off()
  }
  pdf(paste0(dirname3,"/3.3_hubFitline_all_Expr.pdf"),width=ncol_n*2,height=(2*(ceiling(2*length(module_type)/ncol_n-1))))
  grid.arrange(grobs =p_fitline, ncol = ncol_n)
  dev.off()
  
  save.image(file=paste0("","step3.RData"))
  
  ################################################################################################
  ###---------------------------  Step4:    Module Network  -----------------------------###
  # Step5: Module Network 
  # https://smorabit.github.io/hdWGCNA/articles/network_visualizations.html#individual-module-network-plots
  # 1.ModuleNetwork
  dirname4 <- paste0(dirprefix,"/4.ModuleNetwork")
  dir.create(dirname4)
  
  #pdf(paste0(dirname5,"/HubGeneNetwork.pdf"))
  HubGeneNetwork <- ModuleNetworkPlot(
    seurat_obj,
    mods = "all",
    plot_size = c(6, 6),
    wgcna_name = NULL,
    label_center = FALSE,
    edge.alpha = 0.25,
    vertex.label.cex = 1,
    vertex.size = 6,
    outdir = dirname4
  )
  #dev.off()
  
  #hubgene network
  pdf(paste0(dirname4,"/HubGeneNetwork.pdf"))
  HubGeneNetwork <- HubGeneNetworkPlot(
    seurat_obj,
    mods = 'all',  # Names of the modules to plot. If mods = "all", all modules are plotted.
    n_hubs = hub_n_genes,  # The number of hub genes to plot for each module.
    n_other= 5, # The number of non-hub genes to sample from each module
    sample_edges = TRUE,
    edge_prop = 0.75,
    return_graph = FALSE,
    edge.alpha = 0.25,
    vertex.label.cex = 0.5,
    hub.vertex.size = 4,
    other.vertex.size = 1,
    wgcna_name = NULL
  )
  print(HubGeneNetwork)
  dev.off()
  
  ################################################################################################
  ###--------------------Step5:  Differential module eigengene (DME) analysis-------------------------###
  ### DME analysis comparing two groups
  # cluster_select
  dirname5 <- paste0(dirprefix,"/5.diffGroup_in_Module")
  dir.create(dirname5)
  seurat_obj@meta.data$barcode <- rownames(seurat_obj@meta.data)
  
  for (cl in 1:length(all_cluster)){
    diff_cluster_select <- as.character(all_cluster[cl])
    
    group1 <- seurat_obj@meta.data %>% subset(Group == diff_cluster_select) %>% rownames
    group2 <- seurat_obj@meta.data %>% subset(Group != diff_cluster_select) %>% rownames
    
    DMEs <- FindDMEs(
      seurat_obj,
      barcodes1 = group1,
      barcodes2 = group2,
      test.use='wilcox',
      wgcna_name= wgcna_name,
      harmonized = TRUE,
      add_missing = FALSE,
      only.pos = FALSE,
      logfc.threshold = 0,
      min.pct = 0,
      verbose = FALSE,
      pseudocount.use = 0
    )
    
    colnames(DMEs) <- c("pvalue","avg_log2FC",diff_cluster_select,"Others_group","p_val_adj","module")
    
    pdf(paste0(dirname5,"/",diff_cluster_select,"-vs-otherGroups_DMEsLollipopPlot",".pdf"))
    
    DMEsLollipop <- PlotDMEsLollipop(
      seurat_obj, 
      DMEs, 
      wgcna_name = wgcna_name,
      pvalue = "pvalue",
      avg_log2FC = "avg_log2FC"
    )
    print(DMEsLollipop)
    dev.off()
    
    
    DMEs_w <- DMEs
    module <- rownames(DMEs_w)
    DMEs_w <- cbind(module,DMEs_w)
    write_excel(DMEs_w,"DMEs",paste0(dirname5,"/",diff_cluster_select,"-vs-otherGroups_table"))
    
    pdf(paste0(dirname5,"/",diff_cluster_select,"-vs-otherGroups_DMEsVolcanoPlot",".pdf"))
    DMEsVolcano <- PlotDMEsVolcano(
      seurat_obj,
      DMEs,
      wgcna_name = wgcna_name,
      plot_labels = TRUE,
      mod_point_size = 4,
      label_size = 4,
      show_cutoff = TRUE
    )
    print(DMEsVolcano)
    dev.off()
    
  }
  
  
  DMEs_all <- FindAllDMEs(
    seurat_obj,
    group.by = 'Group',
    wgcna_name = wgcna_name
  )
  
  p <- PlotDMEsVolcano(
    seurat_obj,
    DMEs_all,
    wgcna_name = wgcna_name,
    plot_labels=FALSE,
    show_cutoff=FALSE
  )
  
  # facet wrap by each cell type
  pdf(paste0(dirname5,"/all_DMEs_VolcanoPlot",".pdf"))
  plot_all <- p + facet_wrap(~group, ncol=3)
  print(plot_all)
  dev.off()
  
  DMEs_all_w <- DMEs_all
  module_all <- rownames(DMEs_all_w)
  DMEs_all_w <- cbind(module_all,DMEs_all_w)
  write_excel(DMEs_all_w,"all_DMEs",paste0(dirname5,"/all_DMEs_table"))
  
  saveRDS(seurat_obj, file=paste0(dirname5,"/hdWGCNA.rds"))
  
  ################################################################################################
  ###--------------------  step6: Module trait correlation  -------------------------###
  dirname6 <- paste0(dirprefix,"/6.ModuleTraitCorrelation")
  dir.create(dirname6)
  
  all_trais_name <- colnames(seurat_obj@meta.data)[-1]
  for (i in 1:length(all_trais_name)){
    trais <- all_trais_name[i]
    if(is.character(class(seurat_obj@meta.data[,`trais`]))){
      seurat_obj@meta.data[,`trais`] <- as.factor(seurat_obj@meta.data[,`trais`])
    }
    
  }
  
  seurat_obj$nCount_RNA <- as.numeric(seurat_obj$nCount_RNA)
  seurat_obj$nFeature_RNA <- as.numeric(seurat_obj$nFeature_RNA)
  seurat_obj$Group <- as.factor(seurat_obj$Group)
  seurat_obj$Sample <- as.factor(seurat_obj$Sample) 
  
  ## add 2023.08.15
  seurat_obj$nCount_metabolite <- seurat_obj$nCount_RNA
  seurat_obj$nFeature_metabolite <- seurat_obj$nFeature_RNA
  
  # list of traits to correlate
  #cur_traits <- c('nCount_metabolite', 'nFeature_metabolite', 'Group', 'Sample')
  # 判断用于相关性分析的trais，如果为all，则计算Module和所有trais列的相关性。否则只计算输入参数中涉及的trais列。
  if (traiscol == "all"){
    cur_traits <-  colnames(df_group)
    trais_01 <- cur_traits[!(cur_traits %in% c("Group", "Sample"))]
    
  }else{
    cur_traits <- strsplit(traiscol, ",")[[1]]
    trais_01 <- cur_traits
  }
  
  ## 将数值型 0和1转换成数值型
  for (each in trais_01){
    seurat_obj@meta.data[,`each`] <- as.numeric(seurat_obj@meta.data[,`each`])
  }
  
  ## ModuleTraitCorrelation() 函数，稍作改写
  
  seurat_obj <- ModuleTraitCorrelation(
    seurat_obj,
    traits = cur_traits,
    group.by='Group',  ##'Cluster'  ##
    features = "average",
    cor_method = cormethod,  # cor_meth: Which method to use for correlation? Valid choices are pearson, spearman, kendall.
    subset_by = NULL,
    subset_groups = NULL,
    wgcna_name = NULL
  )
  
  mt_cor <- GetModuleTraitCorrelation(seurat_obj)
  allnames <- names(mt_cor$cor)
  
  all_cor <- mt_cor$cor
  
  default_plot <- "n"
  if(default_plot == "y"){
    pdf(paste0(dirname6,"/merged_ModuleTraitCor",".pdf"))
    p <- PlotModuleTraitCorrelation(
      seurat_obj,
      label = 'fdr',
      label_symbol = 'stars',
      text_size = 4,
      text_digits = 4,
      text_color = 'white',
      high_color = '#E89189',
      mid_color = '#EEEEEE',
      low_color = '#4777B4',
      plot_max = 0.2,
      combine=TRUE
    )
    print(p)
    dev.off()
    
    sink(paste0(dirname6,"/Cor_pval_fdr.txt"))
    print(mt_cor)
    sink()
    
    saveRDS(mt_cor, file=paste0(dirname6,"/Cor_pval_fdr.rds"))
    
    p <- PlotModuleTraitCorrelation(
      seurat_obj,
      label = 'fdr',
      label_symbol = 'stars', # 'numeric'
      text_size = 4,
      text_digits = 4,
      text_color = 'white',
      high_color = '#E89189',
      mid_color = '#EEEEEE',
      low_color = '#4777B4',
      plot_max = 0.2,
      combine=FALSE
    )
    
    for (plot in 1:length(p)){
      names_x <- names(p)[plot]
      pdf(paste0(dirname6,"/",names_x,"_ModuleTraitCor",".pdf"))
      print(p[[plot]])
      dev.off()
    }
  }
  
  ## 自己绘制heatmap图
  mt_cor_res <- mt_cor$cor$all_Group
  mt_p_res <- mt_cor$pval$all_Group
  mt_fdr_res <- mt_cor$fdr$all_Group
  
  cmt <- mt_cor_res
  pmt2 <- mt_p_res
  
  cmt.out<-cbind(rownames(cmt),cmt)
  #write.table(cmt.out,file=paste0(dirname7,"/","cor.txt"),sep="\t",row.names=F)
  df <-melt(cmt,value.name="cor")
  df$pvalue <-as.vector(pmt2)
  head(df)
  write.table(df,file=paste0(dirname6,"/","cor-p.xls"),sep="\t",row.names=F)
  
  pmt1<-data.frame(pmt2)
  for(i in 1:length(pmt1[,1]))
  {
    for(j in 1:length(pmt1[1,])){
      corx <- round(cmt[i,j],2)
      if (pmt2[i,j]<=0.001){
        strx <- paste0('***','\n',corx)
        pmt1[i,j]<- strx
        #}else if (pmt2[i,j]>0.0001&pmt2[i,j]<=0.001){
        #  pmt1[i,j]<-'***'
      }else if (pmt2[i,j]>0.001&pmt2[i,j]<=0.01){
        strx <- paste0('**','\n',corx)
        pmt1[i,j]<- strx
      }else if (pmt2[i,j]>0.01&pmt2[i,j]<=0.05){
        strx <- paste0('*','\n',corx)
        pmt1[i,j]<- strx
      }else{
        strx <- paste0('ns','\n',corx)
        pmt1[i,j]<- strx
      }
    }}
  
  mycol<-colorRampPalette(c("navy","white", "firebrick3"))(800)
  pic<-pheatmap(cmt,
                scale = "none",
                cluster_row = F, 
                cluster_col = F, 
                border="white",
                display_numbers = pmt1,
                display_color="white",
                number_color = "lightgray",
                fontsize_number = 18,
                cellwidth = 80,
                cellheight = 80,
                angle_col = c('45'),
                color=mycol,
                fontsize= 25,
                fontsize_row = 20,
                fontsize_col = 20,
                main=paste0(cormethod," - Correlation Heatmap"),
                
                legend_breaks = c(-0.4,-0.2,0,0.2,0.4),
                legend_labels = c("-0.4","-0.2","0","0.2","0.4"))
  ggsave(paste0(dirname6,"/ModuleTrais-cor-p-heatmap.pdf"),plot =pic,width = 12,height = 10,units = c("in"))
  
  save.image(file=paste0("","step6.RData"))
  
  ################################################################################################
  ###--------------------  step7 : 富集分析  -------------------------###
  
  suppressPackageStartupMessages(library(meta))
  merged_anno_col <- "mz"
  
  setwd(wd)
  dirname7 <- paste0(dirprefix,"/7.Enrichment")
  dir.create(dirname7)
  
  df_annotation <- read.xlsx(annotation,sheet = 1)
  
  if(prefix_index == "index-"){
    df_annotation[[`merged_anno_col`]] <- paste0(prefix_index,df_annotation[[`merged_anno_col`]])
  }
  
  merged_all <- merge(all_genesinModule, df_annotation, by.x = "gene_name",by.y=`merged_anno_col`, all = FALSE)
  write_excel(merged_all,"allModule_anno",paste0(dirname7,"/allModule_anno"))
  
  df_annotation_sele <- df_annotation[,c("ID","mz","Metabolites","KEGG")]
  merged_sele <- merge(all_genesinModule, df_annotation_sele, by.x = "gene_name",by.y=`merged_anno_col`, all = FALSE)
  
  all_M <- colnames(seurat_obj@misc[[`wgcna_name`]]$hMEs)
  all_M <- all_M[all_M != "grey"]
  
  setwd(dirname7)
  for(m in 1:length(all_M)){
    m_name <- all_M[m]
    m_df <- merged_sele[merged_sele$module==m_name,]
    m_df <- m_df[,c("gene_name","Metabolites","KEGG")]
    
    
    m_df <- m_df %>%
      mutate(KEGG = gsub(";\\r\\n", ";", KEGG),
             Metabolites = gsub(";\\r\\n", ";", Metabolites))
    m_df_c <- m_df[,c("KEGG")]
    m_df_c <- unlist(str_extract_all(m_df_c, "C[^;]+"))
    m_df_c <- na.omit(m_df_c)
    m_df_c <- unique(m_df_c)
    table_m_df_c <- data.frame(KEGG = m_df_c)
    write_excel(table_m_df_c,paste0(m_name,"_KEGGlist"),paste0(m_name,"_KEGGlist"))
    
    meta::keggrich(name=paste0(m_name,"_KEGGlist.xlsx"),species = species,kegg="KEGG",keggmap=F)
    
  }
  setwd(wd)
  save.image(file=paste0("","step7.RData"))
  
  
  base::print("分析参数存储于 parameter.log")
  sink("parameter.log",append = FALSE)
  base::print(paste0("script_path: ",script_path))
  base::print(paste0("input: ",input))
  base::print(paste0("group: ",group))
  base::print(paste0("annotation: ",annotation))
  base::print(paste0("orders: ",orders))
  base::print(paste0("species: ",species))
  base::print(paste0("hub_n_genes: ",hub_n_genes))
  base::print(paste0("minModuleSize: ",minModuleSize))
  base::print(paste0("deepSplit: ",deepSplit))
  base::print(paste0("networkType: ",networkType))
  base::print(paste0("wgcna_name: ",wgcna_name))
  base::print(paste0("cormethod: ",cormethod))
  base::print(paste0("traiscol: ",traiscol))
  base::print(paste0("kvalue: ",kvalue))
  base::print(paste0("min_cells: ",min_cells))
  base::print(paste0("max_shared: ",max_shared))
  base::print(paste0("dirprefix: ",dirprefix))
  sink()
  
}
