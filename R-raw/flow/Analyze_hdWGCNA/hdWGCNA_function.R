#' @export
mycol <- c("#386CB0","#DB647AB9","#8FBC8F","#7570B3","#c6c386","#925E9F","#F78F27",
           "#c08b74","#4E84C4","#949494","#7FC97F","#BEAED4","#FFC1C1","#FDC086",
           "#F0027F","#BF5B17","#666666","#1B9E77","#7570B3","#66A61E","#E6AB02",
           "#A6761D","#A6CEE3","#B2DF8A","#E89189","#E31A1C","#FF7F00","#6A3D9A",
           "#8DA0CB","#4DAF4A","#984EA3","#999999","#66C2A5","#A49A2E","#FC8D62",
           "#A6D854","#FFD92F","#BEBADA","#DA7EB4","#80B1D3","#FDB462","#52854C",
           "#FB8072","#BC80BD","#B3B3B3","#33A02C","#FC8D62","#B3DE69","#4038b0",
           "#ee7576","#E78AC3","#ff0000","#A65628","#F781BF","#D95F02","#E7298A",
           "#B15928","#1F78B4","#FDBF6F","#CAB2D6","#B15928","#FB9A99","#B3CDE3",
           '#0173b2','#de8f05','#029e73','#d55e00','#cc78bc','#ca9161','#fbafe4',
           '#293352','#ece133','#56b4e9',"#00AFBB","#E7B800","#FC4E07","#FFDB6D",
           "#C4961A","#F4EDCA","#D16103","#C3D7A4","#e94749","#d80172","#A85293")

#' @export
mycol_module <- c("#DA7EB4","#80B1D3","#FDB462","#52854C",
                  "#FB8072","#BC80BD","#33A02C","#FC8D62","#B3DE69","#4038b0",
                  "#ee7576","#E78AC3","#ff0000","#A65628","#F781BF","#D95F02","#E7298A",
                  "#B15928","#1F78B4","#FDBF6F","#CAB2D6","#B15928","#FB9A99","#B3CDE3",
                  '#0173b2','#de8f05','#029e73','#d55e00','#cc78bc','#ca9161','#fbafe4',
                  '#293352','#ece133','#56b4e9',"#00AFBB","#E7B800","#FC4E07","#FFDB6D",
                  "#C4961A","#F4EDCA","#D16103","#C3D7A4","#e94749","#d80172","#A85293")


#================================================================================
# Function
#================================================================================
# 将dataframe转换为excel表格
#' @export
write_excel <- function(
    df,
    sheet_name = "test",
    output_name = "test"
){
  if(file.exists(paste0(output_name,".xlsx"))){
    file.remove(paste0(output_name,".xlsx"))
  }
  wb <- createWorkbook() # 创建 Excel 工作簿
  addWorksheet(wb, sheet_name) # 添加一个名为 "my_sheet" 的工作表
  writeData(wb, sheet_name, df) # 在 "my_sheet" 中写入数据
  saveWorkbook(wb, paste0(output_name,".xlsx")) # 将工作簿保存为 Excel 文件
}

## plot dimension
#' @export
plot_dimension <- function(
    seurat_obj,
    func,# RunPCA/RunTSNE/RunUMAP/RunHarmony
    prefix, # 输出文件的前缀
    reduction_type, # 降维聚类的类型 pca/tsne/umap/harmony
    reduction_title,
    groupby1, # 分组因素1
    groupby2, # 分组因素2
    dirname = dirname1,
    mycol,
    width = 8,
    height = 6
){
  pdf(paste0(dirname1,"/",prefix=prefix,"_",reduction_title,".pdf"),width=width,height=height)
  p_pca_1 <- func(seurat_obj, 
                  group.by=groupby1, 
                  cols = mycol[1:length(unique(seurat_obj$Group))],
                  label=TRUE,
                  reduction=reduction_type) +
    ggtitle(paste0(reduction_title,": ",groupby1)) 
  p_pca_2 <- func(seurat_obj, 
                  group.by=groupby2,
                  cols = mycol[1:length(unique(seurat_obj$Sample))],
                  label=TRUE,
                  reduction=reduction_type) +
    ggtitle(paste0(reduction_title,": ",groupby2)) 
  #NoLegend()
  p_pca <- p_pca_1 | p_pca_2
  print(p_pca)
  dev.off()
}

## raw colour 2 new colour, raw colour with names.
#' @export
colour_raw2new <- function(
    raw_col, # 原始colour, character类型,但是有names
    mycol_module   # 给定的自定义颜色列表, character类型
){
  unique_colors <- unique(raw_col)
  grey_pos <- which(unique_colors == "grey")
  new_col_unique <- mycol_module[1:(length(unique_colors)-1)]
  new_col_unique <- insert(new_col_unique,grey_pos,"grey")
  new_col <- cbind(unique_colors,new_col_unique)
  colnames(new_col) <- c("raw","new")
  return(new_col) ### 生成raw2new的两列矩阵
}

# 从raw clolour转成new colour, 输入和输出都是character类型
#' @export
change_colour <- function(
    raw_col, # 原始colour, character类型
    raw2new_col # raw2new的两列矩阵
){
  new_col <- ifelse(raw_col %in% raw2new_col[,"raw"], raw2new_col[,"new"][match(raw_col, raw2new_col[,"raw"])], raw_col)
  return(new_col)
}


#####  ----------------------------- 绘图 ------------------------------------------
## 
#' @export
featurePlot_plot <- function(seurat_obj,feat_type,ucell_F,oder_type,reduction_type,restrict_range){
  plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features= feat_type, # What to plot? Can select hMEs, MEs, scores, or average
    ucell = ucell_F,  
    order=oder_type, # choose from: TRUE/FALSE/shuffle. order so the points with highest hMEs are on top. shuffle:order so cells are shuffled
    reduction = reduction_type,
    restrict_range = restrict_range,#TRUE,
    point_size = 0.5,
    alpha = 1,
    label_legend = FALSE,
    raster = FALSE,
    raster_dpi = 500,
    raster_scale = 1,
    plot_ratio = 1,
    title = TRUE
  )
  return(plot_list)
}

# plot module correlagram 绘制模块相关性图函数
#' @export
ModuleCorrelogram_plot <- function(seurat_obj,feat_type){
  p <- ModuleCorrelogram(
    seurat_obj,
    MEs2 = NULL,
    features = feat_type, # What to plot? Can select hMEs, MEs, scores, or average
    order = "original",
    method = "ellipse",
    exclude_grey = TRUE,
    type = "upper",
    tl.col = "black",
    tl.srt = 45,
    sig.level = c(1e-04, 0.001, 0.01, 0.05),
    pch.cex = 0.7,
    col = colorRampPalette(c('SkyBlue3',  'white', 'pink'))(200),
    ncolors = 200,
    wgcna_name = NULL,
    wgcna_name2 = NULL,
  )
  return(p)
  
}

#' @export
all_plot_1 <- function(seurat_obj,mode_x,f1,t1,reduction_type,dirnamex,output1,output2,output1_sheetname,restrict_range){
  # # 绘制 hMEs, MEs, scores, or average 的降维投影函数
  plot_umap <- "n"
  if(plot_umap == "y"){
    plot_list <- featurePlot_plot(seurat_obj,mode_x,FALSE,TRUE,reduction_type,restrict_range)
    pdf(paste0(dirname3,output1,reduction_type,".pdf"),width=12)
    p <- wrap_plots(plot_list, ncol=6)
    print(p)
    dev.off()
  }
  # module correlation模块相关性
  pdf(paste0(dirname3,"/",output2,".pdf"))
  p <- ModuleCorrelogram_plot(seurat_obj,mode_x)
  print(p)
  dev.off()
  p_cor <- p$corrPos
  write_excel(p_cor,"corr-p",paste0(dirname3,"/",output2))
  
  
  # harmonized module eigengenes: 使用harmony按照样本进行批次矫正
  if(mode_x == "average"){
    xx_module <- GetAvgModuleExpr(seurat_obj)
  }else if(mode_x == "MEs"){
    xx_module <- GetMEs(seurat_obj, harmonized=FALSE)
  }else if(mode_x == "scores"){
    xx_module <- GetModuleScores(seurat_obj)
  }
  index <- rownames(xx_module)
  xx_module_w <- cbind(index,xx_module)
  write_excel(xx_module_w,output1_sheetname,paste0(dirname3,output1,"table"))
  
}

### plot lines
#' @export
plotline <- function(df_data,legend_TF,xlab_TF,title,mycolour){
  colnames(df_data) <- c("Group","Feature","Value")
  min_v <- min(df_data$Value)
  if(min_v<= 0) {
    plus <- abs(min_v) + 1
  }else{
    plus <- 0
  }
  
  p <- ggplot(data=df_data,aes(x=as.numeric(Group),y=log(Value + plus,10),color=Feature))+
    geom_line(linewidth=1.5)+
    geom_point(size=1.5,color="darkgray")+
    scale_color_manual(values = mycolour) +
    theme_classic()+
    xlab("Group")+
    ylab("log10 (Value) ")+
    scale_x_continuous(
      #breaks = c(1,2,3,4,5),
      breaks = c(1:length_period),
      label = raw_level) +
    guides(color=guide_legend(title = "Feature"))+
    theme(axis.title.x=element_text(vjust=0, size=20,face = "bold"))+  # 坐标轴副标题字体大小
    theme(axis.title.y=element_text(vjust=0, size=20,face = "bold"))+  # 坐标轴副标题字体大小
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+  ## 坐标轴线粗细
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+  ## 坐标轴线粗细
    theme(axis.text.x=element_text(vjust=0.5,size=14,angle=45))+  ## 坐标轴标签大小
    theme(axis.text.y=element_text(vjust=1,size=20))+  ## 坐标轴标签大小
    theme(legend.title=element_text(size=20),legend.text=element_text(size=20))   ##前面改lengend title大小;后面改legend text大小
  if(legend_TF == "without"){
    p <- p + theme(legend.position = "none")
  }
  if(xlab_TF == "without"){
    p <- p+ labs(title = title,x = "",y = "") + theme(plot.title = element_text(hjust = 0.5))
  }
  return(p)
  
}

#' @export
plot_smooth <- function(df_data,colour,legend_TF,xlab_TF,title,mycolour){
  colnames(df_data) <- c("Group","Feature","Value")
  min_v <- min(df_data$Value)
  if(min_v<= 0) {
    plus <- abs(min_v) + 1
  }else{
    plus <- 0
  }
  p <- ggplot(data=df_data,aes(x=as.numeric(Group),y=log(Value+ plus,10)))+
    geom_point(size=2.5,color="darkgray")+
    geom_smooth(method = "loess", show.legend = FALSE, fill = colour,color=colour)+
    scale_color_manual(values = mycolour ) +
    theme_classic()+
    xlab("Group")+
    ylab("log10 (Value) ")+
    scale_x_continuous(
      breaks =  c(1:length_period),
      label = raw_level) +
    theme(axis.title.x=element_text(vjust=0, size=20,face = "bold"))+  # 坐标轴副标题字体大小
    theme(axis.title.y=element_text(vjust=0, size=20,face = "bold"))+  # 坐标轴副标题字体大小
    theme(axis.line.x=element_line(linetype=1,color="black",size=1.2))+  ## 坐标轴线粗细
    theme(axis.line.y=element_line(linetype=1,color="black",size=1.2))+  ## 坐标轴线粗细
    theme(axis.text.x=element_text(vjust=0.5,size=14,angle=45))+  ## 坐标轴标签大小
    theme(axis.text.y=element_text(vjust=1,size=20))+  ## 坐标轴标签大小
    theme(legend.title=element_text(size=20),legend.text=element_text(size=20))   ##前面改lengend title大小;后面改legend text大小
  if(legend_TF == "without"){
    p <- p + theme(legend.position = "none")
  }
  if(xlab_TF == "without"){
    p <- p+ labs(title = title,x = "",y = "") + theme(plot.title = element_text(hjust = 0.5))
  }
  return(p)
}

#' @export
plot_VLNandDot <- function(seurat_obj,wgcna_name,hub_df,mods,modx,modname,output1,output2,output3,dirname3,mycol){
  
  # plot with Seurat's DotPlot function
  p <- DotPlot(seurat_obj, 
               features= paste0(modname,mods),
               group.by = 'Group',
               cols = c("lightgrey", "blue"),
               col.min = -2.5,
               col.max = 2.5,
               dot.min = 0,
               dot.scale = 6,
               idents = NULL,
               split.by = NULL,
               cluster.idents = FALSE,
               scale = TRUE,
               scale.by = "radius",
               scale.min = NA,
               scale.max = NA
  )
  
  # flip the x/y axes, rotate the axis labels, and change color scheme:
  p <- p +
    coord_flip() +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue')
  
  # plot output
  pdf(paste0(dirname3,"/",output1))
  print(p)
  dev.off()
  
  module_type <- unique(hub_df$module)
  
  all_M <- colnames(seurat_obj@misc[[`wgcna_name`]]$hMEs)
  all_M <- all_M[all_M != "grey"]
  
  ln <- length(unique(seurat_obj@meta.data$Group))
  mycol2 <-  alpha(mycol,0.6)[1:ln]
  
  ## Vlnplot的是基于hMEs，而hMEs是存储在 seurat_obj$  中的那些Module特征
  vlnP <- NULL
  for (M in 1:length(all_M)){
    # Plot INH-M4 hME using Seurat VlnPlot function
    M_name <- all_M[M]
    p <- VlnPlot(
      seurat_obj,
      features = paste0(modname,M_name),
      group.by = 'Group',
      cols = mycol2,  #NULL,
      pt.size = 0, # don't show actual data points
      idents = NULL,
      sort = FALSE,
      assay = NULL,
      split.by = NULL,
      adjust = 1,
      y.max = NULL,
      same.y.lims = FALSE,
      log = FALSE,
      ncol = NULL,
      slot = "data",
      split.plot = FALSE,
      stack = FALSE,
      combine = TRUE,
      fill.by = "feature",
      flip = FALSE,
      add.noise = TRUE,
      raster = NULL
    )
    
    # add box-and-whisker plots on top:
    p <- p + geom_boxplot(width=.25, fill='white')
    
    # change axis labels and remove legend:
    p1 <- p + xlab('') + ylab(modx) + NoLegend()
    p2 <- p1 + labs(title = M_name,x = "",y = "") 
    
    if(M == 1){
      vlnP[[1]] <- p2
    }else{
      vlnP[[M]] <- p2
    }
    
    # plot output
    pdf(paste0(dirname3,"/individual/",M_name,output2))
    print(p1)
    dev.off()
  }
  
  ncol_n <- 3  ## 每行排几个
  
  if(length(module_type)%%ncol_n == 0 ){
    height_n <- length(module_type)/ncol_n
  }else{
    height_n <- (length(module_type)+(ncol_n-length(module_type)%%ncol_n))/ncol_n
  }
  pdf(paste0(dirname3,"/",output3),width=(ncol_n*3),height=(3.5*(height_n)))
  grid.arrange(grobs =vlnP, ncol = ncol_n)
  dev.off()
  
}

# plot Module Network
#' @export
ModuleNetworkPlot <- function(
    seurat_obj,
    mods="all",
    outdir = dirname4,
    plot_size = c(6,6),
    wgcna_name=NULL,
    label_center = FALSE, # only label the genes in the middle?
    edge.alpha=0.25,
    vertex.label.cex=1,
    vertex.size=6
){
  # get data from active assay if wgcna_name is not given
  if(is.null(wgcna_name)){wgcna_name <- seurat_obj@misc$active_wgcna}
  
  # get modules, MEs:
  MEs <- GetMEs(seurat_obj, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)
  
  # using all modules?
  if(mods == 'all'){
    mods <- levels(modules$module)
    mods <- mods[mods != 'grey']
  }
  
  # check if we have eigengene-based connectivities:
  if(!all(paste0('kME_', as.character(mods)) %in% colnames(modules))){
    stop('Eigengene-based connectivity (kME) not found. Did you run ModuleEigengenes and ModuleConnectivity?')
  }
  
  # create output folder
  if(!dir.exists(outdir)){dir.create(outdir)}
  
  # tell the user that the output is going to the output dir
  cat(paste0("Writing output files to ", outdir))
  
  # get TOM
  TOM <- GetTOM(seurat_obj, wgcna_name)
  
  # get hub genes:
  hub_n_genes <- 25
  n_hubs <- hub_n_genes #25
  hub_list <- lapply(mods, function(cur_mod){
    cur <- subset(modules, module == cur_mod)
    cur <- cur[,c('gene_name', paste0('kME_', cur_mod))] %>%
      top_n(n_hubs)
    colnames(cur)[2] <- 'var'
    cur %>% arrange(desc(var)) %>% .$gene_name
  })
  names(hub_list) <- mods
  
  # loop over modules
  for(cur_mod in mods){
    xx <- which(mods == cur_mod)
    print(cur_mod)
    cur_color <- modules %>% subset(module == cur_mod) %>% .$color %>% unique
    
    # number of genes, connections
    # might make a setting to change this later but I will also have to change
    # how the graph layout works
    #n_genes = hub_n_genes #25;  ### !!! modified 2023.06.19
    n_conns = 500;
    
    # name of column with current kME info
    cur_kME <- paste0('kME_', cur_mod)
    
    cur_genes <- hub_list[[cur_mod]]
    
    # Identify the columns in the TOM that correspond to these hub genes
    matchind <- match(cur_genes, colnames(TOM))
    reducedTOM = TOM[matchind,matchind]
    orderind <- order(reducedTOM,decreasing=TRUE)
    
    if(dim(as.matrix(reducedTOM))[1] >= 25){
      n_genes <- hub_n_genes
    }else{
      n_genes <- dim(as.matrix(reducedTOM))[1]
    }
    
    # only  keep top connections
    connections2keep <- orderind[1:n_conns];
    reducedTOM <- matrix(0,nrow(reducedTOM),ncol(reducedTOM));
    reducedTOM[connections2keep] <- 1;
    
    # print('here')
    # print(dim(reducedTOM))
    # print(n_genes)
    
    # only label the top 10 genes?
    # if(label_center){cur_genes[11:25] <- ''}
    if(label_center){cur_genes[11:hub_n_genes] <- ''}  ### !!! modified 2023.06.19
    
    # top 10 as center
    gA <- graph.adjacency(as.matrix(reducedTOM[1:10,1:10]),mode="undirected",weighted=TRUE,diag=FALSE)
    gB <- graph.adjacency(as.matrix(reducedTOM[11:n_genes,11:n_genes]),mode="undirected",weighted=TRUE,diag=FALSE)
    layoutCircle <- rbind(layout.circle(gA)/2,layout.circle(gB))
    
    g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
    
    pdf(paste0(outdir, '/', cur_mod,'_top',n_genes,'.pdf'), width=plot_size[1], height=plot_size[2], useDingbats=FALSE);
    plot(g1,
         edge.color=adjustcolor(cur_color, alpha.f=0.25),
         edge.alpha=edge.alpha,
         vertex.color=cur_color,
         vertex.label=as.character(cur_genes),
         vertex.label.dist=1.1,
         vertex.label.degree=-pi/4,
         vertex.label.color="black",
         vertex.label.family='Helvetica',
         vertex.label.font = 3,
         vertex.label.cex=vertex.label.cex,
         vertex.frame.color='black',
         layout= jitter(layoutCircle),
         vertex.size=vertex.size,
         main=paste(cur_mod)
    )
    dev.off();
    
  }
}

#' @export
ModuleTraitCorrelation <- function (seurat_obj, traits, group.by,features,
                                    cor_method, subset_by = NULL, subset_groups = NULL, 
                                    wgcna_name = NULL)
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  if (features == "hMEs") {
    MEs <- GetMEs(seurat_obj, TRUE, wgcna_name)
  }
  else if (features == "MEs") {
    MEs <- GetMEs(seurat_obj, FALSE, wgcna_name)
  }
  else if (features == "scores") {
    MEs <- GetModuleScores(seurat_obj, wgcna_name)
  }
  else if(features == 'average'){
    MEs <- GetAvgModuleExpr(seurat_obj, wgcna_name)
    restrict_range <- FALSE
  }else {
    stop("Invalid feature selection. Valid choices: hMEs, MEs, scores, average")
  }
  if (!is.null(subset_by)) {
    print("subsetting")
    seurat_full <- seurat_obj
    MEs <- MEs[seurat_obj@meta.data[[subset_by]] %in% subset_groups, 
    ]
    seurat_obj <- seurat_obj[, seurat_obj@meta.data[[subset_by]] %in% 
                               subset_groups]
  }
  if (sum(traits %in% colnames(seurat_obj@meta.data)) != length(traits)) {
    stop(paste("Some of the provided traits were not found in the Seurat obj:", 
               paste(traits[!(traits %in% colnames(seurat_obj@meta.data))], 
                     collapse = ", ")))
  }
  if (is.null(group.by)) {
    group.by <- "temp_ident"
    seurat_obj$temp_ident <- Idents(seurat_obj)
  }
  valid_types <- c("numeric", "factor", "integer")
  data_types <- sapply(traits, function(x) {
    class(seurat_obj@meta.data[, x])
  })
  if (!all(data_types %in% valid_types)) {
    incorrect <- traits[!(data_types %in% valid_types)]
    stop(paste0("Invalid data types for ", paste(incorrect, 
                                                 collapse = ", "), ". Accepted data types are numeric, factor, integer."))
  }
  if (any(data_types == "factor")) {
    factor_traits <- traits[data_types == "factor"]
    for (tr in factor_traits) {
      warning(paste0("Trait ", tr, " is a factor with levels ", 
                     paste0(levels(seurat_obj@meta.data[, tr]), collapse = ", "), 
                     ". Levels will be converted to numeric IN THIS ORDER for the correlation, is this the expected order?"))
    }
  }
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != "grey"]
  trait_df <- seurat_obj@meta.data[, traits]
  if (length(traits == 1)) {
    trait_df <- data.frame(x = trait_df)
    colnames(trait_df) <- traits
  }
  if (any(data_types == "factor")) {
    factor_traits <- traits[data_types == "factor"]
    for (tr in factor_traits) {
      trait_df[, tr] <- as.numeric(trait_df[, tr])
    }
  }
  cor_list <- list()
  pval_list <- list()
  fdr_list <- list()
  temp <- Hmisc::rcorr(as.matrix(trait_df), as.matrix(MEs), 
                       type = cor_method)
  cur_cor <- temp$r[traits, mods]
  cur_p <- temp$P[traits, mods]
  p_df <- cur_p %>% reshape2::melt()
  if (length(traits) == 1) {
    tmp <- rep(mods, length(traits))
    tmp <- factor(tmp, levels = mods)
    tmp <- tmp[order(tmp)]
    p_df$Var1 <- traits
    p_df$Var2 <- tmp
    rownames(p_df) <- 1:nrow(p_df)
    p_df <- dplyr::select(p_df, c(Var1, Var2, value))
  }
  p_df <- p_df %>% dplyr::mutate(fdr = p.adjust(value, method = "fdr")) %>% 
    dplyr::select(c(Var1, Var2, fdr))
  cur_fdr <- reshape2::dcast(p_df, Var1 ~ Var2, value.var = "fdr")
  rownames(cur_fdr) <- cur_fdr$Var1
  cur_fdr <- cur_fdr[, -1]
  cor_list[["all_Group"]] <- cur_cor
  pval_list[["all_Group"]] <- cur_p
  fdr_list[["all_Group"]] <- cur_fdr
  trait_df <- cbind(trait_df, seurat_obj@meta.data[, group.by])
  colnames(trait_df)[ncol(trait_df)] <- "group"
  MEs <- cbind(as.data.frame(MEs), seurat_obj@meta.data[, group.by])
  colnames(MEs)[ncol(MEs)] <- "group"
  if (class(seurat_obj@meta.data[, group.by]) == "factor") {
    group_names <- levels(seurat_obj@meta.data[, group.by])
  }else {
    group_names <- levels(as.factor(seurat_obj@meta.data[, 
                                                         group.by]))
  }
  trait_list <- dplyr::group_split(trait_df, group, .keep = FALSE)
  ME_list <- dplyr::group_split(MEs, group, .keep = FALSE)
  names(trait_list) <- group_names
  names(ME_list) <- group_names
  for (i in names(trait_list)) {
    temp <- Hmisc::rcorr(as.matrix(trait_list[[i]]), as.matrix(ME_list[[i]]))
    cur_cor <- temp$r[traits, mods]
    cur_p <- temp$P[traits, mods]
    p_df <- cur_p %>% reshape2::melt()
    if (length(traits) == 1) {
      tmp <- rep(mods, length(traits))
      tmp <- factor(tmp, levels = mods)
      tmp <- tmp[order(tmp)]
      p_df$Var1 <- traits
      p_df$Var2 <- tmp
      rownames(p_df) <- 1:nrow(p_df)
      p_df <- dplyr::select(p_df, c(Var1, Var2, value))
    }
    p_df <- p_df %>% dplyr::mutate(fdr = p.adjust(value, 
                                                  method = "fdr")) %>% dplyr::select(c(Var1, Var2, 
                                                                                       fdr))
    cur_fdr <- reshape2::dcast(p_df, Var1 ~ Var2, value.var = "fdr")
    rownames(cur_fdr) <- cur_fdr$Var1
    cur_fdr <- cur_fdr[, -1]
    cor_list[[i]] <- cur_cor
    pval_list[[i]] <- cur_p
    fdr_list[[i]] <- as.matrix(cur_fdr)
  }
  mt_cor <- list(cor = cor_list, pval = pval_list, fdr = fdr_list)
  if (!is.null(subset_by)) {
    seurat_full <- SetModuleTraitCorrelation(seurat_full, 
                                             mt_cor, wgcna_name)
    seurat_obj <- seurat_full
  }
  else {
    seurat_obj <- SetModuleTraitCorrelation(seurat_obj, mt_cor, 
                                            wgcna_name)
  }
  seurat_obj
}

