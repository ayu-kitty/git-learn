#!/opt/conda/bin/Rscript
#' 激酶-底物相关性网络图
#'
#' @param savepath 保存路径，默认Network
#' @param or 物种缩写，hsa、mmu、rno，默认hsa
#' @param phosfile 输入的磷酸化位点表，默认表达矩阵-phos.xlsx
#' @param profile 输入的蛋白表，默认表达矩阵-pro.xlsx
#' @param corr 激酶底物的相关性卡值，默认0。8
#' @param imagetype 输出图片类型
#' @param height 图片高度
#' @param width 图片宽度
#' @param dpi 图片分辨率
#' @param fontfamily 字体
#' @param ... 
#' @export
KSEA_corr<- function(savepath = "./Network/",species='hsa' ,phosfile = "表达矩阵-ptm.xlsx",profile = "表达矩阵-pro.xlsx",inputpath='rawdata/',
                     corr=0.8,imagetype = c("pdf", "png") , height = 12, width = 14, dpi = 300, fontfamily = "sans", ...){
  pacman::p_load(readxl,stringr,openxlsx,readxl,tidyr,ggpubr,igraph,ggraph,pheatmap,tidygraph,RColorBrewer,ggplot2,colorRamps,R2BayesX,optparse)

  createdir(savepath)
  da_sa_1<- read_excel(paste0(inputpath,"Sample_Group.xlsx"),sheet = 1)
  da_sa_2<- read_excel(paste0(inputpath,"Sample_Group.xlsx"),sheet = 2)
  
  if(species=="hsa"){
    org<-"human"
  }else if(species =="mmu"){
    org<-"mouse"
  }else org<-"rat"
  Kinase <- read_excel("/data/hstore3/public/propip/Kinase_Substrate_Dataset-2023.xlsx")
  kin1<- Kinase[which(Kinase$SUB_ORGANISM==org&Kinase$KIN_ORGANISM==org),1:10]#底物物种筛选后的激酶表格kin1
  kin_Acc<-as.data.frame(unique(kin1[,1:3]))
  name_kin1<- c("KIN_ACC_ID","KINASE","GENE")
  #####读取差异蛋白表，匹配激酶
  pro<- read_excel(paste0(inputpath,profile))%>%as.data.frame()
  pro$KIN_ACC_ID<- pro$Accession
  pro_1<- left_join(x=pro,y=kin_Acc,by = c("KIN_ACC_ID"),multiple = "all")
  pro_result<- pro_1[,-(which("KIN_ACC_ID"==colnames(pro_1)))]
  write.xlsx(pro_result, file = paste0(savepath,"蛋白表达矩阵_ksea.xlsx"))
  
  #####读取差异位点表，匹配激酶
  phos<- read_excel(paste0(inputpath,phosfile))%>%as.data.frame()
  phos <- separate(phos,ID,into = c("SUB_ACC_ID","SUB_MOD_RSD"),sep = ":",remove = F)#根据分列添加基因名和蛋白名列
  phos_1<- left_join(x=phos,y=kin1,by = c("SUB_ACC_ID","SUB_MOD_RSD"),multiple = "all")
  phos_result<- phos_1[,-(2:3)]
  write.xlsx(phos_result, file = paste0(savepath,"位点表达矩阵_ksea.xlsx")) 
  #####可信蛋白对应的激酶
  print(~"..........提取蛋白表对应的激酶..........")
  pro_kin1<- as.data.frame(pro_result)
  pro_kin<- pro_kin1[which(!is.na(pro_kin1$KINASE)),which(colnames(pro_kin1)%in%c('Accession',"Gene Name",da_sa_1$Sample,"KINASE"))]
  
  ######可信位点对应的激酶
  print(~"..........提取位点表对应的激酶..........")
  phos1<- phos_result
  phos_kin<- phos1[which(!is.na(phos1$KINASE)),which(colnames(phos1)%in%c('ID',"Gene Name",da_sa_1$Sample,"KINASE"))]
  phos_kin$`Gene Name`<- sapply(1:nrow(phos_kin),function(...){
    paste0(phos_kin$`Gene Name`[...],":",(str_split(phos_kin$ID[...],":")%>%unlist())[2])
  })
  colnames(phos_kin)[1]<- "Accession"            
  ####建立变量分类表
  ##根据激酶合并蛋白激酶和磷酸化蛋白表
  pro_phos<- merge(pro_kin,phos_kin,by="KINASE")
  colnames(pro_phos)<- str_replace(colnames(pro_phos),".x","_pro")%>%str_replace(.,".y","_phos")
  ####建立蛋白名变量分类表data_type
  type_kin<-as.data.frame(unique(c(pro_phos$Accession_pro,pro_phos$`Gene Name_pro`)))
  type_kin$type<- c('kinase')
  colnames(type_kin)<- c('node',"type")
  type_sub<-as.data.frame(unique(c(pro_phos$Accession_phos,pro_phos$`Gene Name_phos`)))
  type_sub$type<- c('Substrate')
  colnames(type_sub)<- c('node',"type")
  data_type<- rbind(type_kin,type_sub)##变量分类表合并
  #计算蛋白组和磷酸化蛋白的样本相关性
  na_pro<- paste0(da_sa_1$Sample,"_pro")
  na_phos<- paste0(da_sa_1$Sample,"_phos")
  
  pro_phos1<- pro_phos
  gspro<- pro_phos1[,na_pro]
  gsage<- pro_phos1[,na_phos]
  pro_phos1[,"corr"]<-sapply(1:nrow(gspro), function(x){cor(as.numeric(gspro[x,]),as.numeric(gsage[x,]))})
  pro_phos1[,"corr_pvalue"]<-sapply(1:nrow(gspro), function(x){cor.test(as.numeric(gspro[x,]),as.numeric(gsage[x,]))$p.value})
  pro_phos1[,"p-adjust"]=p.adjust(pro_phos1$corr_pvalue, "BH")
  ##导出有激酶、激酶编号、磷酸化蛋白、位点、相关性系数的表格
  dd_out1<- separate(pro_phos1,Accession_phos,into = c("Accession_phos","site"),sep = ":",remove = T)
  dd_out<- separate(dd_out1,`Gene Name_phos`,into = c("Gene Name_phos",NA),sep = ":",remove = T)%>%arrange(.,corr_pvalue)
  write.xlsx(dd_out[,-1], file = paste0(savepath,"激酶底物相关性系数.xlsx"))
  #绘制相关性散点图
  if(nrow(dd_out)>20){
    ks<-head(dd_out,20)
  }else ks<-dd_out
  if (!dir.exists(paste0(savepath,"corrpoint"))){
    dir.create(paste0(savepath,"corrpoint"))
  }
  for(i in 1:nrow(ks)){
    proc<-data.frame(pro=as.numeric(ks[i,na_pro]),phos=as.numeric(ks[i,na_phos]))
    proc[,"lab"]<-paste0(ks$`Gene Name_pro`[i],"-",ks$`Gene Name_phos`[i],"_",ks$site[i])
    ggplot(proc, aes(x = pro, y =phos ))+
      geom_point(size=4)+
      geom_smooth(method = "lm", fullrange = TRUE,se=FALSE,size=3) +
      facet_wrap(~lab) +
      theme(text=element_text(family="ArialMT"))+
      theme(panel.grid.major=element_line(color="gray88",size=2),
            panel.grid.minor=element_line(color="gray88",size=2),
            panel.background = element_rect(fill = "gray99"),
            panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
      stat_cor(data = proc[,1:2],method = "pearson",size = 8, label.y = max(proc$phos)+abs(0.1*max(proc$phos))) +
      labs(y = "Substrate", x = "Kinase") +
      theme(plot.title = element_text(hjust = 0.5, size=20))+
      theme(axis.text = element_text(size = 20,colour = "black"),
            axis.title = element_text(size = 28) ,strip.text = element_text(size = 40),
            strip.background=element_rect(fill="gray88",color="black", size=2, linetype="solid"),
            strip.placement = "outside",
            strip.text.x = element_text(margin=unit(rep(1,4),"lines")))+
      theme(plot.margin = unit(rep(3,4),"lines"))->p
    ggsave(paste0(savepath,"corrpoint/",proc[1,3],".pdf"), height=10, width=10, plot=p)
    ggsave(paste0(savepath,"corrpoint/",proc[1,3],".png"), type="cairo-png", height=10, width=10, plot=p)
  }
  ###整理绘图数据
  file_corr<- pro_phos1[which(as.numeric(pro_phos1$corr)>corr|as.numeric(pro_phos1$corr)<(-corr)),]
  if (nrow(file_corr)<2) {
    print(paste0("相关性系数",corr,"过高，无法筛选激酶，请重新输入相相关性值"))}else{
      print(paste0("相关性系数：",corr))
    }
  dd_pro_phos<- as.data.frame(cbind(file_corr$Accession_pro,file_corr$`Gene Name_pro`,file_corr$Accession_phos,file_corr$`Gene Name_phos`,file_corr$corr))
  savetxt(data = corr,filename = paste0(inputpath,"Pearson_corr.txt"),append = F,row.names = T,)
  
  ######开始绘图
  ######1.节点以蛋白名显示的连接度表
  c(as.character(dd_pro_phos[,1]),as.character(dd_pro_phos[,3]))%>%as_tibble%>%
    group_by(value)%>%
    dplyr::summarize(n=n())->vertices
  colnames(vertices)<- c("node","n")
  ##添加变量分类属性
  vertices%>%select(-n)%>%left_join(data_type,by="node")->vertices
  ##构建graph结构数据
  AC_data<- dd_pro_phos[order(dd_pro_phos$V5),]%>% select(1,3,5)
  colnames(AC_data)<- c("kinase","Substrate","corr")
  AC<- graph_from_data_frame(AC_data,vertices,directed = T)
  print(paste0("节点数目:",vcount(AC)))#节点数目
  print(paste0("链接数:",ecount(AC)))#链接数
  ######2.节点以基因名显示的连接度表
  c(as.character(dd_pro_phos[,2]),as.character(dd_pro_phos[,4]))%>%as_tibble%>%
    group_by(value)%>%
    dplyr::summarize(n=n())->vertices2
  colnames(vertices2)<- c("node","n")
  ##添加变量分类属性
  vertices2%>%select(-n)%>%left_join(data_type,by="node")->vertices2
  ##构建graph结构数据
  GN_data<- dd_pro_phos[order(dd_pro_phos$V5),]%>% select(2,4,5)
  colnames(GN_data)<- c("kinase","Substrate","corr")
  GN<- graph_from_data_frame(GN_data,vertices2,directed = T)
  ##导出表
  net_query<- igraph::as_data_frame(AC,what = "both")$edges
  vertices_query<- igraph::as_data_frame(AC,what = "both")$vertices
  net_gene<- igraph::as_data_frame(GN,what = "both")$edges
  vertices_gene<- igraph::as_data_frame(GN,what = "both")$vertices
  list_vertices <- list("vertices_query" = vertices_query, "vertices_gene" = vertices_gene)
  write.xlsx(list_vertices, file = paste0(savepath,"nodes.xlsx"))
  list_net <- list("net_query" = net_query, "net_gene" = net_gene)
  write.xlsx(list_net, file = paste0(savepath,"network.xlsx"))
  
  ##节点形状
  ##创建三角形
  mytriangle <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]}
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v]}
    vertex.frame.width <- params("vertex", "frame.width")
    if (length(vertex.frame.width) != 1 && !is.null(v)) {
      vertex.frame.width <- vertex.frame.width[v]}
    vertex.size <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]}
    vertex.size <- rep(vertex.size, length.out = nrow(coords))
    vertex.frame.color[vertex.frame.width <= 0] <- NA
    vertex.frame.width[vertex.frame.width <= 0] <- 1
    if (length(vertex.frame.width) == 1) {
      symbols(x = coords[, 1], y = coords[, 2], bg = vertex.color,fg = vertex.frame.color, stars=cbind(vertex.size, vertex.size, vertex.size), lwd = vertex.frame.width,add = TRUE, inches = FALSE)
    }else {
      mapply(coords[, 1], coords[, 2], vertex.color, vertex.frame.color,vertex.size, vertex.frame.width, FUN = function(x,y, bg, fg, size, lwd) {
        symbols(x = x, y = y, bg = bg, fg = fg, lwd = lwd,stars=cbind(vertex.size, vertex.size, vertex.size),add = TRUE, inches = FALSE)})
    }
  }
  # 给网络图节点赋形状值
  add_shape("triangle", clip=shapes("circle")$clip,plot=mytriangle)
  shape<- c("circle","triangle")
  names(shape)<- c("kinase","Substrate")
  V(AC)$point.shape<- shape[match(V(AC)$type,names(shape))]
  ##给网络图的边赋颜色值
  col_up=c("#FF0000","#8B0000")
  col1 <- colorRampPalette(col_up)(length(which(E(AC)$corr>0)))
  col_down=c("#63B8FF","#00008B")
  col2 <- colorRampPalette(col_down)(length(which(!E(AC)$corr>0)))
  E(AC)$color<- c(col2,col1)
  color<- colorRampPalette(c("#00008B","#63B8FF","#FFFAFA","#FF0000","#8B0000"))(100)
  ##开始绘图
  ##
  png(file=paste0(savepath,"network_query.png"),width = 4200,height = 4200,res = 300)
  opar <- par(no.readonly = TRUE)
  par(fig = c(0, 1, 0.1, 1))
  plot.igraph(AC,layout=layout_with_fr,
              vertex.frame.color=F,vertex.shape=V(AC)$point.shape,vertex.color='#BEBEBE',vertex.frame.color='#BEBEBE',vertex.size=8,vertex.label.color="black",vertex.label.family="sans",vertex.label.cex=0.85,vertex.label.dist=1,vertex.label.degree=0,
              edge.width=8,edge.arrow.size=0,edeg.color=E(AC)$color,edge.lty=1)
  par(fig = c(0, 0.5, 0, 0.3), new = TRUE)
  colorlegend (title = 'Pearson correlation coefficient',color = color,x = c(-1, 0,1),plot = TRUE)
  legend(title = "type","top",xpd=TRUE,legend = c("kinase","Substrate"),pch=c(16,17),pt.cex=2,col=c('#BEBEBE', '#BEBEBE'))
  par(opar)
  dev.off()
  ##保存pdf格式
  pdf(file=paste0(savepath,"network_query.pdf"),width =15,height = 15)
  opar <- par(no.readonly = TRUE)
  par(fig = c(0, 1, 0.2, 1))
  plot.igraph(AC,layout=layout_with_fr,
              vertex.frame.color=F,vertex.shape=V(AC)$point.shape,vertex.color='#BEBEBE',vertex.frame.color='#BEBEBE',vertex.size=8,vertex.label.color="black",vertex.label.family="sans", vertex.label.cex=0.85, vertex.label.dist=1,vertex.label.degree=0,
              edge.width=8,edge.arrow.size=0, edeg.color=E(AC)$color, edge.lty=1)
  par(fig = c(0, 0.5, 0, 0.3), new = TRUE)
  colorlegend (title = 'Pearson correlation coefficient',color = color,x = c(-1,0,1),plot = TRUE)
  legend(title = "type","top",xpd=TRUE, legend = c("kinase","Substrate"), pch=c(16,17), pt.cex=2,col=c('#BEBEBE','#BEBEBE'))
  par(opar)
  dev.off()
  
  ##节点形状
  # 给网络图节点赋形状值
  V(GN)$point.shape<- shape[match(V(GN)$type,names(shape))]
  ##给网络图的边赋颜色值
  col1 <- colorRampPalette(col_up)(length(which(E(GN)$corr>0)))
  col2 <- colorRampPalette(col_down)(length(which(!E(GN)$corr>0)))
  E(GN)$color<- c(col2,col1)
  color<- colorRampPalette(c("#00008B","#63B8FF","#FFFAFA","#FF0000","#8B0000"))(100)
  ##开始绘图
  png(file=paste0(savepath,"network_gene.png"),width = 4200,height = 4200,res = 300)
  opar2 <- par(no.readonly = TRUE)
  par(fig = c(0, 1, 0.1, 1))
  plot.igraph(GN,layout=layout_with_fr,
              vertex.frame.color=F,vertex.shape=V(GN)$point.shape,vertex.color='#BEBEBE',vertex.frame.color='#BEBEBE',vertex.size=8,vertex.label.color="black",vertex.label.family="sans",vertex.label.cex=0.85,vertex.label.dist=1,vertex.label.degree=0,
              edge.width=8,edge.arrow.size=0,edeg.color=E(GN)$color,edge.lty=1)
  par(fig = c(0, 0.5, 0, 0.3), new = TRUE)
  colorlegend (title = 'Pearson correlation coefficient',color = color,x = c(-1, 0,1),plot = TRUE)
  legend(title = "type","top",xpd=TRUE,legend = c("kinase","Substrate"),pch=c(16,17),pt.cex=2,col=c('#BEBEBE', '#BEBEBE'))
  par(opar2)
  dev.off()
  ##保存pdf格式
  pdf(file=paste0(savepath,"network_gene.pdf"),width =15,height = 15)
  opar2 <- par(no.readonly = TRUE)
  par(fig = c(0, 1, 0.2, 1))
  plot.igraph(GN,layout=layout_with_fr,
              vertex.frame.color=F,vertex.shape=V(GN)$point.shape,vertex.color='#BEBEBE',vertex.frame.color='#BEBEBE',vertex.size=8,vertex.label.color="black",vertex.label.family="sans", vertex.label.cex=0.85, vertex.label.dist=1,vertex.label.degree=0,
              edge.width=8,edge.arrow.size=0, edeg.color=E(GN)$color, edge.lty=1)
  par(fig = c(0, 0.5, 0, 0.3), new = TRUE)
  colorlegend (title = 'Pearson correlation coefficient',color = color,x = c(-1,0,1),plot = TRUE)
  legend(title = "type","top",xpd=TRUE, legend = c("kinase","Substrate"), pch=c(16,17), pt.cex=2,col=c('#BEBEBE','#BEBEBE'))
  par(opar2)
  dev.off()  
  }

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 14, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 9, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-s","--savepath",default = "./Network/", help = "绘图存放分析文件夹路径，默认./Network/")
  parser$add_argument("-ip","--inputpath",default = "rawdata/", help = "输入文件存放路径，默认rawdata/")
  parser$add_argument("-or","--species",default = "hsa", help = "物种英文名，hsa，mmu，rno,默认hsa")
  parser$add_argument("-phos","--phosfile",default = "表达矩阵-ptm.xlsx", help = "磷酸化表达矩阵,默认表达矩阵-ptm.xlsx")
  parser$add_argument("-pro","--profile",default = "表达矩阵-pro.xlsx", help = "蛋白表达矩阵,默认表达矩阵-pro.xlsx")
  parser$add_argument("-c","--corr",default=0.8, type= "double",help="激酶-底物间相关性，默认为0.8")
  args <- parser$parse_args()
  KSEA_corr <- do.call(KSEA_corr,args = args)
}
  