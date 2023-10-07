#!/opt/conda/bin/Rscript
#' 蛋白互作网络图
#'
#' @param savepath 
#' @param inputpath 
#' @param inputfile 
#' @param col 
#' @param noncol 
#' @param compare 
#' @param backpath 
#' @param backfile 
#' @param number 
#' @param fc 
#' @param imagetype 
#' @param height 
#' @param width 
#' @param dpi 
#' @param fontfamily 
#'
#' @export
#'
ppiplot <- function(savepath = "./ppi/",inputpath="./",inputfile = "A-B-diff-protein.xls",
                col=NULL,noncol="#4f94cd",df="D",
                backpath="../background/",backfile="protein2protein_network.xls",
                number=25,imagetype = c("pdf","png") ,height = 10,width = 15,dpi=300,fontfamily="sans",...){
  pacman::p_load(dplyr,tidyr,igraph,ggraph,stringr,patchwork)
  #互作网络关系表格
  compare<-gsub("-diff-.*","",inputfile)
  ppi_all <- readdata(paste0(backpath,"/",backfile))
  diff <- readdata(paste0(inputpath,"/",inputfile))
  diff<-diff[!(duplicated(diff[,1])),] %>% as.data.frame()
  if(!is.na(grep(":",diff[1,1])[1])){
    ppi_fc<-data.frame("Accession"=gsub(":.*","",diff[,1]),"Gene Name"=gsub(".*:","",diff[,1]),check.names = F)
  }else  ppi_fc<-data.frame("Accession"=diff[,1])
  ppi_fc["Label"]<-diff[,1]
  if(ncol(diff)!=1){
    if("FC" %in% names(diff)){
      ppi_fc[,"FoldChange"]<-diff[,"FC"]
    }else{
      ppi_fc[,"FoldChange"]<-diff[,"FoldChange"]
    }
    
  }
  ll<-ppi_all[ppi_all[,1]%in%ppi_fc[,1],]
  proppi<-ll[ll[,2]%in%ppi_fc[,1],]
  if(nrow(proppi)!=0){
    ppi_fc[,"Degree"]<-sapply(1:nrow(ppi_fc),function(x){
      grep(pattern=ppi_fc[x,1],x = c(proppi[,1],proppi[,2]))%>%length()
    })
    colnames(proppi) <- c("node1","node2","combined_score")
    savexls(proppi,paste0(savepath,compare,"/ppi_network.xls"))
    savexls(ppi_fc,paste0(savepath,compare,"/ppi_nodes.xls"))
    #有无蛋白
    if("FoldChange" %in% names(ppi_fc)){
      tmp <- ppi_fc[which(ppi_fc$FoldChange != "inf" & ppi_fc$FoldChange != "0"), ]
      max <- max(as.numeric(tmp$FoldChange))
      min <- min(as.numeric(tmp$FoldChange))
      ppi_fc$FoldChange[ppi_fc$FoldChange=="0"] <- min
      ppi_fc$FoldChange[ppi_fc$FoldChange=="inf"] <- max
      ppi_fc$FoldChange <- as.numeric(ppi_fc$FoldChange)
    }
    
    #选取连接度前number(默认25)展示——修饰去重后
    if(nrow(ppi_fc[ppi_fc$Degree!=0,])>=number){
      ppid<-ppi_fc[ppi_fc$Degree!=0,]
      ppi_fc_top <- ppid[order(ppid$Degree,decreasing = T),][1:number,]
      ll_top<-proppi[proppi[,1]%in%ppi_fc_top[,1],]
      proppi_top<-ll_top[ll_top[,2]%in%ppi_fc_top[,1],]
    }else{
      ppid <- ppi_fc[ppi_fc$Degree!=0,]
      ppi_fc_top <- ppid[order(ppid$Degree,decreasing = T),]
      proppi_top <- proppi
    }
    
    #绘图
    #没有基因名的蛋白填充
    if("Gene Name" %in% names(ppi_fc_top)){
      ppi_fc_top[,"Gene Name"] <- apply(ppi_fc_top,1,function(x){ifelse(x[2]==""|x[2]=="  "|is.na(x[2]),x[2] <- x[1],x[2] <- x[2])})
    }
    
    colnames(proppi_top) <- c("node1","node2","combined_score")
    if(is.null(col)){
      my_color <-c("#0000A1","#003399","#1F6ED4","azure", "#f47d8c","#DE3F2E", "#A50026")
    }else{
      my_color <- col
    }
    
    #基因名和蛋白名展示
    ppilabel <- ppi_fc_top[,2]
    names(ppilabel) <- ppi_fc_top[,1]
    #展示蛋白重新计算连接度
    list1 <- as.character(proppi_top$node1)
    list2 <- as.character(proppi_top$node2)
    list <- as.character(cbind(list1, list2))
    df_num <- as.data.frame(table(list))
    names(df_num) <- c("Accession", "degree")
    gene2Num <- df_num[order(df_num[, 2], decreasing = T), ]
    gene2FC <- select(ppi_fc_top,c("Accession",grep("FoldChange|Gene Name",names(ppi_fc_top),value = T),"Label"))
    
    mergr_nodes <- merge(
      gene2Num, gene2FC,
      by= "Accession", all.x = T)
    
    if("FoldChange" %in% names(ppi_fc)){
      #取上下调
      df_Up <- mergr_nodes[which(mergr_nodes[, 3] > 1), ]
      df_Down <- mergr_nodes[which(mergr_nodes[, 3] < 1), ]
      df_Up <- df_Up[order(df_Up[, 2], decreasing = F), ]
      df_Down <- df_Down[order(df_Down[, 2], decreasing = F), ]
      my_nodes <- rbind(df_Up, df_Down)
    }else{
      #无上下调
      df_all <- mergr_nodes[order(mergr_nodes[, 2], decreasing = F), ]
      my_nodes <- df_all
      my_color <- noncol
    }
    
    if(length(colnames(proppi_top))>=3){
      my_link <- proppi_top[,c("node1","node2","combined_score")]
      my_link <- my_link %>% drop_na()
    }else{
      my_link <- proppi_top[,c("node1","node2")]
      my_link <- my_link %>% drop_na()
    }
    
    net <- graph.data.frame(my_link, my_nodes, directed = F)
    
    deg <- igraph::degree(net, mode = "all")
    
    V(net)$deg <- deg
    p0 <- ggraph(net, layout = "kk")
    if("FoldChange" %in% names(my_nodes)){
      pp <- p0 + geom_edge_arc(
        color = "grey",
        strength = 0.1, alpha = 1,width=0.5,
        end_cap = circle(3, "mm"), start_cap = circle(3, "mm")
      ) +
        geom_node_point(
          pch = 19,
          aes(color = log2(V(net)$FoldChange+0.001), size = deg)
        ) +
        geom_node_text(
          aes(label = V(net)$name),
          size = 4,
          repel = TRUE,
          max.overlaps = Inf
        ) +
        scale_size_continuous(breaks = bubblescale(my_nodes$degree)[[1]], range = bubblescale(my_nodes$degree)[[2]]+4,name = "Degree")+
        scale_colour_gradientn(
          expression(paste(log[2], " FoldChange")),
          limits = c(
            -max(abs(log2(V(net)$FoldChange+0.001))),
            max(abs(log2(V(net)$FoldChange+0.001)))
          ),
          colours = my_color
        ) +
        guides(size = guide_legend(order = 1)) +
        theme_void() + 
        labs(title = compare)+
        theme(legend.text=element_text(size=10),legend.title=element_text(size=12),
              plot.margin = margin(t = 20,  # 顶部边缘距离
                                   r = 0,  # 右边边缘距离
                                   b = 40,  # 底部边缘距离
                                   l = 20)) # 左边边缘距离
    }else{
      pp <- p0 + geom_edge_arc(
        color = "grey",
        strength = 0.1, alpha = 1,width=0.5,
        end_cap = circle(3, "mm"), start_cap = circle(3, "mm")
      ) +
        geom_node_point(
          pch = 19,
          color = noncol, aes(size = deg)
        ) +
        geom_node_text(
          aes(label = V(net)$name),
          size = 4,
          repel = TRUE,
          max.overlaps = Inf
        ) +
        scale_size_continuous(breaks = bubblescale(my_nodes$degree)[[1]], range = bubblescale(my_nodes$degree)[[2]]+4,name = "Degree")+
        guides(size = guide_legend(order = 1)) +
        theme_void() + 
        labs(title = compare)+
        theme(legend.text=element_text(size=10),legend.title=element_text(size=12),
              plot.margin = margin(t = 20,  # 顶部边缘距离
                                   r = 0,  # 右边边缘距离
                                   b = 40,  # 底部边缘距离
                                   l = 20)) # 左边边缘距离
    }
    if("Gene Name" %in% names(my_nodes)){
      if("FoldChange" %in% names(my_nodes)){
        pg <- p0 + geom_edge_arc(
          color = "grey",
          strength = 0.1, alpha = 1,width=0.5,
          end_cap = circle(3, "mm"), start_cap = circle(3, "mm")
        ) +
          geom_node_point(
            pch = 19,
            aes(color = log2(V(net)$FoldChange+0.001), size = deg)
          ) +
          geom_node_text(
            aes(label = ppilabel[V(net)$name]),
            size = 4,
            repel = TRUE,
            max.overlaps = Inf
          ) +
          scale_size_continuous(breaks = bubblescale(my_nodes$degree)[[1]], range = bubblescale(my_nodes$degree)[[2]]+4,name = "Degree")+
          scale_colour_gradientn(
            expression(paste(log[2], " FoldChange")),
            limits = c(
              -max(abs(log2(V(net)$FoldChange+0.001))),
              max(abs(log2(V(net)$FoldChange+0.001)))
            ),
            colours = my_color
          ) +
          guides(size = guide_legend(order = 1)) +
          theme_void() + 
          labs(title = compare)+
          theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.title=element_text(size=15),
                plot.margin = margin(t = 20,  # 顶部边缘距离
                                     r = 0,  # 右边边缘距离
                                     b = 40,  # 底部边缘距离
                                     l = 20)) # 左边边缘距离
      }else{
        pg <- p0 + geom_edge_arc(
          color = "grey",
          strength = 0.1, alpha = 1,width=0.5,
          end_cap = circle(3, "mm"), start_cap = circle(3, "mm")
        ) +
          geom_node_point(
            pch = 19,
            aes(size = deg),color=noncol
          ) +
          geom_node_text(
            aes(label = ppilabel[V(net)$name]),
            size = 4,
            repel = TRUE,
            max.overlaps = Inf
          ) +
          scale_size_continuous(breaks = bubblescale(my_nodes$degree)[[1]], range = bubblescale(my_nodes$degree)[[2]]+4,name = "Degree")+
          guides(size = guide_legend(order = 1)) +
          theme_void() + 
          labs(title = compare)+
          theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.title=element_text(size=15),
                plot.margin = margin(t = 20,  # 顶部边缘距离
                                     r = 0,  # 右边边缘距离
                                     b = 40,  # 底部边缘距离
                                     l = 20)) # 左边边缘距离
      }
    }
    
    #前number(默认25)连接度绘制棒棒糖图
    #top25 <- df_num
    if(ncol(diff)!=1){
      scale_color <- scale_color_gradient2(low=my_color[1:3], high=my_color[5:7], midpoint = 0)
    }else{
      ppi_fc_top$FoldChange <- 1
      scale_color <- scale_colour_gradientn(colours = noncol)
    }
    if(df=="D"){
      ppi_top<-left_join(ppi_fc_top,df_num,by="Accession")
      ppi_top <- ppi_top[order(ppi_top$degree,decreasing = T),]
      ppi_top$Label<-factor(ppi_top$Label,levels =rev(ppi_top$Label))
      
      p2 <- ggplot(ppi_top, aes(y=Label, x=degree, color = log2(FoldChange))) +
        geom_segment(aes(x = 0,
                         y = Label,
                         xend = degree,
                         yend = Label),
                     size=2) +
        geom_point(aes(size=degree),stat='identity') +
        #geom_text(aes(label=Degree),hjust=-0.5, vjust=0.5,size=3,color="black") +
        labs(title="Top protein rank") + xlab("Degree")+ylab("")+
        scale_color+
        xlim(0,max(ppi_top$degree)*1.2)+
        theme_bw()+
        theme(panel.grid=element_blank())+
        theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
        theme(axis.text.y=element_text(size=10,color="black"),
              axis.text.x=element_text(size=10,color="black"))+
        theme(plot.title = element_text(hjust = 0.5, size=12))+
        guides(colour = "none",size="none")+
        theme(plot.margin=unit(c(8,2,8,0), "lines"))
    }else{
      ppi_top<-left_join(ppi_fc_top,df_num,by="Accession")
      ppi_top <- ppi_top[order(ppi_top$FoldChange,decreasing = T),]
      ppi_top$Label<-factor(ppi_top$Label,levels =rev(ppi_top$Label))
      
      p2 <- ggplot(ppi_top, aes(y=Label, x=log2(FoldChange), color = log2(FoldChange))) +
        geom_segment(aes(x = 0,
                         y = Label,
                         xend = log2(FoldChange),
                         yend = Label),
                     size=2) +
        geom_point(aes(size=log2(FoldChange)),stat='identity') +
        labs(title="Top protein rank") + xlab("log2(FoldChange)")+ylab("")+
        scale_color+
        xlim(0,max(abs(log2(ppi_top$FoldChange)))*1.2)+
        theme_bw()+
        theme(panel.grid=element_blank())+
        theme(panel.border = element_rect(fill=NA,color="black", linewidth=1))+
        theme(axis.text.y=element_text(size=10,color="black"),
              axis.text.x=element_text(size=10,color="black"))+
        theme(plot.title = element_text(hjust = 0.5, size=12))+
        guides(colour = "none",size="none")+
        theme(plot.margin=unit(c(8,2,8,0), "lines"))
    }
    
    savexls(ppi_top,paste0(savepath,compare,"/Top_nodes.xls"))
    savexls(my_link,paste0(savepath,compare,"/Top_network.xls"))
    proplot <- cowplot::plot_grid(pp, p2, nrow = 1,ncol=3,rel_widths = c(2.5,1,0.1))
    ggplotsave(plot = proplot,savepath=paste0(savepath,compare,"/"),
               mapname = paste0("ppi_query"),
               width = width,
               height = height,
               imagetype = imagetype,
               family=fontfamily,...)
    if("Gene Name" %in% names(my_nodes)){
      geneplot <- cowplot::plot_grid(pg, p2, nrow = 1,ncol=3,rel_widths = c(2.5,1,0.1))
      ggplotsave(plot = geneplot,savepath=paste0(savepath,compare,"/"),
                 mapname = paste0("ppi_gene"),
                 width = width,
                 height = height,
                 imagetype = imagetype,
                 family=fontfamily,...)
    }
    
  }else{
    savetxt(data = "所选蛋白没有在数据库中匹配到信息",
            filename = paste0(savepath,compare,"/说明.txt"),append = T)
    return()
  }
  
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("png","pdf"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 15, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 10, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-if","--inputfile",default = "A-B-diff-protein.xls", help = "输入数据文件，数据文件名必须包含-diff-")
  parser$add_argument("-c","--col", default=NULL, help="有上下调时颜色")
  parser$add_argument("-df","--df", default="D", help="ppi右侧棒棒糖图横坐标D/F（degree/FoldChange)，默认为D")
  parser$add_argument("-nc","--noncol", default="#4f94cd", help="无上下调时颜色，默认蓝色")
  parser$add_argument("-bp","--backpath", default="../background/", help="背景路径，默认../background/")
  parser$add_argument("-bf","--backfile", default="protein2protein_network.xls", help="背景文件，默认protein2protein_network.xls")
  parser$add_argument("-n","--number", type="integer",default=25, help="展示连接度高的蛋白互作图，默认连接度前25")
  parser$add_argument("-s","--savepath",default = "./ppi/", help = "结果保存路径,默认./ppi/")
  parser$add_argument("-ip","--inputpath", default="./", help="分析路径")
  args <- parser$parse_args()
  ppi_plot <- do.call(ppiplot,args = args)
}
