#!/opt/conda/bin/Rscript
#' 富集环状图
#' @param savepath 保存总富集路径
#' @param type 数据库类型缩写
#' @param incompare 比较组名称
#' @param filt 数据是否进行筛选
#' @param number 绘图top数量
#' @param imagetype 保存图片类型
#' @param height 高度
#' @param width 宽度
#' @param dpi 分辨率
#' @param fontfamily 字体类型
#' @export
circleplot<-function(savepath = "./enrich/",type = "K",inputpath="./",inputfile = NULL,incompare = "A_B", filt = "T",
                     height = 13,width = 15,dpi=300,fontfamily="sans",...){
  # 加载包
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(ComplexHeatmap))
  
  mycol <- c("#FF7F0E","#9068be","#2CA03A","#FF7F0E","#9068be","#2CA03A","#D62728","#17BECF","#D11250","#8C564B")
  gklevel1=c("biological_process","cellular_component","molecular_function",
             "Metabolism","Genetic Information Processing","Environmental Information Processing",
             "Cellular Processes","Organismal Systems","Human Diseases","Drug Development")
  names(mycol)<-gklevel1
  comenrichpath<-paste0(savepath,getenrichtype(type)$filn,"/",incompare)
  enrichfile<-dir(path =comenrichpath, pattern = "enrichment-.*-Total.*")
  
  if(is.null(inputfile)){
    degFile<-dir(path =inputpath, pattern = paste0("^",incompare,"-diff-.*"))[1]
  }else{
    degFile<-basename(inputfile)
  }
  
  enrich_dat <- read.delim(paste0(comenrichpath,"/",enrichfile), sep="\t", header=T, quote="",check.names=FALSE)
  enrich_dat$`p-value`[enrich_dat$`p-value`==0]<-min(enrich_dat$`p-value`[enrich_dat$`p-value`!=0])
  deg_dat <- read.delim(paste0(inputpath,"/",degFile), sep="\t", header=T, quote="",check.names=FALSE)
  deg_dat <- deg_dat[!duplicated(deg_dat[,1]),]
  deg_dat[,1]<-gsub(":.*","",deg_dat[,1]) %>% as.data.frame()
  names(deg_dat)[1]<-"id"
  pg<-unique(unlist(strsplit(as.character(enrich_dat$Substances),",")))
  pgn<-data.frame("id"=gsub(":.*","",pg),pg)
  deg_dat<-left_join(pgn,deg_dat,by="id")[,-1]%>% as.data.frame()
  names(deg_dat)[1]<-"id"
  rownames(deg_dat)<-deg_dat[,1]
  if (ncol(deg_dat) < 3 || nrow(enrich_dat) < 10) {
    savetxt(data = "Term数量低于10或没有上下调信息时，不提供富集条目环状图",
            filename = paste0(comenrichpath,"/说明.txt"),append = T)
    return()
  }else{
    if (filt == "T"){
      enrich_dat <- head(enrich_dat, 20)
      savexls(enrich_dat,paste0(comenrichpath,"/",getenrichtype(type)$filn,".circos.xls"))
    }
    
    enrich_dat <- enrich_dat[
      order(enrich_dat$ListHits / enrich_dat$PopHits, decreasing = T), , drop = F
    ]
    
    reg_count <- function(x, deg_dat, flag) {
      length(
        which(
          unlist(strsplit(x, ",")) %in%
            rownames(deg_dat)[deg_dat$type == flag]
        )
      )
    }
    
    
    dat <- data.frame(
      category = enrich_dat[
        ,
        grepl("Category|Classification_level1", colnames(enrich_dat)),
        drop = T
      ],
      gene_num.min = 0,
      gene_num.max = max(enrich_dat$PopHits),
      gene_num.rich = enrich_dat$PopHits,
      log.p = -log10(enrich_dat[["p-value"]]),
      rich.factor = enrich_dat$ListHits / enrich_dat$PopHits,
      up.regulated = vapply(
        enrich_dat$Substances,
        reg_count,
        deg_dat = deg_dat, flag = "Up",
        0,
        USE.NAMES = F
      ),
      down.regulated = vapply(
        enrich_dat$Substances,
        reg_count,
        deg_dat = deg_dat, flag = "Down",
        0,
        USE.NAMES = F
      )
    )
    circplot<-function(...){
      circos.par(
        gap.degree = 2, start.degree = 90, circle.margin = c(1e-5, .8, 1e-5, 1e-5)
      )
      
      ## 第一圈
      # 选择作图数据集，定义区块的基因总数量范围
      plot_data <- dat[c("id", "gene_num.min", "gene_num.max")]
      # 定义分组颜色
      cat_color <- mycol[dat$category]
      cat_col <-cat_color[!duplicated(cat_color)]
      # 一个总布局
      circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)
      circos.track(
        # 圈图的高度、颜色等设置
        ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = cat_color,
        panel.fun = function(x, y) {
          # ylim、xlim 用于指定 id 文字标签添加的合适坐标
          ylim <- get.cell.meta.data("ycenter")
          xlim <- get.cell.meta.data("xcenter")
          # sector.name 用于提取 id 名称
          sector.name <- get.cell.meta.data("sector.index")
          # 绘制外周的刻度线
          circos.axis(
            h = "top", labels.cex = 0.8, major.tick.length = 0.8,
            labels.niceFacing = FALSE
          )
          # 将 id 文字标签添加在图中指定位置处
          circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)
        }
      )
      
      ## 第二圈，绘制富集的基因和富集 p 值
      # 选择作图数据集，包括富集基因数量以及 p 值等信息
      plot_data <- dat[c("id", "gene_num.min", "gene_num.rich", "log.p")]
      # 标签数据集，仅便于作图时添加相应的文字标识用
      label_data <- dat["gene_num.rich"]
      # 定义一个 p 值的极值，以方便后续作图
      p_max <- round(max(dat$log.p)) + 1
      # 这两句用于定义 p 值的渐变颜色
      RdYlBu <- rev(SelectColors("brewer_celsius"))
      colorsChoice <- colorRampPalette(rev(RdYlBu))
      color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))
      
      circos.genomicTrackPlotRegion(
        plot_data,
        # 圈图的高度、颜色等设置
        track.height = 0.08, bg.border = NA, stack = TRUE,
        panel.fun = function(region, value, ...) {
          # 区块的长度反映了富集基因的数量，颜色与 p 值有关
          circos.genomicRect(
            region, value,
            col = color_assign(value[[1]]), border = NA, ...
          )
          # 指定文字标签（富集基因数量）添加的合适坐标
          ylim <- get.cell.meta.data("ycenter")
          xlim <- label_data[get.cell.meta.data("sector.index"), 1] / 2
          sector.name <- label_data[get.cell.meta.data("sector.index"), 1]
          # 将文字标签添（富集基因数量）加在图中指定位置处
          circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)
        }
      )
      
      ## 第三圈，绘制上下调基因
      # 首先基于表格中上下调基因的数量，计算它们的占比
      dat$all.regulated <- dat$up.regulated + dat$down.regulated
      dat$up.proportion <- dat$up.regulated / dat$all.regulated
      dat$down.proportion <- dat$down.regulated / dat$all.regulated
      
      # 随后，根据上下调基因的相对比例，分别计算它们在作图时的“区块坐标”和“长度”
      dat$up <- dat$up.proportion * dat$gene_num.max
      plot_data_up <- dat[c("id", "gene_num.min", "up")]
      names(plot_data_up) <- c("id", "start", "end")
      plot_data_up$type <- 1 # 分配 1 指代上调基因
      
      dat$down <- dat$down.proportion * dat$gene_num.max + dat$up
      plot_data_down <- dat[c("id", "up", "down")]
      names(plot_data_down) <- c("id", "start", "end")
      plot_data_down$type <- 2 # 分配 2 指代下调基因
      
      # 选择作图数据集、标签数据集，并分别为上下调基因赋值不同颜色
      plot_data <- rbind(plot_data_up, plot_data_down)
      label_data <- dat[c("up", "down", "up.regulated", "down.regulated")]
      color_assign <- colorRamp2(breaks = c(1, 2), col = c("#f47d8c", "#727fb5"))
      
      # 继续绘制圈图
      suppressMessages(circos.genomicTrackPlotRegion(
        plot_data,
        # 圈图的高度、颜色等设置
        track.height = 0.08, bg.border = NA, stack = TRUE,
        panel.fun = function(region, value, ...) {
          circos.genomicRect(
            # 区块的长度反映了上下调基因的相对占比
            region, value,
            col = color_assign(value[[1]]), border = NA, ...
          )
          # 指定文字标签（上调基因数量）添加的合适坐标
          ylim <- get.cell.meta.data(
            "cell.bottom.radius"
          ) - 0.5
          xlim <- label_data[get.cell.meta.data("sector.index"), 1] / 2
          sector.name <- label_data[get.cell.meta.data("sector.index"), 3]
          # 将文字标签（上调基因数量）添加在图中指定位置处
          circos.text(
            xlim, ylim, sector.name,
            cex = 0.8, niceFacing = FALSE
          )
          xlim <- (label_data[get.cell.meta.data("sector.index"), 2] +
                     label_data[get.cell.meta.data("sector.index"), 1]) / 2
          sector.name <- label_data[get.cell.meta.data("sector.index"), 4]
          # 类似的操作，将下调基因数量的标签也添加在图中
          circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)
        }
      ))
      
      ## 第四圈，绘制富集得分
      # 选择作图数据集，标准化后的富集得分
      plot_data <- dat[
        c("id", "gene_num.min", "gene_num.max", "rich.factor")
      ]
      # 将通路的分类信息提取出，和下一句一起，便于作图时按分组分配颜色
      label_data <- dat["category"]
      
      circos.genomicTrack(
        plot_data,
        # 圈图的高度、颜色等设置
        ylim = c(0, 1), track.height = 0.4, bg.col = "gray95", bg.border = NA,
        panel.fun = function(region, value, ...) {
          # sector.name 用于提取 id 名称，并添加在下一句中匹配的高级分类，以分配颜色
          sector.name <- get.cell.meta.data("sector.index")
          # 等位线
          # circos.lines(c(0, max(region)), c(0.1, 0.1), col = "gray", lwd = 0.3)
          circos.lines(c(0, max(region)), c(0.2, 0.2), col = "gray", lwd = 0.3)
          # circos.lines(c(0, max(region)), c(0.3, 0.3), col = "gray", lwd = 0.3)
          circos.lines(c(0, max(region)), c(0.4, 0.4), col = "gray", lwd = 0.3)
          # circos.lines(c(0, max(region)), c(0.5, 0.5), col = "gray", lwd = 0.3)
          circos.lines(c(0, max(region)), c(0.6, 0.6), col = "gray", lwd = 0.3)
          # circos.lines(c(0, max(region)), c(0.7, 0.7), col = "gray", lwd = 0.3)
          circos.lines(c(0, max(region)), c(0.8, 0.8), col = "gray", lwd = 0.3)
          # circos.lines(c(0, max(region)), c(0.9, 0.9), col = "gray", lwd = 0.3)
          # 绘制矩形区块，高度代表富集得分，颜色代表分类
          circos.genomicRect(
            region, value,
            col = cat_col[label_data[sector.name, 1]],
            border = NA, ytop.column = 1, ybottom = 0, ...
          )
        }
      )
      
      polygon(
        x = c(-0.2, -0.2, -0.05, -0.05),
        y = c(0.12, 0.18, 0.18, 0.12),
        col = "grey", border = NA
      )
      text(-0.075, 0.15, "  Number ", cex = 1, pos = 4)
      text(-0.125, 0.15, "", cex = 1)
      
      polygon(
        x = c(-0.2, -0.2, -0.05, -0.05),
        y = c(0.02, 0.08, 0.08, 0.02),
        col = "#f47d8c", border = NA
      )
      text(-0.075, 0.05, "  Up-regulated", cex = 1, pos = 4)
      
      polygon(
        x = c(-0.2, -0.2, -0.05, -0.05),
        y = c(-0.08, -0.02, -0.02, -0.08),
        col = "#727fb5", border = NA
      )
      text(-0.075, -0.05, "  Down-regulated", cex = 1, pos = 4)
      
      polygon(
        x = c(-0.2, -0.2, -0.05, -0.05),
        y = c(-0.125, -0.175, -0.19, -0.11),
        col = "grey", border = NA
      )
      text(-0.075, -0.15, "  Rich Factor [0-1)", cex = 1, pos = 4)
      
      ## 绘图完毕后，不要忘了清除痕迹，以免影响下一次作图
      circos.clear()
      
      
      category_legend <- Legend(
        labels = names(cat_col),
        type = "points", pch = NA,
        background = cat_col,
        labels_gp = gpar(fontsize = 14),
        grid_height = unit(0.5, "cm"),
        grid_width = unit(0.5, "cm")
      )
      
      pvalue_legend <- Legend(
        col_fun = colorRamp2(
          round(seq(0, p_max, length.out = 6), 0),
          colorRampPalette(rev(RdYlBu))(6)
        ),
        legend_height = unit(3, "cm"),
        labels_gp = gpar(fontsize = 14),
        title_gp = gpar(fontsize = 16),
        title_position = "topleft",
        title = expression(paste("-", log[10], "p-value"))
      )
      
      lgd_list_vertical <- packLegend(category_legend, pvalue_legend)
      pushViewport(viewport(x = 0.85, y = 0.5))
      grid.draw(lgd_list_vertical)
      upViewport()
    }
    
    # 默认按照原表格中的排列顺序
    dat$id <- factor(enrich_dat[["id"]], levels = enrich_dat[["id"]])
    rownames(dat) <- dat$id
    
    #### 创建一个 pdf 画板
    pdf(
      file.path(comenrichpath, paste0(getenrichtype(type)$filn, ".circos.pdf")),
      width = width, height = height,family=fontfamily
    )
    ## 整体布局
    circplot()
    dev.off()
    #### 创建一个 png 画板
    png(
      file.path(comenrichpath, paste0(getenrichtype(type)$filn, ".circos.png")),
      width = width*200, height = height*200,res=200,family=fontfamily
    )
    ## 整体布局
    circplot()
    dev.off()
  }
}


if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 15, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 13, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-t","--type",type="character", default="K", help="enrich database class,such as G/K/W/R/I", metavar="character")
  parser$add_argument("-ic","--incompare",default = "A_B", help = "比较组名称，默认A_B")
  parser$add_argument("-f","--filt",default = "T", type= "character",help = "数据是否需要筛选,默认T")
  parser$add_argument("-ip","--inputpath",default = "./", help = "-diff-文件所在路径，默认当前路径")
  parser$add_argument("-sp","--savepath",default = "./enrich/", help = "富集结果保存路径,默认./enrich/")

  args <- parser$parse_args()
  enrich_circleplot <- do.call(circleplot,args = args)
}
