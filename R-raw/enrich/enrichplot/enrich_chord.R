#!/opt/conda/bin/Rscript

#' GO_chord plot
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
#' @param ... 
#' @export
mainchord <- function(savepath = "./enrich/",type = "K",enrichfile=NULL ,inputpath="./",inputfile = NULL,incompare = NULL, filt = "T",number = 10,imagetype = c("pdf","png") ,
                      height = 10,width = 8,dpi=300,fontfamily="sans",...) {
  pacman::p_load(dplyr,ggplot2,GOplot,stringr,Hmisc)
  environment(enrich_chord) <- environment(GOChord)
  if (is.null(enrichfile)){
    comenrichpath<-paste0(savepath,getenrichtype(type)$filn,"/",incompare)
  }else{
    comenrichpath<-'./'
  }
  enrichfile<-dir(path =comenrichpath, pattern = "enrichment-.*-Total.*")
  
  if(is.null(inputfile)){
    degFile<-dir(path =inputpath, pattern = paste0("^",incompare,"-diff-.*"))[1]
  }else{
    degFile<-basename(inputfile)
  }
  topDf <- readdata(paste0(comenrichpath,"/",enrichfile))
  
  if (getenrichtype(type)$filn == "GO") {
    ctg <- "Category"
    term <- "Term"
  } else if (getenrichtype(type)$filn == "KEGG") {
    ctg <- "Classification_level1"
    term <- "Term"
  } else {
    ctg <- "Term"
    term <- "Term"
  }
  if(!'p-value'%in%colnames(topDf)){
    print('~请检查p-value列及列名！')
    return()
  }else if(!'ListHits'%in%colnames(topDf)){
    print('~请检查ListHits列及列名！')
    return()
  }else if(!'Substances'%in%colnames(topDf)){
    print('~请检查Substances列及列名！')
    return()
  }else if(!term%in%colnames(topDf)){
    print('~请检查Term列及列名！')
    return()
  }else if(!ctg%in%colnames(topDf)){
    print(paste0('~请检查',ctg,'列及列名！'))
    return()
  }else{
    topDf$`p-value`[topDf$`p-value`==0]<-min(topDf$`p-value`[topDf$`p-value`!=0])
    if(is.na(degFile)){
      degDf<- unique(unlist(strsplit(as.character(topDf$Substances),",")))%>%as.data.frame()
    }else{
      degDf <- readdata(paste0(inputpath,"/",degFile))
    }
    degDf <- degDf[!duplicated(degDf[,1]),]%>%as.data.frame()
    degDf[,1]<-gsub(":.*","",degDf[,1]) %>% as.data.frame()
    names(degDf)[1]<-"id"
    pg<-unique(unlist(strsplit(as.character(topDf$Substances),",")))
    pgn<-data.frame("id"=gsub(":.*","",pg),pg)
    degDf<-left_join(pgn,degDf,by="id")[,-1]%>% as.data.frame()
    names(degDf)[1]<-"id"
    if (filt == "T"){
      topDf <- filter(topDf,ListHits>=3)[1:number,] %>% na.omit()
      savexls(topDf,paste0(comenrichpath,"/",getenrichtype(type)$filn,".chord.xls"))
    }
    if (nrow(topDf) < 3) {
      savetxt(data = "Term数量低于3，不提供富集条目和弦图",
              filename = paste0(comenrichpath,"/说明.txt"),append = T)
      return()
    }else{
      chord <- make_chord_data(topDf, degDf, db = getenrichtype(type)$filn,number = number)
      chordPlot <- enrich_chord(chord,degDf=degDf,gene.order = "logFC")
      ggplotsave(plot = chordPlot,savepath=comenrichpath,
                 mapname = paste0(getenrichtype(type)$filn,".chord"),
                 width = width,
                 height = height,
                 imagetype = imagetype,
                 family=fontfamily)
  }
      }
}
# parse input table, make chord data for GOplot
#' @export
make_chord_data <- function(topDf, degDf, db = "KEGG",number) {
  if (db == "GO") {
    ctg <- "Category"
    n <- unique(topDf[, ctg, drop = T])
    term <- "Term"
  } else if (db == "KEGG") {
    ctg <- "Classification_level1"
    n <- gsub(" ", "_", unique(topDf[, ctg, drop = T]))
    term <- "Term"
  } else {
    ctg <- "Term"
    n <- unique(topDf[, ctg, drop = T])
    term <- "Term"
  }
  if("FoldChange" %in% colnames(degDf)){
    geneDf <- data.frame(ID = degDf$id, logFC = degDf$FoldChange)
  }else if("FC" %in% colnames(degDf)){
    geneDf <- data.frame(ID = degDf$id, logFC = degDf$FC)
  }else geneDf <- data.frame(ID = degDf$id, logFC=runif(length(degDf$id),0.1,10))
  
  termDf <- topDf[, c(ctg, "id", term, "p-value", "Substances")]
  colnames(termDf) <- c(
    "category", "ID", "term", "adj_pval", "genes"
  ) 
  termDf$term <- Hmisc::capitalize(as.character(termDf$term))
  termDf$genes <- gsub(";", ", ", termDf$genes)
  termDf$category <- gsub(" ", "_", termDf$category) # KEGG
  
  # p值排序，最多 10 term
  termDf <- head(termDf[order(termDf$adj_pval, decreasing=F), , drop = F], number)
  
  circDf <- circle_dat(termDf, geneDf)
  
  # fix Inf value
  geneDf$logFC[geneDf$logFC == "0"] <-
    min(as.numeric(geneDf$logFC[geneDf$logFC!="inf"&geneDf$logFC!="0"]))
  geneDf$logFC[geneDf$logFC == "inf"] <-
    max(as.numeric(geneDf$logFC[geneDf$logFC!="inf"]))
  geneDf$logFC<-log2(geneDf$logFC)
  geneDf$ID <- toupper(geneDf$ID) # bug fix, same with circle_dat
  
  # most 10 genes in each term
  subgene <- function(geneDf, circDf) {
    select_genes <- function(x) {
      subg <- geneDf[
        which(geneDf$ID %in% circDf$genes[circDf$term == x]), ,
        drop = F
      ]
      subg <- subg[
        head(order(abs(subg$logFC), decreasing = T), 10), , drop = F
      ]
    }
    
    subg <- do.call(rbind, lapply(unique(circDf$term), select_genes))
    subg <- subg[!duplicated(subg), ]
    
    subc <- circDf[which(circDf$genes %in% subg$ID), , drop = F]
    
    return(list(subg, subc))
  }
  
  chordDf <- chord_dat(
    data = subgene(geneDf, circDf)[[2]],
    genes = subgene(geneDf, circDf)[[1]],
    process = unique(
      subgene(geneDf, circDf)[[2]]$term
    )
  )
  
  return(chordDf)
}
#' @export
check_chord <- function(mat, limit){
  
  if(all(colSums(mat) >= limit[2]) & all(rowSums(mat) >= limit[1])) return(mat)
  
  tmp <- mat[(rowSums(mat) >= limit[1]),]
  mat <- tmp[,(colSums(tmp) >= limit[2])]
  
  mat <- check_chord(mat, limit)
  return(mat)
}

# modified GOChord. colour, border, label default size ...
#' @export
enrich_chord <- function(data,degDf=degDf, title, space, gene.order, gene.size, gene.space,
                         nlfc = 1, lfc.col, lfc.min, lfc.max, ribbon.col,
                         border.size, process.label, limit) {
  y <- id <- xpro <- ypro <- xgen <- ygen <- lx <- ly <- ID <- logFC <- NULL
  Ncol <- dim(data)[2]
  if (missing(title)) {
    title <- ""
  }
  if (missing(space)) {
    space <- 0
  }
  if (missing(gene.order)) {
    gene.order <- "none"
  }
  if (missing(gene.size)) {
    gene.size <- 2
  }
  if (missing(gene.space)) {
    gene.space <- 0.2
  }
  if (missing(lfc.col)) {
    lfc.col <- c("#f47d8c", "azure", "#727fb5")
  }
  if (missing(lfc.min)) {
    lfc.min <- -Inf
  }
  if (missing(lfc.max)) {
    lfc.max <- Inf
  }
  if (missing(border.size)) {
    border.size <- 0
  }
  if (missing(process.label)) {
    process.label <- 10
  }
  if (missing(limit)) {
    limit <- c(0, 0)
  }
  if (gene.order == "logFC") {
    data <- data[order(data[, Ncol], decreasing = T), , drop = F]
  }
  
  if (sum(!is.na(match(colnames(data), "logFC"))) > 0) {
    if (nlfc == 1) {
      # bug fix, do not drop
      cdata <- check_chord(data[, 1:(Ncol - 1), drop = F], limit)
      lfc <- sapply(rownames(cdata), function(x) {
        data[match(
          x,
          rownames(data)
        ), Ncol, drop = F]
      })
    } else {
      cdata <- check_chord(data[, 1:(Ncol - nlfc)], limit)
      lfc <- sapply(rownames(cdata), function(x) {
        data[
          ,
          (Ncol - nlfc + 1), drop = F
        ]
      })
    }
  } else {
    cdata <- check_chord(data, limit)
    lfc <- 0
  }
  if (missing(ribbon.col)) {
    # 36 colors
    default.col <- c("#1F78B4","#A6CEE3","#33A02C","#B2DF8A","#E31A1C","#FB9A99","#FF7F00" ,"#FDBF6F" ,"#6A3D9A", "#CAB2D6" ,"#B15928" ,"#FFFF99") 
    #default.col <- c("#EDB11A","#88ABDA","#ED665D","#AD8BC9","#FABB6E","#719F65","#ED97CA","#84A3A9","#BA5B6A","#D1D266")
    colRib <- default.col[1:dim(cdata)[2]]
  } else {
    colRib <- ribbon.col
  }
  nrib <- colSums(cdata)
  ngen <- rowSums(cdata)
  Ncol <- dim(cdata)[2]
  Nrow <- dim(cdata)[1]
  colRibb <- c()
  for (b in seq_len(length(nrib))) {
    colRibb <- c(colRibb, rep(
      colRib[b],
      202 * nrib[b]
    ))
  }
  r1 <- 1
  r2 <- r1 + 0.1
  xmax <- c()
  x <- 0
  for (r in seq_len(length(nrib))) {
    perc <- nrib[r] / sum(nrib)
    xmax <- c(xmax, (pi * perc) - space)
    if (length(x) <= Ncol - 1) {
      x <- c(x, x[r] + pi * perc)
    }
  }
  xp <- c()
  yp <- c()
  l <- 50
  for (s in 1:Ncol) {
    xh <- seq(x[s], x[s] + xmax[s], length = l)
    xp <- c(
      xp, r1 * sin(x[s]), r1 * sin(xh), r1 * sin(x[s] +
                                                   xmax[s]), r2 * sin(x[s] + xmax[s]), r2 * sin(rev(xh)),
      r2 * sin(x[s])
    )
    yp <- c(
      yp, r1 * cos(x[s]), r1 * cos(xh), r1 * cos(x[s] +
                                                   xmax[s]), r2 * cos(x[s] + xmax[s]), r2 * cos(rev(xh)),
      r2 * cos(x[s])
    )
  }
  df_process <- data.frame(x = xp, y = yp, id = rep(c(1:Ncol),
                                                    each = 4 + 2 * l
  ))
  xp <- c()
  yp <- c()
  logs <- NULL
  x2 <- seq(0 - space, -pi - (-pi / Nrow) - space, length = Nrow)
  xmax2 <- rep(-pi / Nrow + space, length = Nrow)
  for (s in 1:Nrow) {
    xh <- seq(x2[s], x2[s] + xmax2[s], length = l)
    if (nlfc <= 1) {
      xp <- c(
        xp, (r1 + 0.05) * sin(x2[s]), (r1 + 0.05) *
          sin(xh), (r1 + 0.05) * sin(x2[s] + xmax2[s]),
        r2 * sin(x2[s] + xmax2[s]), r2 * sin(rev(xh)),
        r2 * sin(x2[s])
      )
      yp <- c(
        yp, (r1 + 0.05) * cos(x2[s]), (r1 + 0.05) *
          cos(xh), (r1 + 0.05) * cos(x2[s] + xmax2[s]),
        r2 * cos(x2[s] + xmax2[s]), r2 * cos(rev(xh)),
        r2 * cos(x2[s])
      )
    } else {
      tmp <- seq(r1, r2, length = nlfc + 1)
      for (t in 1:nlfc) {
        logs <- c(logs, data[s, (dim(data)[2] + 1 - t)])
        xp <- c(
          xp, (tmp[t]) * sin(x2[s]), (tmp[t]) *
            sin(xh), (tmp[t]) * sin(x2[s] + xmax2[s]),
          tmp[t + 1] * sin(x2[s] + xmax2[s]), tmp[t +
                                                    1] * sin(rev(xh)), tmp[t + 1] * sin(x2[s])
        )
        yp <- c(
          yp, (tmp[t]) * cos(x2[s]), (tmp[t]) *
            cos(xh), (tmp[t]) * cos(x2[s] + xmax2[s]),
          tmp[t + 1] * cos(x2[s] + xmax2[s]), tmp[t +
                                                    1] * cos(rev(xh)), tmp[t + 1] * cos(x2[s])
        )
      }
    }
  }
  if (lfc[1] != 0) {
    if (nlfc == 1) {
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow),
                                                      each = 4 + 2 * l
      ), logFC = rep(lfc, each = 4 +
                       2 * l))
    } else {
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:(nlfc *
                                                             Nrow)), each = 4 + 2 * l), logFC = rep(logs,
                                                                                                    each = 4 + 2 * l
                                                             ))
    }
  } else {
    df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow),
                                                    each = 4 + 2 * l
    ))
  }
  aseq <- seq(0, 180, length = length(x2))
  angle <- c()
  for (o in aseq) {
    if ((o + 270) <= 360) {
      angle <- c(angle, o + 270)
    } else {
      angle <- c(angle, o - 90)
    }
  }
  df_texg <- data.frame(
    xgen = (r1 + gene.space) * sin(x2 +
                                     xmax2 / 2), ygen = (r1 + gene.space) * cos(x2 + xmax2 / 2),
    labels = rownames(cdata), angle = angle
  )
  df_texp <- data.frame(
    xpro = (r1 + 0.15) * sin(x + xmax / 2),
    ypro = (r1 + 0.15) * cos(x + xmax / 2),
    labels = str_wrap(colnames(cdata), 70),
    stringsAsFactors = FALSE
  )
  cols <- rep(colRib, each = 4 + 2 * l)
  x.end <- c()
  y.end <- c()
  processID <- c()
  for (gs in seq_len(length(x2))) {
    val <- seq(x2[gs], x2[gs] + xmax2[gs], length = ngen[gs] +
                 1)
    pros <- which((cdata[gs, ] != 0) == T)
    for (v in 1:(length(val) - 1)) {
      x.end <- c(x.end, sin(val[v]), sin(val[v + 1]))
      y.end <- c(y.end, cos(val[v]), cos(val[v + 1]))
      processID <- c(processID, rep(pros[v], 2))
    }
  }
  df_bezier <- data.frame(x.end = x.end, y.end = y.end, processID = processID)
  df_bezier <- df_bezier[order(df_bezier$processID, -df_bezier$y.end), ]
  x.start <- c()
  y.start <- c()
  for (rs in seq_len(length(x))) {
    val <- seq(x[rs], x[rs] + xmax[rs], length = nrib[rs] +
                 1)
    for (v in 1:(length(val) - 1)) {
      x.start <- c(x.start, sin(val[v]), sin(val[v + 1]))
      y.start <- c(y.start, cos(val[v]), cos(val[v + 1]))
    }
  }
  df_bezier$x.start <- x.start
  df_bezier$y.start <- y.start
  df_path <- bezier(df_bezier, colRib)
  if (length(df_genes$logFC) != 0) {
    tmp <- sapply(df_genes$logFC, function(x) {
      ifelse(x >
               lfc.max, lfc.max, x)
    })
    logFC <- sapply(tmp, function(x) {
      ifelse(x < lfc.min,
             lfc.min, x
      )
    })
    df_genes$logFC <- logFC
  }
  
  if (ncol(degDf) > 1) {
    bks <- unique(c(min(df_genes$logFC), max(df_genes$logFC)))
    lbs <- unique(c(round(min(df_genes$logFC)), round(max(df_genes$logFC))))
    if (length(bks) > length(lbs)) {
      lbs <- unique(
        c(round(min(df_genes$logFC), 3), round(max(df_genes$logFC), 3))
      )
    }
    g <- ggplot()+geom_polygon(data = df_genes, aes(x, y,
                                                    group = id,
                                                    fill = logFC
    ), inherit.aes = F, color = "white", size = 0.5) +
      scale_fill_gradient2("logFC",
                           space = "Lab", low = lfc.col[3],
                           mid = lfc.col[2], high = lfc.col[1], guide = guide_colorbar(
                             title.position = "top",
                             title.hjust = 0.5,title.vjust = 1,
                             order = 1
                           ), breaks = bks, labels = lbs
      ) + 
      theme(
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.box = "vertical", legend.direction = "horizontal",
        legend.spacing.y = unit(0, "line"),
        legend.spacing.x = unit(0, "line"),
        legend.key = element_rect(fill = NA, color = NA),
        legend.key.height = unit( # symbol size
          1.2,
          "line"
        ), legend.key.width = unit(
          .6,
          "line"
        )
      ) + theme(
        plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(2, 2, 2, 2), 'cm'),
        legend.box.margin=margin(50, 0, 0, 0)
      )
  } else {
    g <- ggplot()+geom_polygon(
      data = df_genes, aes(x, y, group = id),
      fill = "gray50", inherit.aes = F, color = "white", size = 0.5
    ) +
      theme(
        legend.position = "bottom",
        # legend.background = element_blank(),
        legend.box = "vertical", legend.direction = "horizontal",
        legend.spacing.y = unit(0, "pt"),
        legend.spacing.x = unit(0, "line"),
        legend.key = element_rect(fill = NA, color = NA),
        legend.key.height = unit( # symbol size
          1.2,
          "line"
        ), legend.key.width = unit(
          .6,
          "line"
        )
      ) + theme(
        plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(2, 2, 2, 2), 'cm'),
        legend.box.margin=margin(50, 0, 0, 0)
      )
  }
  g+geom_polygon(data = df_process, aes(x, y,
                                        group = id
  ), fill = "gray70", inherit.aes = F, color = NA, size = border.size) +
    geom_polygon(
      data = df_process, aes(x, y, group = id),
      fill = cols, inherit.aes = F, alpha = 0.6, color = NA,
      size = border.size
    ) +
    geom_point(aes(x = xpro, y = ypro, size = factor(labels,
                                                     levels = labels
    ), shape = NA), data = df_texp) +
    scale_size_discrete(labels=function(x)ifelse(nchar(x)>50,paste0(substr(x,1,50),"..."),x)) +
    guides(size = guide_legend("",
                               ncol = 2, byrow = T,
                               override.aes = list(
                                 shape = 22, fill = unique(cols), colour = NA,
                                 size = 10
                               )
    )) +
    theme(legend.key = element_blank(),legend.text = element_text(
      size = process.label,margin = margin(r = 0, b = 0, l = 0, t = 0, unit = "pt")
    )) +
    geom_text(aes(xgen * .95, ygen * .95, label = labels, angle = angle),
              hjust = 1,
              data = df_texg, size = gene.size
    ) +
    geom_polygon(aes(
      x = lx,
      y = ly, group = ID
    ),
    data = df_path, fill = colRibb, alpha = 0.8,
    color = NA, size = border.size, inherit.aes = F
    ) +
    coord_fixed(clip = "off") +
    labs(title = title) +
    theme_blank
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("png","pdf"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 8, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 10, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  
  # 此图参数
  parser$add_argument("-t","--type",type="character", default="K", help="enrich database class,such as G/K/W/R/I", metavar="character")
  parser$add_argument("-ic","--incompare",default = "A_B", help = "比较组名称，默认A_B")
  parser$add_argument("-if","--inputfile",default = NULL, help = "输入-diff-文件，默认不输入")
  parser$add_argument("-enf","--enrichfile",default = NULL, help = "输入-enrichment-文件，默认不输入")
  parser$add_argument("-f","--filt",default = "T", type= "character",help = "数据是否需要筛选,默认T")
  parser$add_argument("-ip","--inputpath",default = "./", help = "-diff-文件所在路径，默认当前路径")
  parser$add_argument("-s","--savepath",default = "./enrich/", help = "富集结果保存路径,默认./enrich/GO/")
  parser$add_argument("-n","--number",default = 10,help = "top number of term,默认10")
  args <- parser$parse_args()
  enrich_chordplot <- do.call(mainchord,args = args)
}