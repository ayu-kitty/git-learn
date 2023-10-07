#!/opt/conda/bin/Rscript

#' ggplot图片保存
#'
#' @param plot ggplot数据
#' @param savepath 保存路径
#' @param mapname 图片名称
#' @param imagetype 图片格式
#' @param other 其他参数
#' @param height 图片长度
#' @param width 图片宽度
#' @param dpi 图片分辨率
#' @param compression tiff格式压缩模式
#' @param family 字体
#'
#' @export
ggplotsave <- function(plot,
                       savepath = "./",
                       mapname = "test",
                       imagetype = c("png", "pdf"),
                       other = ggplot2::theme(),
                       height = 6,
                       width = 5,
                       compression = "zip",
                       dpi = 300,
                       family = "sans",
                       bg = "white",
                       classfile = "",
                       ...){
  suppressMessages(library("ggplot2"))
  args <- list(plot = plot,
               savepath = savepath,
               mapname = mapname,
               imagetype = imagetype,
               other = theme(),
               height = height,
               width = width,
               compression = compression,
               dpi = dpi,
               family = family)
  newargs <- list(...)
  args <- c(args,newargs)
  
  if("ggplot" %in% class(plot)){
    plot <- plot+
      theme(text = element_text(family = family),
            title = element_text(family = family))
    
    if ("list" %in% class(other)) {
      for (i in 1:length(other)) {
        if("gg" %in% class(other[[i]])){
          plot <- plot + other[[i]]
        }
      }
    } else {
      if("gg" %in% class(other)){
        plot <- plot + other
      }
    }
  }
  
  args$plot <- plot
  
  for (type in imagetype) {
    filename <- paste0(savepath,"/",mapname,".",type)
    if (is.na(type)) {
      if("ggplot" %in% class(plot)){
        print(args$plot)
      }else{
        grid::grid.draw(args$plot) 
      }
      # return(args)
      return()
    } else {
      createdir(filename = savepath)
      # print(paste0(filename,"保存中"))
      if (type == "tiff") {
        ggsave(filename = filename, 
               plot = plot, 
               dpi = dpi, 
               width = width, height = height,
               bg = bg,
               compression = compression,
               device = tiff,
               family = family)
      }else if (type == "png") {
        ggsave(filename = filename, 
               plot = plot, 
               dpi = dpi, 
               width = width, height = height,
               bg = bg,
               device = png,
               family = family)
      }else if (type == "jpg") {
        ggsave(filename = filename, 
               plot = plot, 
               dpi = dpi, 
               width = width, height = height,
               bg = bg,
               device = jpeg,
               family = family)
      }else if(type == "html"){
        suppressMessages(library("plotly"))
        htmlplot <- ggplotly(p = plot,
                             width = width * 100,
                             height = height * 100 )
        # htmltools::save_html(html = htmlplot, 
        #                      file = filename)
        # htmlwidgets::saveWidget(widget = htmlplot,
        #                         file = filename, 
        #                         selfcontained = TRUE)
        logfile <- file(tempfile(), open = "wt")
        sink(file = logfile)
        sink(file = logfile,type = "message")
        runinpath(path = savepath,
                  moudle = htmlwidgets::saveWidget,
                  moudlename = "html保存", 
                  widget = htmlplot,
                  file = paste0(mapname,".html"),
                  selfcontained = TRUE)
        sink(type = "message")
        sink()
      } else {
        ggsave(filename = filename, 
               plot = plot,  
               dpi = dpi, 
               width = width, height = height,
               bg = bg,
               family = family)
      }
      Sys.chmod(paths = filename,mode = "0777",use_umask = F)
    }
  }
  
  # return(args)
  return()
}
