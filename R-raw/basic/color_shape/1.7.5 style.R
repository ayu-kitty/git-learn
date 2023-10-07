
# 测试
# setwd("/data/hstore1/database/test/2023-03-21多元统计分析/")
# stylefun_group()
# stylefun_group(styletype = "shape")

#' @export
stylefun_group <- function(classfile = "classtype.xlsx",
                           subset = NULL,
                           styletype = "fill",
                           n = NULL,
                           fillpalette = "procol",
                           filloutpalette = NULL,
                           colourpalette = fillpalette,
                           colouroutpalette = NULL,
                           shapepalette = "common",
                           shapeoutpalette = NULL,
                           value = "Group",
                           object = NULL,
                           sample = F,
                           samplefile = paste0(dirname(classfile),"/classfile.yaml"),
                           ...){
  classtry <- try({
    class <- readdata(classfile)
    
    if(sample){
      
      colors <- stylefun_group(classfile = classfile,
                               styletype = styletype,
                               fillpalette = fillpalette,
                               filloutpalette = filloutpalette,
                               colourpalette = colourpalette,
                               colouroutpalette = colouroutpalette,
                               shapepalette = shapepalette,
                               shapeoutpalette = shapeoutpalette,
                               sample = F)
      
      sampledata <- readdata(samplefile)
      
      class <- NULL
      
      for ( i in 1:length(sampledata)) {
        
        class2 <- data.frame("Group" = sampledata[[i]],
                             "type" = unname(colors[names(sampledata)[i]]))
        colnames(class2)[2] <- styletype
        
        class <- rbind(class,class2)
      }
      class <- class[!duplicated(class[,1]),]
    }
    
  },silent = T)
  
  if("try-error" %in% class(classtry)){
    warning(paste0("未找到",classfile,"分组样式模板文件,将自动生成样式模板"))
    class <- NULL
  }
  
  if ( !is.null(object) ){
    if ("data.frame" %in% class(object)){
      colid <- ifelse( is.null(value), colnames(object)[1], value )
      subset <- unique(object[[colid]])
    }else if ( is.factor(object) ){
      subset <- unique(object)
    }else if( is.list(object) ){
      subset <- names(object)
    }else if ( is.vector(object) ){
      subset <- unique(object)
    }
  }
  
  if(!is.null(subset)){
    subset <- as.character(subset)
  }

  
  if (styletype == "fill") {
    
    if("fill" %in%  colnames(class)){
      colors <- SelectColors(palette = fillpalette,object = class,subset = subset,n = n,outpalette = class[,"fill"],value = value)
    }else{
      colors <- SelectColors(palette = fillpalette,object = class,subset = subset,n = n,outpalette = filloutpalette,value = value)
    }
    
  }else if (styletype == "colour") {
    
    if("colour" %in%  colnames(class)){
      colors <- SelectColors(palette = colourpalette,object = class,subset = subset,n = n,outpalette = class[,"colour"],value = value)
    }else{
      colors <- SelectColors(palette = colourpalette,object = class,subset = subset,n = n,outpalette = colouroutpalette,value = value)
    }
    
  }else if (styletype == "shape") {
    
    if("shape" %in%  colnames(class)){
      colors <- SelectShape(palette = shapepalette,object = class,subset = subset,n = n,outpalette = class[,"shape"],value = value)
    }else{
      colors <- SelectShape(palette = shapepalette,object = class,subset = subset,n = n,outpalette = shapeoutpalette,value = value)
    }
    
  }else{
    
    stop("styletype风格类型返回错误，填写fill、colour、shape")
    
  }
  
  return(colors)
  
}

#' @export
stylefun_group_build <- function(classfile = "./oecloud/rawdata/classfile.yaml",
                                 subset = NULL,
                                 styletype = "fill",
                                 n = NULL,
                                 fillpalette = "procol",
                                 filloutpalette = NULL,
                                 colourpalette = "procol",
                                 colouroutpalette = NULL,
                                 shapepalette = "common",
                                 shapeoutpalette = NULL,
                                 value = "Group",
                                 object = NULL,
                                 sample = F,
                                 samplefile = paste0(dirname(classfile),"/classfile.yaml")){
  suppressMessages(library("pryr"))
  
  tmpfun <- make_function(args = c(list(classfile = classfile,
                                        subset = subset,
                                        styletype = styletype,
                                        n = n,
                                        fillpalette = fillpalette,
                                        filloutpalette = filloutpalette,
                                        colourpalette = colourpalette,
                                        colouroutpalette = colouroutpalette,
                                        shapepalette = shapepalette,
                                        shapeoutpalette = shapeoutpalette,
                                        value = value,
                                        object = object,
                                        sample = sample,
                                        samplefile = samplefile),alist(... = )), 
                          body = body(stylefun_group))
  
  if("package:lmbio" %in% search()){
    detach("package:lmbio")
    assignInNamespace(x = "stylefun_group",value = tmpfun,ns = asNamespace("lmbio"))
    suppressMessages(library("lmbio"))
  }else{
    assignInNamespace(x = "stylefun_group",value = tmpfun,ns = asNamespace("lmbio"))
  }
  
}
