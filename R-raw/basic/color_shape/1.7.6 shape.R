#' @export
SelectShape <- function(palette = "common",
                        n = NULL,
                        subset = NULL,
                        object = NULL,
                        value = "Group",
                        outpalette = NULL,
                        ...){
  if ( !is.null(object) ){
    if ("data.frame" %in% class(object)){
      colid <- ifelse( is.null(value), colnames(object)[1], value )
      lnames <- unique(object[[colid]])
    }else if ( is.factor(object) ){
      lnames <- unique(object)
    }else if( is.list(object) ){
      lnames <- names(object)
    }else if ( is.vector(object) ){
      lnames <- unique(object)
    }
    n <- length(lnames)
  }else if ( !is.null(n) ) {
    lnames <- NULL
    if (!is.null(subset)){ 
      n <- length(subset)
      lnames <- subset
    }
  }else{
    lnames <- NULL
    if (!is.null(subset)){ 
      n <- length(subset)
      lnames <- subset
    }
  }
  
  if(is.null(outpalette)){
    palettes <- NULL
    for ( i in 1:length(palette)) {
      if ( !is.null(palette_shape[[palette[i]]]) ){
        palettes <- c(palettes,palette_shape[[palette[i]]])
      }else{
        stop("NO specified discrete palette FOUND!")
      }
    }
  }else{
    palettes <- outpalette
  }
  
  if ( is.null(n) ){
    shape2pick <- palettes
  }else{
    if(length(palettes) >= n){
      shape2pick <- palettes[1:n]
    }else{
      shape2pick <- rep(palettes,1000)[1:n]
    }
  }
  
  shape_use <- unname(shape2pick)
  
  if (!is.null(lnames) ){ 
    names(shape_use) <- lnames
    if (!is.null(subset)){
      if(is.character(subset)){
        if(all(subset %in% lnames)){
          shape_use <- shape_use[subset]
        }else{
          stop(paste0("以下",paste(subset[!(subset %in% lnames)],collapse = ";"),"分组不存在")) 
        }
      }else{
        shape_use <- shape_use[subset]
      }
    }
  }
  
  return(shape_use)
}

#' @export
palette_shape <- list(common = c(21, 22, 24, 23, 25))
