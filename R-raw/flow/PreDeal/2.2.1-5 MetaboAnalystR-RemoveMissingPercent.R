#'Data processing: remove variables with missing values
#'@description Remove variables based upon a user-defined percentage cut-off of missing values.
#'If a user specifies a threshold of 20% (0.2), it will remove variables that are missing
#'in at least 20% of all samples.
#'@usage RemoveMissingPercent(mSetObj, percent)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param percent Input the percentage cut-off you wish to use. For instance, 50 percent is represented by percent=0.5. 
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@import qs
#' @export
#'
RemoveMissingPercent2 <- function(mSetObj=NA, percent=percent,
                                  missvarremovebygroup = T){
  
  mSetObj <- .get.mSet(mSetObj);
  if(!.on.public.web & !is.null(mSetObj$dataSet$norm)){    
    int.mat <- mSetObj$dataSet$norm;
    good.inx <- apply(is.na(int.mat), 2, sum)/nrow(int.mat)<percent;
    mSetObj$dataSet$norm <- as.data.frame(int.mat[,good.inx, drop=FALSE]);
  }else if(missvarremovebygroup){
    int.mat <- qs::qread("preproc.qs");
    cls <- mSetObj$dataSet$proc.cls;
    int.mat[int.mat == 0] <- NA;
    good.inx <- apply(int.mat, 2, mulnaratio,group = cls)<=percent;
    preproc <- as.data.frame(int.mat[,good.inx, drop=FALSE]);
    qs::qsave(preproc, "preproc.qs");
  }else{  
    int.mat <- qs::qread("preproc.qs");
    int.mat[int.mat == 0] <- NA;
    good.inx <- apply(is.na(int.mat), 2, sum)/nrow(int.mat)<=percent;
    preproc <- as.data.frame(int.mat[,good.inx, drop=FALSE]);
    qs::qsave(preproc, "preproc.qs");
  }
  mSetObj$msgSet$replace.msg <- c(mSetObj$msgSet$replace.msg, paste(sum(!good.inx), "variables were removed for threshold", round(100*percent, 2), "percent."));
  return(.set.mSet(mSetObj));
}

#' @export
mulnaratio <- function(x,group){
  
  singlegroup <- unique(group)
  rationum <- NULL
  
  for ( i in singlegroup) {
    
    data <- x[group == i]
    rationum2 <- sum(is.na(data))/length(data)
    rationum <- c(rationum,rationum2)
    
  }
  
  return(min(rationum))
  
}
