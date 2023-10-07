#'Data processing: Replace missing variables
#'@description Replace missing variables by min/mean/median/KNN/BPCA/PPCA/svdImpute.
#'@usage ImputeMissingVar(mSetObj, method)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param method Select the option to replace missing variables, either 
#'replacement based on the minimum ("min), the mean ("mean"), or the median ("median") value of each feature columns,
#'or several options to impute the missing values, using k-nearest neighbour ("KNN"), probabilistic PCA ("PPCA"), 
#'Bayesian PCA ("BPCA") method, or Singular Value Decomposition ("svdImpute") 
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@import qs
#' @export
#'
ImputeMissingVar2 <- function(mSetObj=NA, method="min",percent = 0.5){
  if(.on.public.web){
    # make this lazy load
    if(!exists("my.impute.missing")){ # public web on same user dir
      .load.scripts.on.demand("util_missing.Rc");    
    }
    return(my.impute.missing2(mSetObj, method,percent));
  }else{
    return(my.impute.missing2(mSetObj, method,percent));
  }
}

#' @export
my.impute.missing2 <- function(mSetObj=NA, method="min",percent = 0.5){
  
  mSetObj <- .get.mSet(mSetObj);
  
  int.mat <- qs::qread("preproc.qs");
  new.mat <- NULL;
  msg <- mSetObj$msgSet$replace.msg;
  
  if(method=="none"){
    new.mat<-int.mat;
  }else if(method=="exclude"){
    good.inx<-apply(is.na(int.mat), 2, sum)==0
    new.mat<-int.mat[,good.inx, drop=FALSE];
    msg <- c(msg,"Variables with missing values were excluded.");
  }else if(method=="halfmin"){
    minvalue <- min(int.mat[int.mat!=0],na.rm = T)/2;
    int.mat[is.na(int.mat)] <- minvalue;
    int.mat[int.mat==0] <- minvalue;
    new.mat<- int.mat;
    msg <- c(msg,"Missing variables were replaced by halfmin (1/2 of the min positive value for all variable)");
  }else if(method=="valuemin"){
    minvalue <- min(int.mat[int.mat!=0],na.rm = T);
    int.mat[is.na(int.mat)] <- minvalue;
    int.mat[int.mat==0] <- minvalue;
    msg <- c(msg,"Missing variables were replaced by valuemin (the min positive value for all variable)");
  }else if(method=="mean_half"){
    cls <- mSetObj$dataSet$proc.cls;
    new.mat<- ReplaceMissingBymean_half(int.mat,group = cls,precent = percent);
    msg <- c(msg,"Missing variables were replaced by mean_half_na (the min positive value for all variable)");
  }else if(method=="na_half"){
    cls <- mSetObj$dataSet$proc.cls;
    new.mat<- ReplaceMissingByna_half(int.mat,group = cls,precent = percent);
    msg <- c(msg,"Missing variables were replaced by mean_half (the min positive value for all variable)");
  }else if(method=="min"){
    new.mat<- ReplaceMissingByLoD2(int.mat);
    msg <- c(msg,"Missing variables were replaced by LoDs (1/5 of the min positive value for each variable)");
  }else if(method=="colmin"){
    new.mat<-apply(int.mat, 2, function(x){
      if(sum(is.na(x))>0){
        x[is.na(x)]<-min(x,na.rm=T)/2;
      }
      x;
    });
    msg <- c(msg,"Missing variables were replaced by 1/2 of min values for each feature column.");
  }else if (method=="mean"){
    new.mat<-apply(int.mat, 2, function(x){
      if(sum(is.na(x))>0){
        x[is.na(x)]<-mean(x,na.rm=T);
      }
      x;
    });
    msg <- c(msg,"Missing variables were replaced with the mean value for each feature column.");
  }else if (method == "median"){
    new.mat<-apply(int.mat, 2, function(x){
      if(sum(is.na(x))>0){
        x[is.na(x)]<-median(x,na.rm=T);
      }
      x;
    });
    msg <- c(msg,"Missing variables were replaced with the median for each feature column.");
  }else{
    if(method == "knn_var"){
      new.mat<-t(impute::impute.knn(t(int.mat))$data);
    }else if(method == "knn_smp"){
      new.mat<-impute::impute.knn(data.matrix(int.mat))$data;
    }else{
      if(method == "bpca"){
        new.mat<-pcaMethods::pca(int.mat, nPcs =5, method="bpca", center=T)@completeObs;
      }else if(method == "ppca"){
        new.mat<-pcaMethods::pca(int.mat, nPcs =5, method="ppca", center=T)@completeObs;
      }else if(method == "svdImpute"){
        new.mat<-pcaMethods::pca(int.mat, nPcs =5, method="svdImpute", center=T)@completeObs;
      }else{
        new.mat<-int.mat
        new.mat[is.na(new.mat)]<-0
      }
    }
    msg <- c(msg, paste("Missing variables were imputated using", toupper(gsub("_", "-", method))));
  }
  
  mSetObj$dataSet$proc.feat.num <- ncol(int.mat);
  qs::qsave(as.data.frame(new.mat), file="data_proc.qs");
  mSetObj$msgSet$replace.msg <- msg;
  return(.set.mSet(mSetObj))
}

#' @export
ReplaceMissingByLoD2 <- function(int.mat){
  int.mat <- as.matrix(int.mat);
  
  rowNms <- rownames(int.mat);
  colNms <- colnames(int.mat);
  int.mat <- apply(int.mat, 2, .replace.by.lod2);
  rownames(int.mat) <- rowNms;
  colnames(int.mat) <- colNms;
  return (int.mat);
}

.replace.by.lod2 <- function(x){
  lod <- min(x[x>0], na.rm=T)/5;
  x[x==0|is.na(x)] <- lod;
  return(x);
}

ReplaceMissingBymean_half <- function(int.mat,group,precent = 0.5) {
  int.mat <- as.matrix(int.mat);
  
  rowNms <- rownames(int.mat);
  colNms <- colnames(int.mat);
  # minvalue <- min(int.mat[int.mat != 0],na.rm = T)/2
  minvalue <- apply(int.mat,1,function(x){min(x[x != 0],na.rm = T)/2})
  int.mat <- apply(int.mat, 2, .replace.by.mean_half,group = group,minvalue = minvalue,precent = precent);
  rownames(int.mat) <- rowNms;
  colnames(int.mat) <- colNms;
  return (int.mat);
}

.replace.by.mean_half <- function(x,group,minvalue,precent = 0.5){
  
  singlegroup <- unique(group)
  
  for ( i in singlegroup) {
    if(length(x[group == i & (x==0|is.na(x))])/length(x[group == i]) <= precent){
      x[group == i & (x==0|is.na(x))] <- mean(x[group == i & (x > 0)],na.rm=T)
    }else{
      x[group == i & (x==0|is.na(x))] <- minvalue[group == i & (x==0|is.na(x))]
    }
  }
  
  return(x);
}

ReplaceMissingByna_half <- function(int.mat,group,precent = 0.5) {
  int.mat <- as.matrix(int.mat);
  
  rowNms <- rownames(int.mat);
  colNms <- colnames(int.mat);
  int.mat <- apply(int.mat, 2, .replace.by.na_half,group = group,precent = precent);
  rownames(int.mat) <- rowNms;
  colnames(int.mat) <- colNms;
  return (int.mat);
}

.replace.by.na_half <- function(x,group,precent = 0.5){
  
  singlegroup <- unique(group)
  
  for ( i in singlegroup) {
    if(length(x[group == i & (x==0|is.na(x))])/length(x[group == i]) <= precent){
      x[group == i & (x==0|is.na(x))] <- mean(x[group == i & (x > 0)],na.rm=T)
    }
  }
  
  return(x);
}
