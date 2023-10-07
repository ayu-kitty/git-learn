#'Normalization
#'@description This function performs row-wise normalization, transformation, and 
#'scaling of your metabolomic data. 
#'@usage Normalization(mSetObj, rowNorm, transNorm, scaleNorm, ref=NULL, ratio=FALSE, ratioNum=20)
#'@param mSetObj Input the name of the created mSetObj (see InitDataObjects)
#'@param rowNorm Select the option for row-wise normalization, "QuantileNorm" for Quantile Normalization, 
#'"CompNorm" for Normalization by a reference feature,
#'"SumNorm" for Normalization to constant sum, 
#'"MedianNorm" for Normalization to sample median, and 
#'"SpecNorm" for Normalization by a sample-specific factor.
#'@param transNorm Select option to transform the data, "LogNorm" for Log Normalization,
#'and "CrNorm" for Cubic Root Transformation. 
#'@param scaleNorm Select option for scaling the data, "MeanCenter" for Mean Centering,
#'"AutoNorm" for Autoscaling, "ParetoNorm" for Pareto Scaling, amd "RangeNorm" for Range Scaling.
#'@param ref Input the name of the reference sample or the reference feature, use " " around the name.  
#'@param ratio This option is only for biomarker analysis.
#'@param ratioNum Relevant only for biomarker analysis.  
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}, Jasmine Chong
#'McGill University, Canada
#'@import qs
#' @export
#'
Normalization2 <- function(mSetObj=NA, rowNorm, transNorm, scaleNorm, ref=NULL, ratio=FALSE, ratioNum=20){
  
  mSetObj <- .get.mSet(mSetObj);
  
  # PreparePrenormData() called already
  data <- qs::qread("prenorm.qs");
  
  print(dim(data));
  
  cls <- mSetObj$dataSet$prenorm.cls;
  
  # note, setup time factor
  if(substring(mSetObj$dataSet$format,4,5)=="mf"){
    if(is.null(mSetObj$dataSet$prenorm.facA)){
      nfacA <- mSetObj$dataSet$facA;
      nfacB <- mSetObj$dataSet$facB;
    }else{
      nfacA <- mSetObj$dataSet$prenorm.facA;
      nfacB <- mSetObj$dataSet$prenorm.facB;
    }
    
    mSetObj$dataSet$facA <- nfacA;
    mSetObj$dataSet$facB <- nfacB;
    if(mSetObj$dataSet$design.type =="time" | mSetObj$dataSet$design.type =="time0"){
      # determine time factor and should order first by subject then by each time points
      if(tolower(mSetObj$dataSet$facA.lbl) == "time"){ 
        time.fac <- nfacA;
        exp.fac <- nfacB;
      }else{
        time.fac <- nfacB;
        exp.fac <- nfacA;
      }
      # now make sure time fac is ordered
      lvls <- levels(time.fac);
      time.points <- as.numeric(as.character(lvls));
      ord.lvls <- lvls[order(time.points)];
      time.fac <- ordered(time.fac, levels = ord.lvls);
      mSetObj$dataSet$time.fac <- time.fac;
      mSetObj$dataSet$exp.fac <- exp.fac;
    }
  }
  
  colNames <- colnames(data);
  rowNames <- rownames(data);
  
  # row-wise normalization
  if(rowNorm=="QuantileNorm"){
    data<-QuantileNormalize(data);
    # this can introduce constant variables if a variable is 
    # at the same rank across all samples (replaced by its average across all)
    
    varCol <- apply(data, 2, var, na.rm=T);
    constCol <- (varCol == 0 | is.na(varCol));
    constNum <- sum(constCol, na.rm=T);
    if(constNum > 0){
      print(paste("After quantile normalization", constNum, "features with a constant value were found and deleted."));
      data <- data[,!constCol, drop=FALSE];
      colNames <- colnames(data);
      rowNames <- rownames(data);
    }
    rownm<-"Quantile Normalization";
  }else if(rowNorm=="GroupPQN"){
    grp.inx <- cls == ref;
    ref.smpl <- apply(data[grp.inx, , drop=FALSE], 2, mean);
    data<-t(apply(data, 1, ProbNorm, ref.smpl));
    rownm<-"Probabilistic Quotient Normalization by a reference group";
  }else if(rowNorm=="SamplePQN"){
    ref.smpl <- data[ref, , drop=FALSE];
    data<-t(apply(data, 1, ProbNorm, ref.smpl));
    rownm<-"Probabilistic Quotient Normalization by a reference sample";
  }else if(rowNorm=="CompNorm"){
    data<-t(apply(data, 1, CompNorm, ref));
    rownm<-"Normalization by a reference feature";
  }else if(rowNorm=="SumNorm"){
    data<-t(apply(data, 1, SumNorm));
    rownm<-"Normalization to constant sum";
  }else if(rowNorm=="MedianNorm"){
    data<-t(apply(data, 1, MedianNorm));
    rownm<-"Normalization to sample median";
  }else if(rowNorm=="MedianMeanNorm"){
    mean<-mean(apply(data, 1 ,median,na.rm = T))
    data<-t(apply(data, 1, MedianMeanNorm,mean=mean));
    rownm<-"Normalization to sample median * mean";
  }else if(rowNorm=="SpecNorm"){
    if(!exists("norm.vec")){
      norm.vec <- rep(1,nrow(data)); # default all same weight vec to prevent error
      print("No sample specific information were given, all set to 1.0");
    }
    rownm<-"Normalization by sample-specific factor";
    data<-data/norm.vec;
  }else{
    # nothing to do
    rownm<-"N/A";
  }

  # use apply will lose dimension info (i.e. row names and colnames)
  rownames(data)<-rowNames;
  colnames(data)<-colNames;
  
  # if the reference by feature, the feature column should be removed, since it is all 1
  if(rowNorm=="CompNorm" && !is.null(ref)){
    inx<-match(ref, colnames(data));
    data<-data[,-inx, drop=FALSE];
    colNames <- colNames[-inx];
  }
  
  # record row-normed data for fold change analysis (b/c not applicable for mean-centered data)
  row.norm <- as.data.frame(CleanData(data, T, T)); #moved below ratio 
  qs::qsave(row.norm, file="row_norm.qs");
  # this is for biomarker analysis only (for compound concentration data)
  if(ratio){
    min.val <- min(abs(data[data!=0]))/2;
    norm.data <- log2((data + sqrt(data^2 + min.val))/2);
    transnm<-"Log2 Normalization";
    ratio.mat <- CalculatePairwiseDiff(norm.data);
    
    fstats <- Get.Fstat(ratio.mat, cls);
    hit.inx <- rank(-fstats) < ratioNum;  # get top n
    
    ratio.mat <- ratio.mat[, hit.inx, drop=FALSE];
    
    data <- cbind(norm.data, ratio.mat);
    
    colNames <- colnames(data);
    rowNames <- rownames(data);
    mSetObj$dataSet$use.ratio <- TRUE;
    mSetObj$dataSet$proc.ratio <- data;
    
  }else{
    mSetObj$dataSet$use.ratio <- FALSE;
    # transformation
    # may not be able to deal with 0 or negative values
    if(transNorm=='LogNorm'){
      if(min(data,na.rm = T) <= 0){
        stop("transNorm的LogNorm方法不允许矩阵中有负值")
      }
      min.val <- min(abs(data[data!=0]),na.rm = T)/10;
      data<-apply(data, 2, LogNorm, min.val);
      transnm<-"Log10 Normalization";
    }else if(transNorm=='Log2minNorm'){
      if(min(data,na.rm = T) <= 0){
        stop("transNorm的Log2minNorm方法不允许矩阵中有负值")
      }
      min.val <- -log2(min(data,na.rm = T))+10
      data<-apply(data, 2, Log2minNorm, min.val);
      transnm<-"Log2min Transformation";
    }else if(transNorm=='SrNorm'){
      min.val <- min(abs(data[data!=0]),na.rm = T)/10;
      data<-apply(data, 2, SquareRootNorm, min.val);
      transnm<-"Square Root Transformation";
    }else if(transNorm=='CrNorm'){
      norm.data <- abs(data)^(1/3);
      norm.data[data<0] <- - norm.data[data<0];
      data <- norm.data;
      transnm<-"Cubic Root Transformation";
    }else{
      transnm<-"N/A";
    }
  }
  
  # scaling
  if(scaleNorm=='MeanCenter'){
    data<-apply(data, 2, MeanCenter);
    scalenm<-"Mean Centering";
  }else if(scaleNorm=='AutoNorm'){
    data<-apply(data, 2, AutoNorm);
    scalenm<-"Autoscaling";
  }else if(scaleNorm=='ParetoNorm'){
    data<-apply(data, 2, ParetoNorm);
    scalenm<-"Pareto Scaling";
  }else if(scaleNorm=='RangeNorm'){
    data<-apply(data, 2, RangeNorm);
    scalenm<-"Range Scaling";
  }else{
    scalenm<-"N/A";
  }
  
  # note after using "apply" function, all the attribute lost, need to add back
  rownames(data)<-rowNames;
  colnames(data)<-colNames;
  
  # need to do some sanity check, for log there may be Inf values introduced
  data <- CleanData(data, F, F);
  
  if(ratio){
    mSetObj$dataSet$ratio <- CleanData(ratio.mat, T, F)
  }
  
  mSetObj$dataSet$norm <- as.data.frame(data);
  if(substring(mSetObj$dataSet$format,4,5)=="mf"){
    if(!identical(rownames(mSetObj$dataSet$norm), rownames(mSetObj$dataSet$meta.info))){
      mSetObj$dataSet$meta.info <- mSetObj$dataSet$meta.info[match(rownames(mSetObj$dataSet$norm), rownames(mSetObj$dataSet$meta.info)), ]
      print("Metadata and data norm order synchronized.")
    }
    mSetObj$dataSet$meta.info <- mSetObj$dataSet$meta.info[rownames(data),]  
  }
  
  qs::qsave(mSetObj$dataSet$norm, file="complete_norm.qs");
  mSetObj$dataSet$cls <- cls;
  
  mSetObj$dataSet$rownorm.method <- rownm;
  mSetObj$dataSet$trans.method <- transnm;
  mSetObj$dataSet$scale.method <- scalenm;
  mSetObj$dataSet$norm.all <- NULL; # this is only for biomarker ROC analysis
  
  return(.set.mSet(mSetObj));
}

#'Row-wise Normalization
#'@description Row-wise norm methods, when x is a row.
#'Normalize by a sum of each sample, assume constant sum (1000).
# Return: normalized data.
#'Options for normalize by sum median, reference sample,
#'reference reference (compound), or quantile normalization
#'@param x Input data to normalize
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'
SumNorm<-function(x){
  1000*x/sum(x, na.rm=T);
}

# normalize by median
MedianNorm<-function(x){
  x/median(x, na.rm=T);
}

# normalize by median*mean
MedianMeanNorm<-function(x,mean){
  x/median(x, na.rm=T)*mean;
}

# normalize by a reference sample (probability quotient normalization)
# ref should be the name of the reference sample
ProbNorm<-function(x, ref.smpl){
  x/median(as.numeric(x/ref.smpl), na.rm=T)
}

# normalize by a reference reference (i.e. creatinine)
# ref should be the name of the cmpd
CompNorm<-function(x, ref){
  1000*x/x[ref];
}

# perform quantile normalization on the raw data (can be log transformed later by user)
# https://stat.ethz.ch/pipermail/bioconductor/2005-April/008348.html
QuantileNormalize <- function(data){
  return(t(preprocessCore::normalize.quantiles(t(data), copy=FALSE)));
}

#'Column-wise Normalization
#'@description Column-wise norm methods, when x is a column
#'Options for log, zero mean and unit variance, and
#'several zero mean and variance/SE 
#'@param x Input data
#'@param min.val Input minimum value
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada

# generalize log, tolerant to 0 and negative values
LogNorm<-function(x, min.val){
  log10((x + sqrt(x^2 + min.val^2))/2)
}

Log2minNorm<-function(x, min.val){
  log2(x)+min.val
}

# square root, tolerant to negative values
SquareRootNorm<-function(x, min.val){
  ((x + sqrt(x^2 + min.val^2))/2)^(1/2);
}

# normalize to zero mean and unit variance
AutoNorm<-function(x){
  (x - mean(x))/sd(x, na.rm=T);
}

# normalize to zero mean but variance/SE
ParetoNorm<-function(x){
  (x - mean(x))/sqrt(sd(x, na.rm=T));
}

# normalize to zero mean but variance/SE
MeanCenter<-function(x){
  x - mean(x);
}

# normalize to zero mean but variance/SE
RangeNorm<-function(x){
  if(max(x) == min(x)){
    x;
  }else{
    (x - mean(x))/(max(x)-min(x));
  }
}

#' @export
CleanData <-function(bdata, removeNA=T, removeNeg=T, removeConst=T){
  
  if(sum(bdata==Inf, na.rm=TRUE)>0){
    inx <- bdata == Inf;
    bdata[inx] <- NA;
    bdata[inx] <- max(bdata, na.rm=T)*2
  }
  if(sum(bdata==-Inf, na.rm=TRUE)>0){
    inx <- bdata == -Inf;
    bdata[inx] <- NA;
    bdata[inx] <- min(bdata, na.rm=T)/2
  }
  if(removeNA){
    if(sum(is.na(bdata))>0){
      bdata[is.na(bdata)] <- min(bdata, na.rm=T)/2
    }
  }
  if(removeNeg){
    if(sum(as.numeric(bdata<=0)) > 0){
      inx <- bdata <= 0;
      bdata[inx] <- NA;
      bdata[inx] <- min(bdata, na.rm=T)/2
    }
  }
  if(removeConst){
    varCol <- apply(data.frame(bdata), 2, var, na.rm=T); # getting an error of dim(X) must have a positive length, fixed by data.frame
    constCol <- (varCol == 0 | is.na(varCol));
    constNum <- sum(constCol, na.rm=T);
    if(constNum > 0){
      bdata <- data.frame(bdata[,!constCol, drop=FALSE], check.names = F); # got an error of incorrect number of dimensions, added drop=FALSE to avoid vector conversion
    }
  }
  bdata;
}
