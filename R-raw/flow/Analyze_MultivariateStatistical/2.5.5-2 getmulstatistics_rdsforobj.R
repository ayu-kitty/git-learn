#!/opt/conda/bin/Rscript

#' getmulstatistics_rdsforobj
#'
#' 运行多元统计分析
#' 
#' @param obj 数据
#' @param mulstatistics_rds_savepath 结果保存路径
#'
#' @export
getmulstatistics_rdsforobj <- function(obj,
                                       mulstatistics_rds_savepath = "oecloud/mulstatisticsanalyst",
                                       all = T){
  data <- obj
  
  if(all){
    group <- unique(data[["info"]][["sample"]][["class"]][[1]])
    data$statistics$ALL$rawdif <- "ALL"
    data$statistics$ALL$difname <- "ALL"
    data$statistics$ALL$class <- group
    data$statistics$ALL$mode <- "PCA"
    data$statistics$ALL$scale <- data$info$handle$PCAscaleC
    data$statistics$ALL$log10L <- data[["info"]][["dif"]][["log10L"]][[1]][1]
    data$statistics$ALL$predI <- NA
    data$statistics$ALL$orthoI <- 0
    data$statistics$ALL$permI <- 0
    data$statistics$ALL$shape <- whto(data.frame(data$info$class$group,
                                                 data$info$class$shape,
                                                 stringsAsFactors = F),
                                      group)
    data$statistics$ALL$colour <- whto(data.frame(data$info$class$group,
                                                  data$info$class$colour,
                                                  stringsAsFactors = F),
                                       group)
    data$statistics$ALL$fill <- whto(data.frame(data$info$class$group,
                                                data$info$class$fill,
                                                stringsAsFactors = F),
                                     group)
    data$statistics$ALL$data <- readdata(paste0(mulstatistics_rds_savepath,"/mulstatistic_PCA-All.rds"))
    
    group <- group[group!="QC"]
    data$statistics$ALL2$rawdif <- "Allsample"
    data$statistics$ALL2$difname <- "Allsample"
    data$statistics$ALL2$class <- group
    data$statistics$ALL2$mode <- "PCA"
    data$statistics$ALL2$scale <- data$info$handle$PCAscaleC
    data$statistics$ALL2$log10L <- data[["info"]][["dif"]][["log10L"]][[1]][1]
    data$statistics$ALL2$predI <- NA
    data$statistics$ALL2$orthoI <- 0
    data$statistics$ALL2$permI <- 0
    data$statistics$ALL2$shape <- whto(data.frame(data$info$class$group,
                                                  data$info$class$shape,
                                                  stringsAsFactors = F),
                                       group)
    data$statistics$ALL2$colour <- whto(data.frame(data$info$class$group,
                                                   data$info$class$colour,
                                                   stringsAsFactors = F),
                                        group)
    data$statistics$ALL2$fill <- whto(data.frame(data$info$class$group,
                                                 data$info$class$fill,
                                                 stringsAsFactors = F),
                                      group)
    data$statistics$ALL2$data <- readdata(paste0(mulstatistics_rds_savepath,"/mulstatistic_PCA-Allsample.rds"))
    
    data$statistics$ALL3$rawdif <- "Allsample"
    data$statistics$ALL3$difname <- "Allsample"
    data$statistics$ALL3$class <- group
    data$statistics$ALL3$mode <- "PLS"
    data$statistics$ALL3$scale <- data[["info"]][["dif"]][["PLSscaleC"]][[1]][1]
    data$statistics$ALL3$log10L <- data[["info"]][["dif"]][["log10L"]][[1]][1]
    data$statistics$ALL3$predI <- NA
    data$statistics$ALL3$orthoI <- 0
    data$statistics$ALL3$permI <- 200
    data$statistics$ALL3$shape <- whto(data.frame(data$info$class$group,
                                                  data$info$class$shape,
                                                  stringsAsFactors = F),
                                       group)
    data$statistics$ALL3$colour <- whto(data.frame(data$info$class$group,
                                                   data$info$class$colour,
                                                   stringsAsFactors = F),
                                        group)
    data$statistics$ALL3$fill <- whto(data.frame(data$info$class$group,
                                                 data$info$class$fill,
                                                 stringsAsFactors = F),
                                      group)
    data$statistics$ALL3$data <- readdata(paste0(mulstatistics_rds_savepath,"/mulstatistic_PLS-DA-Allsample.rds"))
  }
  
  for (i in 1:length(data$info$dif$compare)) {
    
    data$statistics$PCA[[i]] <- list()
    data$statistics$PLS[[i]] <- list()
    data$statistics$OPLS[[i]] <- list()
    
    for (j in 1:length(data$info$dif$compare[[i]])) {
      data$statistics$PCA[[i]][[j]] <- list()
      data$statistics$PLS[[i]][[j]] <- list()
      data$statistics$OPLS[[i]][[j]] <- list()
      
      data$statistics$PCA[[i]][[j]]$rawdif <- data$info$dif$compare[[i]][j]
      data$statistics$PCA[[i]][[j]]$difname <- gsub("/", "_", data$info$dif$compare[[i]][j])
      data$statistics$PCA[[i]][[j]]$class <- unlist(strsplit(data$info$dif$compare[[i]][j], split = "/"))
      data$statistics$PCA[[i]][[j]]$mode <- "PCA"
      data$statistics$PCA[[i]][[j]]$scale <- data$info$dif$PCAscaleC[[i]][j]
      data$statistics$PCA[[i]][[j]]$predI <- NA
      data$statistics$PCA[[i]][[j]]$orthoI <- 0
      data$statistics$PCA[[i]][[j]]$log10L <- data$info$dif$log10L[[i]][j]
      data$statistics$PCA[[i]][[j]]$shape <- whto(data.frame(data$info$class$group,
                                                             data$info$class$shape,
                                                             stringsAsFactors = F),
                                                  data$statistics$PCA[[i]][[j]]$class)
      data$statistics$PCA[[i]][[j]]$colour <- whto(data.frame(data$info$class$group,
                                                              data$info$class$colour,
                                                              stringsAsFactors = F),
                                                   data$statistics$PCA[[i]][[j]]$class)
      data$statistics$PCA[[i]][[j]]$fill <- whto(data.frame(data$info$class$group,
                                                            data$info$class$fill,
                                                            stringsAsFactors = F),
                                                 data$statistics$PCA[[i]][[j]]$class)
      
      
      data$statistics$OPLS[[i]][[j]]$rawdif <- data$info$dif$compare[[i]][j]
      data$statistics$OPLS[[i]][[j]]$difname <- gsub("/", "_", data$info$dif$compare[[i]][j])
      data$statistics$OPLS[[i]][[j]]$class <- unlist(strsplit(data$info$dif$compare[[i]][j], split = "/"))
      if (length(data$statistics$OPLS[[i]][[j]]$class) > 2) {
        data$info$dif$OPLSpermI[[i]][j] <- 0
        data$info$dif$PLSpermI[[i]][j] <- 200
      }
      data$statistics$OPLS[[i]][[j]]$mode <- "OPLS"
      data$statistics$OPLS[[i]][[j]]$scale <- data$info$dif$OPLSscaleC[[i]][j]
      data$statistics$OPLS[[i]][[j]]$predI <- 1
      data$statistics$OPLS[[i]][[j]]$orthoI <- NA
      data$statistics$OPLS[[i]][[j]]$permI <- data$info$dif$OPLSpermI[[i]][j]
      data$statistics$OPLS[[i]][[j]]$log10L <- data$info$dif$log10L[[i]][j]
      data$statistics$OPLS[[i]][[j]]$shape <- data$statistics$PCA[[i]][[j]]$shape
      data$statistics$OPLS[[i]][[j]]$colour <- data$statistics$PCA[[i]][[j]]$colour
      data$statistics$OPLS[[i]][[j]]$fill <- data$statistics$PCA[[i]][[j]]$fill
      
      data$statistics$PLS[[i]][[j]]$rawdif <- data$info$dif$compare[[i]][j]
      data$statistics$PLS[[i]][[j]]$difname <- gsub("/", "_", data$info$dif$compare[[i]][j])
      data$statistics$PLS[[i]][[j]]$class <- unlist(strsplit(data$info$dif$compare[[i]][j], split = "/"))
      data$statistics$PLS[[i]][[j]]$mode <- "PLS"
      data$statistics$PLS[[i]][[j]]$scale <- data$info$dif$PLSscaleC[[i]][j]
      data$statistics$PLS[[i]][[j]]$predI <- NA
      data$statistics$PLS[[i]][[j]]$orthoI <- 0
      data$statistics$PLS[[i]][[j]]$permI <- data$info$dif$PLSpermI[[i]][j]
      data$statistics$PLS[[i]][[j]]$log10L <- data$info$dif$log10L[[i]][j]
      data$statistics$PLS[[i]][[j]]$shape <- data$statistics$PCA[[i]][[j]]$shape
      data$statistics$PLS[[i]][[j]]$colour <- data$statistics$PCA[[i]][[j]]$colour
      data$statistics$PLS[[i]][[j]]$fill <- data$statistics$PCA[[i]][[j]]$fill
      
      comparename <- gsub(pattern = "/",replacement = "-vs-",x = data$info$dif$compare[[i]][[j]])
      group <- unlist(strsplit(x = data$info$dif$compare[[i]][[j]],split = "/"))
      
      filename <- paste0(mulstatistics_rds_savepath,"/mulstatistic_PCA-",comparename,".rds")
      if(file.exists(filename)){
        data$statistics$PCA[[i]][[j]]$data <- readdata(filename)
      }
      
      filename <- paste0(mulstatistics_rds_savepath,"/mulstatistic_PLS-DA-",comparename,".rds")
      if(file.exists(filename)){
        data$statistics$PLS[[i]][[j]]$data <- readdata(filename)
      }
      
      if(length(group) > 2){
      }else{
        filename <- paste0(mulstatistics_rds_savepath,"/mulstatistic_OPLS-DA-",comparename,".rds")
        if(file.exists(filename)){
          data$statistics$OPLS[[i]][[j]]$data <- readdata(filename)
        }
      }
    }
  }
  
  return(data)
}