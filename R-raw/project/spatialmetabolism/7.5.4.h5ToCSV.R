# counts.csv

# setwd("D:/段光前-研发-202306/4.空间转录组和空间代谢组联合/5.项目/1. DZLM2023061613-20230703/spaceranger/S2209151T_2209011T/outs/") #S2209151T_2209011T S2206301T_2207211P

#' @export
h5ToCSV <- function(path=NULL, filename="filtered_feature_bc_matrix.h5", 
                    samplename, savename="counts.csv",
                    ...){
  #BiocManager::install("hdf5r")
  library(hdf5r)
  library(Seurat) 
  library(vroom)
  library(magrittr)
  library(rlang)
  
  path <- path %||% getwd()
  file <- file.path(path, filename)
  data_sample <- Read10X_h5(file) 
  # data_seurat <- CreateSeuratObject(data_sample,project = "data_sample")
  expr_df <- data_sample %>% as.matrix() %>% data.frame()
  colnames(expr_df) <- paste0(samplename, ".", gsub(pattern = "\\.[0-9]", replacement = "", colnames(expr_df)))
  expr_df %<>% cbind(data.frame(Gene=rownames(.)), .)
  
  # save
  vroom::vroom_write(expr_df, file = file.path(path, savename), delim = ",", ...)
  
}




