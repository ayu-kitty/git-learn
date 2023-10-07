#!/opt/conda/bin/Rscript

#' 读取fcs文件（或修改）并保存为rds
#' 
#' @param originpath 原始存放fcs文件的路径
#' @param newpath 存放修改后的fcs文件的路径
#' @param saveflowdata 逻辑，是否保存flowSet对象为rds文件
#' @param RDSname 保存的rds文件名
#' 
#' @export
read_fcs <- function(originpath = "./rawdata",
                     newpath = "ordered_fcs",
                     rdsname = "flowSet_object.rds",
                     ...){
  suppressMessages(library("flowCore"))
  # 单个fcs文件读取
  
  if(length(list.files(originpath, '*\\.fcs$')) == 1){
    fs1 = list.files(originpath,'*\\.fcs$' )
    samp <- read.flowSet(files = fs1 ,path = originpath)
  }else{
    fs1 = list.files(originpath,'*\\.fcs$' )
    # 多个fcs文件分别读取 存入一个list
    cytof_list <- lapply(fs1, function(x){read.flowSet(files = x ,path = originpath)}) 
    # fcs文件数
    sample_number <- length(list.files(originpath, '*\\.fcs$'))
    # fcs文件的table(colnames())
    ids <- table(unlist(lapply(cytof_list,colnames)))
    # 列名
    id <- names(ids[ids==sample_number])
    # 新建一个目录
    createdir(filename = paste0(originpath,"/",newpath))
    # 修改fcs文件并保存到新目录中
    lapply(fs1, function(x){
      tmp = read.FCS(file.path(originpath, x))
      tmp@exprs = tmp@exprs[, id]
      tmp@parameters@data = tmp@parameters@data[match(id, tmp@parameters@data$name ), ]
      write.FCS(tmp, file.path(originpath, newpath, x))
    })
    # 读取修改后的fcs文件
    fs2 <- list.files(file.path(originpath, newpath), '*\\.fcs$' )
    samp <- read.flowSet(files = fs2, path = file.path(originpath, newpath))
  }
  
  saverds(data = samp,filename = file.path(originpath,rdsname))

  return(samp)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  
  # 基本参数
  parser$add_argument("-op","--originpath",default = "./rawdata/", 
                      help = "fcs文件所在路径,件以`/DLM20225678/rawdata/`形式传参")
  parser$add_argument("-np","--newpath",default = "ordered_fcs", 
                      help = "存放修改后的fcs文件的路径")
  parser$add_argument("-sn","--rdsname", default = "flowSet_object.rds", 
                      help = "保存的rds文件名")

  args <- parser$parse_args()
  
  readfcsargs <- do.call(what = read_fcs, args = args)
}
## 注意：每个目录和文件必须严格存在！rawdata是固定的下机数据存放目录 flowSet_object.rds是固定的文件名