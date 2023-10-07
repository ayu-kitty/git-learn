#!/opt/conda/bin/Rscript

#' 文件统一复制
#'
#' @param srcpath 图片、文件所在的源目录
#' @param dstpath 默认当前目录下的src/
#'
#' @export
srccopy <- function(srcpath,dstpath){
  
  if(dir.exists(file.path(srcpath))){
    if(!dir.exists(file.path(dstpath))){dir.create(dstpath)}

    # copy源图片文件
    src_file = list.files(file.path(srcpath))
    file.copy(from = file.path(srcpath, src_file),
              to = file.path(dstpath, src_file))
  
  }
}

#' 复制目录
#'
#' @param srcpath 图片、文件所在的源目录
#' @param dstpath 默认当前目录下的src/
#'
#' @export
copydir <- function(from,
                    to,
                    rmfile = NULL){
  
  if(file.exists(paths = from)){
    if(dir.exists(paths = from)){
      if(length(list.files(from))>0){
        if(!dir.exists(paths = to)){
          warning(paste0("未找到",to,"目录，自动创建此目录"),immediate. = T)
          createdir(filename = to) 
        }
        system(paste0("cp -vrf ",from,"/* ",to),
               ignore.stdout = T, ignore.stderr = T)
        system(paste0("chmod -R 777 ",to),
               ignore.stdout = T, ignore.stderr = T)
        if(!is.null(rmfile)){
          unlink(dir(path = to,pattern = rmfile,full.names = T,recursive = T))
        }
      }else{
        warning(paste0(from,"目录下无文件，不进行文件复制"),immediate. = T)
      }
    }else{
      if(!dir.exists(paths = to)){
        warning(paste0("未找到",to,"目录，自动创建此目录"),immediate. = T)
        createdir(filename = to) 
      }
      system(paste0("cp -vrf ",from," ",to),
             ignore.stdout = T, ignore.stderr = T)
      system(paste0("chmod -R 777 ",to),
             ignore.stdout = T, ignore.stderr = T)
    }
  }else{
    warning(paste0("未找到",from,"目录，跳过复制目录操作"),immediate. = T)
  }
  
}
