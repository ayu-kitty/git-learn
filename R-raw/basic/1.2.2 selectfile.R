
#' 根据重要性、时间、测试等进行文件选取
#'
#' @param path 相对数据库查找路径
#' @param database 数据库路径
#' @param file 查找文件
#' @param selectmode 选择模式，默认!(重要性)，可填写"!"(重要性)、"#"(测试)、""(时间),"目录名"
#'
#' @export
selectfile <- function(path = "database/kegg",
                       database = lmbio::databasepath(path = path),
                       file = "hsa",
                       selectmode = "!",
                       version = lmbio::lmbioversion){
  
  filepath <- list.files(path = database,pattern = paste0("^",file,"$","|!",file,"$|^#",file,"$"),
                         recursive = T,full.names = T,include.dirs = T)
  
  filepath2 <- list.files(path = database,pattern = paste0("^",file,"$","|!",file,"$|^#",file,"$"),
                         recursive = T,full.names = F,include.dirs = T)
  
  filepath3 <- gsub(pattern = "/.*",replacement = "",filepath2)
  filepath3 <- gsub(pattern = "^!",replacement = "",filepath3)
  filepath3 <- gsub(pattern = "^#",replacement = "",filepath3)
  
  if(any(grepl(pattern = paste0("~",version,"$"),x = filepath3))){
    filepath <- filepath[grepl(pattern = paste0("~",version,"$"),x = filepath3)]
    filepath2 <- filepath2[grepl(pattern = paste0("~",version,"$"),x = filepath3)]
  }else{
    warning("未找到version:",version,"相关版本，将使用所有版本检索",immediate. = T)
  }
  
  if(length(filepath) > 0){
    cat(paste(c("查找到以下文件:",filepath),collapse = "\n"),"\n")
  }else{
   stop(paste0("在",database,"目录下未找到",file)) 
  }
  
  if(selectmode == "!"){
    
    cat("使用!进行重要性筛选\n")
    
    if(any(grepl(pattern = paste0("/!",file,"$"),x = filepath))){
      
      if(any(grepl(pattern = paste0(database,"/!"),x = filepath))){
       warning("目录及文件均查询到!,可能出错请注意",immediate. = T) 
      }
      
      filepath <- filepath[grepl(pattern = paste0("/!",file,"$"),x = filepath)]
      
    }else if(any(grepl(pattern = paste0(database,"/!"),x = filepath))){
      
      filepath <- filepath[grepl(pattern = paste0(database,"/!"),x = filepath)]
      
    }else{
      
      print("未查询到!，将参数selectmode修改为空")
      
      filepath <- selectfile(database = database,
                             file = file,
                             selectmode = "",
                             version = version)
      
      return(filepath)
    }
    
  }else if(selectmode == "#"){
    
    cat("使用#进行测试筛选\n")
    
    if(any(grepl(pattern = paste0("/#",file,"$"),x = filepath))){
      
      if(any(grepl(pattern = paste0(database,"/#"),x = filepath))){
        warning("目录及文件均查询到#,可能出错请注意",immediate. = T) 
      }
      
      filepath <- filepath[grepl(pattern = paste0("/#",file,"$"),x = filepath)]
      
    }else if(any(grepl(pattern = paste0(database,"/#"),x = filepath))){
      
      filepath <- filepath[grepl(pattern = paste0(database,"/#"),x = filepath)]
      
    }else{
      
      stop(paste0("在",database,"目录下未找到#模式的",file)) 
      
    }
    
  }else if(selectmode == ""){
    
    cat("除#外使用时间进行筛选\n")
    
    filepath3 <- gsub(pattern = "/.*",replacement = "",filepath2)
    filepath3 <- gsub(pattern = "^!",replacement = "",filepath3)
    filepath3 <- gsub(pattern = "^#.*",replacement = "0",filepath3)
    
    filepath <- filepath[filepath3 == max(filepath3) & filepath3 !="0"]
    
    if(length(filepath) == 0){
      
      stop(paste0("在",database,"目录下未找到",file)) 
      
    }
    
  }else{
    
    cat(paste0("使用",selectmode,"进行筛选\n"))
    
    if(any(grepl(pattern = paste0(database,"/",selectmode,"/"),x = filepath))){
      
      filepath <- filepath[grepl(pattern = paste0(database,"/",selectmode,"/"),x = filepath)]
      
    }else{
      stop(paste0("在",database,"目录下未找到",selectmode,"的",file)) 
    }
    
  }
  
  if(length(filepath) > 1){
    cat(paste(c("####查找到多个文件，出错中####:",filepath),collapse = "\n"))
    cat("\n")
    stop("请处理后在运行")
  }
  
  cat("最终选择文件:\n")
  cat(filepath,"\n")
  return(filepath)
}
