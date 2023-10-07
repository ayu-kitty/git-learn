#!/opt/conda/bin/Rscript

#' GetAnalystInfo
#'
#' 获取线上登记单
#'
#' @param analysis_id 项目分析编号
#' @param token token
#'
#' @export
GetAnalystInfo <- function(analysis_id,
                           savepath =ifelse(grepl(pattern = "^http",x = analysis_id),"./",analysis_id),
                           filename = "分析确认单.xlsx",
                           token = "oebiotech.confrim.T34u",
                           overwrite = F) {
  
  if(grepl(pattern = "^http",x = analysis_id)){
    
    runinpath(path = savepath,
              moudle = download.file,
              url = analysis_id,
              destfile = filename)
    
  }else{
    print(paste0("分析编号：", analysis_id, " 线上登记单读取中"))
    suppressMessages(library("httr"))
    suppressMessages(library("jsonlite"))
    suppressMessages(library("openxlsx"))
    
    # analysis_id <- "LM2021-19113-aa"
    # token <- "oebiotech.confrim.T34u"
    
    if(analysis_id == ""){
      stop("分析编号为空")
    }
    
    PostData <- POST( url = "https://cloud.oebiotech.com/oeconfirm/get_record",
                      body = list(analysis_id = analysis_id, token = token))
    
    if (PostData$status_code == 200) {
      Postcontent <- httr::content(PostData, "raw")
      JsonData <- rawToChar(Postcontent)
      JsonData <- gsub(pattern = "NaN", replacement = '\"NaN\"', x = JsonData)
      alldata <- fromJSON(JsonData)
      
      if(length(alldata) == 1){
        print(alldata[[1]])
        return()
      }
      
      wb <- openxlsx::createWorkbook()
      
      for (i in 1:length(alldata)) {
        if (names(alldata)[i] == "complete") {
          break
        }
        
        listdata <- alldata[[names(alldata)[i]]]
        
        if (is.null(unlist(listdata))) {
          listdata <- data.frame(matrix(names(listdata), nrow = 1, byrow = TRUE))
          colnames(listdata) <- t(listdata)[, 1, drop = T]
          listdata <- listdata[-1, ]
        } else {
          for (j in 1:length(listdata)) {
            for (k in 1:length(listdata[[j]])) {
              if (is.null(listdata[[j]][[k]])) {
                listdata[[j]][[k]] <- NA
              } else if (is.na(listdata[[j]][[k]])) {
              } else if (listdata[[j]][[k]] == "" | listdata[[j]][[k]] == "NaN") {
                listdata[[j]][[k]] <- NA
              }
            }
          }
          
          listdata <- t(data.frame(matrix(unlist(listdata), nrow = length(listdata), byrow = TRUE)))
          listdata <- as.data.frame(listdata, stringsAsFactors = F)
          colnames(listdata) <- names(alldata[[names(alldata)[i]]])
          listdata[listdata == "NaN"] <- NA
          
          if (names(alldata)[i] == "样本实验信息") {
            # colnames(listdata) <- t(listdata)[, 3, drop = T]
            # listdata <- listdata[-1:-3, ]
          }
          
          if (names(alldata)[i] == "分析基本信息") {
            listdata["add", ] <- c("分析编号", analysis_id)
          }
        }
        
        # print(listdata)
        
        wb <- addsheet1(data = listdata, wb = wb, sheet = names(alldata)[i])
      }
      
      try({
        savewb(wb = wb,filename = paste0(savepath,"/",filename),overwrite = overwrite)
      })
      
    } else {
      stop("未查询到此分析编号对应分析确认单", call. = F)
    }
  }
  
  print("complete!")
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  
  suppressMessages(library("lmbio"))
  
  parser <- ArgumentParser()
  
  parser$add_argument("-ad","--analysis_id", help = "分析编号",required = T)
  parser$add_argument("-sp","--savepath",default = "", help = "保存路径")
  parser$add_argument("-fn","--filename",default = "分析确认单.xlsx",help = "保存文件名")
  parser$add_argument("-ow","--overwrite",default = F,help = "是否覆盖原始分析单",action='store_true')
  
  args <- parser$parse_args()
  
  if(args$savepath == ""){args$savepath <- NULL}
  
  result <- do.call(what = GetAnalystInfo,args = args) 
  
}
