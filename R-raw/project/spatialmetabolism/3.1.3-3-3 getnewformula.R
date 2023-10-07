#!/opt/conda/bin/Rscript

#' @export
getnewformula <- function(formula = "C2H2O",
                          adductElement = "-H2+Na",
                          nM = 1,
                          ...){
  element <- stringr::str_extract_all(string = formula,pattern = "[A-Z][a-z]*[0-9]*")[[1]]
  data <- data.frame(rawelement = element)
  data[,"element"] <- apply(data[,"rawelement",drop = F],1,function(x,...){stringr::str_extract_all(x,...)[[1]]},pattern = "[A-Z][a-z]*")
  data[,"num"] <- apply(data[,"rawelement",drop = F],1,function(x,...){ifelse(length(stringr::str_extract_all(x,...)[[1]])==0,1,stringr::str_extract_all(x,...)[[1]])},pattern = "[0-9]+")
  data[,"num"] <- as.numeric(data[,"num"])
  data[,"num"] <- data[,"num"]*nM
  
  adductElement <- unlist(strsplit(x = adductElement,split = "\\+"))
  for ( i in 1:length(adductElement)) {
    element <- stringr::str_extract_all(string = adductElement[i],pattern = "[A-Z][a-z]*[0-9]*")[[1]]
    data2 <- data.frame(rawelement = element)
    data2[,"element"] <- apply(data2[,"rawelement",drop = F],1,function(x,...){stringr::str_extract_all(x,...)[[1]]},pattern = "[A-Z][a-z]*")
    data2[,"num"] <- apply(data2[,"rawelement",drop = F],1,function(x,...){ifelse(length(stringr::str_extract_all(x,...)[[1]])==0,1,stringr::str_extract_all(x,...)[[1]])},pattern = "[0-9]+")
    data2[,"num"] <- as.numeric(data2[,"num"])
    if(grepl(pattern = "\\-",x = adductElement[i])){
      data2[,"num"] <- -data2[,"num"]
    }
    data <- rbind(data,data2)
  }
  
  data <- data %>% group_by(element) %>% summarise(num = sum(num))
  
  if(any(data$num < 0)){
    return(NA)
  }else{
    data <- as.data.frame(data[data$num!=0,])
    data$num <- as.character(data$num)
    data[data$num == "1","num"] <- ""
    newformula <- paste(apply(data, 1,paste, collapse = ""),collapse = "")
    return(newformula)
  }
}
