#!/opt/conda/bin/Rscript

# p值及adj.p计算
#' @export
p_ana <- function(data,
                  class,
                  paired = F,
                  p_method = "t.test",
                  p_adjust_methods = "BH",
                  ...){
  
  if(any(is.na(data))){
    for ( i in 1:dim(data)[1]) {
      for ( group in unique(class[,1])) {
        if(all(is.na(data[i,row.names(class[class[,1] == group,,drop = F])]))){
          data[i,row.names(class[class[,1] == group,,drop = F])] <- min(as.matrix(data),na.rm = T)
        }
      }
    }
  }
  
  p.value <- apply(data,1,calp,
                   class = class,
                   paired = paired,
                   p_method = p_method,
                   ...)
  
  if(all(p.value == 1)){
    return(NULL)
  }else if(all(is.na(p.value))){
    return(NULL)
  }
  
  diffdata <- as.data.frame(p.value)
  colnames(diffdata) <- "p-value"
  diffdata[,"q-value"] <- p.adjust(diffdata$`p-value`,method = p_adjust_methods)
  row.names(diffdata) <- row.names(data)
  
  return(diffdata)
}


# p值计算
calp <- function(data,
                 class,
                 p_method = "t.test",
                 ...){
  data <- as.data.frame(data)
  colnames(data) <- "value"
  
  groupdata <- merge(data,class,by=0)
  
  if(length(unique(groupdata$Group))>2){
    pfun <- switch (p_method,
                    "oneway.test" = oneway.test.method,
                    "kruskal.test" = kruskal.test.method,
                    oneway.test.method)
  }else if(length(unique(groupdata$Group))==2){
    pfun <- switch (p_method,
                    "t.test"= t.test.method,
                    "f.test"= f.test.method,
                    "ft.test"= ft.test.method,
                    "wilcox.test" = wilcox.test.method,
                    "oneway.test" = oneway.test.method,
                    "kruskal.test" = kruskal.test.method,
                    t.test.method)
  }else{
    stop()
  }
  
  p.value <- tryCatch(
    pfun(data = groupdata,formula = value~Group,...),
    error = function(e){
      # print(e)
      # warning("p值运算报错，返回1",immediate. = T)
      return(1)
    }
  )

  return(p.value)
}

# T检验
t.test.method <- function(data,
                          formula,
                          paired = F,
                          ...){
  data <- t.test(formula,
                 data = data,
                 var.equal = T,
                 paired = paired,
                 ...)
  return(data$p.value)
}

# F检验
f.test.method <- function(data,
                          formula,
                          paired = F,
                          ...){
  data <- t.test(formula,
                 data = data,
                 var.equal = F,
                 paired = paired,
                 ...)
  return(data$p.value)
}

# 根据方差齐性进行F检验或T检验
ft.test.method <- function(data,
                           formula,
                           paired = F,
                           ...){
  # 方差齐性检验
  sdi <- car::leveneTest(y = formula,
                         data = data)$`Pr(>F)`[1]
  
  if (sdi >= 0.05) {
    p.value <- t.test.method(data = data,
                             formula = formula,
                             paired = paired,
                             ...)
  } else {
    p.value <- f.test.method(data = data,
                             formula = formula,
                             paired = paired,
                             ...)
  }
  
  return(p.value)
}


# 秩和检验
wilcox.test.method <- function(data,
                               formula,
                               paired = F,
                               ...){
  data <- wilcox.test(formula,
                      data = data,
                      paired = paired,
                      ...)
  return(data$p.value)
}


# 单变量统计
oneway.test.method <- function(data,
                               formula,
                               paired = F,
                               var.equal = T,
                               ...){
  data <- oneway.test(formula,
                      data = data,
                      var.equal = var.equal,
                      ...)
  return(data$p.value)
}

# Kruskal-Wallis检验
kruskal.test.method <- function(data,
                                formula,
                                paired = F,
                                ...){
  data <- kruskal.test(formula,
                       data = data,
                       ...)
  return(data$p.value)
}
