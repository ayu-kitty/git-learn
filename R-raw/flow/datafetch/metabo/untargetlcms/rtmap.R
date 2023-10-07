
#' RT预测曲线
#'
#' @param lmdata 拟合数据
#'
#' @export
rtmap <- function(lmdata,summarydata){
  suppressMessages(library("ggplot2"))
  pp <- ggplot(lmdata)+
    geom_point(mapping = aes(x=x,y=y))+
    geom_abline(slope = lmdata$coefficients[2],intercept = lmdata$coefficients[1])+
    labs(x="",y="",title = paste0("y = ",
                                  round(lmdata$coefficients[2],digits = 4),"x+",
                                  round(lmdata$coefficients[1],digits = 4),"   R2 =",
                                  round(summarydata$r.squared,digits = 4)))+
    theme_bw()+
    theme(aspect.ratio = 1,
          text = element_text(family = "sans"),
          title = element_text(family = "sans"))
  
  ggsave(filename = "标准曲线.jpg",plot = pp)
}

