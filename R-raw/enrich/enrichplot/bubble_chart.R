#!/opt/conda/bin/Rscript
#' @export
bubblescale<-function(data){
  maxn <- max(data)
  minn <- min(data)
  if ((maxn - minn) == 0) {
    x <- maxn
    lab <- c(x)
    ran <- c(3, 3)
  } else if ((maxn - minn) == 1) {
    lab <- c(minn, minn + 1)
    ran <- c(3, 4)
  } else if ((maxn - minn) == 2) {
    lab <- c(minn, minn + 1, minn + 2)
    ran <- c(3, 5)
  } else if ((maxn - minn) %% 5 == 0 || (maxn - minn + 1) %% 5 == 0) {
    x <- round((maxn - minn) / 5)
    lab <- c(minn, minn + x, minn + 2 * x, minn + 3 * x, minn + 4 * x, minn + 5 * x)
    ran <- c(1, 5)
  } else if ((maxn - minn) %% 4 == 0 || (maxn - minn + 1) %% 4 == 0) {
    x <- round((maxn - minn) / 4)
    lab <- c(minn, minn + x, minn + 2 * x, minn + 3 * x, minn + 4 * x)
    ran <- c(1, 5)
  } else if ((maxn - minn) %% 3 == 0 || (maxn - minn + 1) %% 3 == 0) {
    x <- round((maxn - minn) / 3)
    lab <- c(minn, minn + x, minn + 2 * x, minn + 3 * x)
    ran <- c(2, 5)
  } else {
    x <- floor((maxn - minn) / 4)
    lab <- c(minn, minn + x, minn + 2 * x, minn + 3 * x, minn + 4 * x)
    ran <- c(1, 5)
  }
  return(list(lab,ran))
}
#' 富集气泡图
#'
#' @param savepath 保存路径
#' @param type 数据库类型
#' @param indata 富集表格名称，默认不输入为None
#' @param incompare 比较组名称
#' @param filt 数据是否筛选
#' @param grid 是否进行分面
#' @param number 绘制条目数量
#' @param height 高度
#' @param width 宽度
#' @param dpi 分辨率
#' @param fontfamily 字体类型 
#' @param ... 
#'
#' @export
bubble<-function(savepath = "./", type = NULL, indata="None",incompare = NULL, filt = "T", grid = "T", number =
                   20,imagetype = c("pdf","png") ,height = 8,width = 14,dpi=300,fontfamily="sans",col="solarExtra",...){
  pacman::p_load(dplyr,ggplot2,stringr,Hmisc)
  mycol<- rev(SelectColors(col))
  
  if(indata=="None"){
    comenrichpath<-paste0(savepath,getenrichtype(type)$filn,"/",incompare,"/")
    enrichfiles<-paste0(comenrichpath,dir(path =comenrichpath, pattern = "enrichment-*"))
    savepath = "./enrich/"
  }else{
    comenrichpath<-paste0(dirname(indata),"/")
    enrichfiles<-indata
  }
  typef <- ifelse(is.null(type),"enrichment",getenrichtype(type)$filn)
  fixed_col <- c("p-value","ListHits","Enrichment_score","Term")
  for(enrichfile in enrichfiles){
    tudd<-ifelse(grepl("-Total",enrichfile),"Total",ifelse(grepl("-Up",enrichfile),"Up",ifelse(grepl("-Down",enrichfile),"Down","Total")))
    endat<-readdata(enrichfile)
    file_col <- colnames(endat)
    stopname <- setdiff(fixed_col,file_col)
    if(length(stopname)>0){stop(paste0("缺少",paste0(stopname,collapse = "、"),"列或大小写名称不对应"))}
    endat$`p-value`[endat$`p-value`==0]<-min(endat$`p-value`[endat$`p-value`!=0])
    if(filt == "T"){
      if(grid == "T"){
        if(!"Category" %in% file_col){stop("选择分面时，输入表需存在Category列！请检查Category列是否存在或名称是否对应")}
        endata_filt<-rbind(filter(endat,ListHits!=1,Category=="molecular_function")[1:5,],
                           filter(endat,ListHits!=1,Category=="cellular_component")[1:5,],
                           filter(endat,ListHits!=1,Category=="biological_process")[1:5,])
        endata<-endata_filt %>% arrange(Category, Enrichment_score)
      }else{
        endata_filt<-endat[order(endat$`p-value`),][1:number,]
        endata<-endata_filt %>% arrange(Enrichment_score)
      }
      na.omit(endata)->endata
      savexls(endata,paste0(comenrichpath,typef,".Bubble.",tudd,".xls"))      
    }else {
      endata<-endat
      na.omit(endata)->endata
    }
    if(nrow(endata)>2){
      endata$Term<-capitalize(as.character(endata$Term))
      endata$Term <- factor(endata$Term, levels=endata$Term)
      if(grid == "T"){
        endata$Category<-str_replace_all(endata$Category,c("biological_process"="BP","cellular_component"="CC","molecular_function"="MF"))
        
        p<-ggplot(endata,aes(Enrichment_score,Term))+
          geom_point(aes(size=ListHits,color=endata$`p-value`))+
          facet_grid(Category ~ .,scales="free")+
          guides(size = guide_legend(order=1))+
          labs(x = "Enrichment_score", y = "",
               title = ifelse(is.null(incompare),paste0("Top ",typef," Terms"),paste0(incompare,"(",tudd,")  \n Top ",typef," Terms")),
               color="P-value",size="Count") +
          scale_colour_gradientn(colours = mycol)+
          theme_bw()+
          theme(aspect.ratio=19/21)+
          theme(strip.background=element_rect(fill="gray60",color="gray60", linewidth=1, linetype="solid")) +
          theme(panel.border = element_rect(color = "black", linewidth = 1, fill = NA)) +
		  scale_size_continuous(breaks = bubblescale(endata$ListHits)[[1]], range = bubblescale(endata$ListHits)[[2]])+
          scale_y_discrete(labels=function(x)ifelse(nchar(x)>60,paste0(substr(x,1,60),"..."),x))+
          theme(plot.title = element_text(size=18,hjust=0.5)) +
          theme(plot.margin=unit(c(2,0,2,0), "lines"))+
          theme(axis.text.x = element_text(size=15,color="black"),strip.text = element_text(size = 15,color="black"),axis.text.y = element_text(size = 15,color="black"),axis.title = element_text(size=18,color="black"),
                legend.text=element_text(size=15),legend.title=element_text(size=18))
        
      }else{
        p <- ggplot(endata,aes(Enrichment_score,Term))+
          geom_point(aes(size=ListHits,color=endata$`p-value`)) +
          scale_y_discrete(labels=function(x)ifelse(nchar(x)>60,paste0(substr(x,1,60),"..."),x))+
          theme_bw()+
          theme(panel.border = element_rect(color = "black", linewidth = 1, fill = NA)) +
          theme(axis.line = element_blank(), strip.background = element_blank()) +
          theme(aspect.ratio=16/9)+
          theme(axis.text.x=element_text(size=15,color="black"),axis.text.y=element_text(size=15,color="black"),axis.title=element_text(size=18,color="black"),
                legend.text=element_text(size=15),legend.title=element_text(size=18))+
          scale_colour_gradientn(colours = mycol)+
		  scale_size_continuous(breaks = bubblescale(endata$ListHits)[[1]], range = bubblescale(endata$ListHits)[[2]])+
          guides(size = guide_legend(order = 1)) + 
          theme(plot.title = element_text(hjust = 0.5,size=18))+
          theme(plot.margin=unit(c(2,2,2,2), "lines"))+
          labs(y ='',
               title = ifelse(is.null(incompare),paste0("Top ",typef," Terms"),paste0(incompare,"(",tudd,")  \n Top ",typef," Terms")),
               color="P-value",size="Count") 	
        
      }
      ggplotsave(plot = p,savepath=comenrichpath,
                 mapname = paste0(typef,".Bubble.",tudd),
                 width = width,
                 height = height,
                 imagetype = imagetype,
                 family=fontfamily,
                 ...)
      
      
    }else{
      warning("Term数量低于3，不提供富集气泡图")
      savetxt(data = "Term数量低于3，不提供富集气泡图",
              filename = paste0(comenrichpath,"/说明.txt"),append = T)
      return()
    }
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-i","--imagetype",default = c("pdf","png"), help = "图片格式")
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-wi","--width",default = 14, type= "double",help = "图片宽度")
  parser$add_argument("-he","--height",default = 8, type= "double",help = "图片高度")
  parser$add_argument("-d","--dpi",default = 300, type= "double",help = "分辨率")
  parser$add_argument("-col","--col",default = "solarExtra" ,help = "气泡图色系")
  # 此图参数
  
  parser$add_argument("-t","--type",type="character", default=NULL, help="enrich database class,such as G/K/W/R/I", metavar="character")
  parser$add_argument("-ic","--incompare",default = NULL, help = "比较组名称，默认NULL")
  parser$add_argument("-id","--indata",default = "None", help = "绘图数据表格文件（包含路径），默认None不输入")
  parser$add_argument("-f","--filt",default = "T", type= "character",help = "数据是否需要筛选,默认T")
  parser$add_argument("-g","--grid", type = "character", default = "F",help = "是否需要分面，默认为F")
  parser$add_argument("-n","--number",default = 20,help = "top number each category,默认20")
  parser$add_argument("-s","--savepath",default = "./", help = "富集结果保存路径,默认./,流程中默认./enrich/")
  
  args <- parser$parse_args()

  bubbleplot <- do.call(bubble,args = args)
}