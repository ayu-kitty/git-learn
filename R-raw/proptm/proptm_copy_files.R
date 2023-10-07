#!/opt/conda/bin/Rscript
#' 复制联合分析所需文件
#'
#' @param ptmfile 修饰项目文件路径
#' @param profile 蛋白项目文件路径
#' @param path 运行结果保存路径
#' @export
copy_files<- function(ptmfile = './',profile= './',path = 'rawdata/',project='LM0001'){
  pacman::p_load(lmbio,dplyr,openxlsx)
  createdir(path)
  
  prodir <- list.files(path = profile,pattern = "^表达矩阵|差异表达矩阵\\(未筛选\\)",recursive = T)
  ptmdir <- list.files(path = ptmfile,pattern = "^表达矩阵|差异表达矩阵\\(未筛选\\)|classtype",recursive = T)
  
  for (r in prodir) {
    namem = strsplit(r,'/')[[1]][length(strsplit(r,'/')[[1]])]%>%gsub('.xlsx','-pro.xlsx',.)
    base::file.copy(from = paste0(profile,"/",r), to = paste0(path,namem))
  }
  for (t in ptmdir) {
    namem = strsplit(t,'/')[[1]][length(strsplit(t,'/')[[1]])]%>%gsub('.xlsx','-ptm.xlsx',.)
    base::file.copy(from = paste0(ptmfile,"/",t), to = paste0(path,namem))
  }
  ##FC筛选条件
  ratio_pro = readxlsx(paste0(profile,'/',
                              list.files(path = profile,pattern = "overview_information",recursive = T)),
                       sheet = '筛选条件')
  ratio_ptm = readxlsx(paste0(ptmfile,'/',
                              list.files(path = ptmfile,pattern = "overview_information",recursive = T)),
                       sheet = '筛选条件')
  colnames(ratio_pro)<- paste0(colnames(ratio_pro),'_pro')
  colnames(ratio_ptm)<- paste0(colnames(ratio_ptm),'_ptm')
  ratio = cbind(ratio_pro,ratio_ptm)
  savexlsx(ratio,paste0(path,'/ratio.xlsx'))
  
  initfile<- c('init','project_title')
  initdir<- ifelse(project%in%(strsplit(ptmfile,'/')%>%unlist),ptmfile,profile)
  init = readxlsx(paste0(initdir,'/',
                              list.files(path = profile,pattern = "overview_information",recursive = T)),
                       sheet = 'init')
  rownames(init)<- init$目录
  if(grepl('磷酸化',init$信息[1],ignore.case = F)){
    init['项目类型','信息'] = '蛋白组和磷酸化蛋白组联合分析'
  }else if(grepl('乙酰化',init$信息[1],ignore.case = F)){
    init['项目类型','信息'] = '蛋白组和乙酰化蛋白组联合分析'
  }else if(grepl('泛素化',init$信息[1],ignore.case = F)){
    init['项目类型','信息'] = '蛋白组和泛素化蛋白组联合分析'
  }else {
    init['项目类型','信息'] = '蛋白组和乙酰化蛋白组联合分析'
  }
  init['完成时间','信息'] = format(Sys.Date(),"%Y年%m月")
  savetxt(init,filename = paste0(path,'init.txt'),append = T,col.names=F,quote = F)
  
  name<-ifelse('联系人'%in%init[,1],paste0('-',init['联系人',2]),'')
  project_title <- paste0(init["项目编号",2],"-",init["客户名称",2],name,"-",init[1,'信息'],'结题报告')
  write(x = project_title,file = paste0(path,"/project_title.txt"))
  
  pro_title = read.xlsx(paste0(initdir,'/',
                              list.files(path = ptmfile,pattern = "overview_information",recursive = T)),
                       sheet = 'project_title',colNames = F)%>%as.character()
  
  grptm<- list.files(path = ptmfile,pattern = "Sample_Group",recursive = T)
  grpro<- list.files(path = profile,pattern = "Sample_Group",recursive = T)
  group<- list()
  for (i in 1:2) {
    group_pro<- readdata(paste0(profile,'/',grpro),sheet = i)
    group_ptm<- readdata(paste0(ptmfile,'/',grptm),sheet = i)
    group[[i]]<- merge(group_pro,group_ptm)
  }
  names(group)<- c('样品','比较组')
  #savexlsx(group,paste0(path,'/Sample_Group.xlsx'))
  pro<- getSheetNames(paste0(path,'差异表达矩阵(未筛选)-pro.xlsx'))
  ptm<- getSheetNames(paste0(path,'差异表达矩阵(未筛选)-ptm.xlsx'))
  if(all(pro==ptm)){
    gr<- gsub("-vs-",'/',pro)
    gro<- strsplit(pro,'-vs-')%>%unlist()%>%unique()
    new_gr<- setdiff(gr,group$比较组$比较组)
    project_gr<- sapply(new_gr,function(x){
      group$比较组$比较组[agrep(x,group$比较组$比较组)]
    })%>%as.data.frame()
    project_gr$new_group<- rownames(project_gr)
    group$比较组<- merge(group$比较组,project_gr,by.x='比较组',by.y='.',all.x=T)
  }
  savexlsx(group,paste0(path,'/Sample_Group.xlsx'))
}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-pt","--ptmfile",default = './', help = "修饰项目文件路径，例如：/data/hstore3/Projects/202306/ZLM2023010890")
  parser$add_argument("-pr","--profile",default = './', help = "蛋白项目文件路径，例如：/data/hstore3/Projects/202306/ZLM2023010891")
  parser$add_argument("-ip","--path",default = 'rawdata/' ,help = "运行结果保存路径,默认：'rawdata/")
  parser$add_argument("-i","--project",default = NULL ,help = "项目编号")
  args <- parser$parse_args()
  copy_files <- do.call(copy_files,args = args)
}
