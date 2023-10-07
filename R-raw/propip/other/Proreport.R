#!/opt/conda/bin/Rscript

#' 项目记录
#' @export
prolog<-function(savepath="./"){
  prouser<-c("葛淑婷","康倩倩","潘静","李阿雨")
  names(prouser)<-c("shuting.ge","qianqian.kang","jing.pan","ayu.li")
  if(file.exists(paste0(savepath,"/项目执行记录.txt"))){
    log<-readdata(paste0(savepath,"/项目执行记录.txt"),header=F)
  }else log<-data.frame(id=c(1,2,3))
  
  if(file.exists("report/result/1.Project_information/overview_information.xlsx")){
    comn<-readxlsx("report/result/1.Project_information/overview_information.xlsx",sheet="diff_number") %>% nrow()
    oin<-readxlsx("report/result/1.Project_information/overview_information.xlsx",sheet=2)
    proid<-oin[3,2]
    if(is.na(grep(proid,log[,1])[1])){
      file_user <- prouser[file.info("report/")$uname]  # 获取文件的所有者用户名
      complete_time<-file.info(dir("./",pattern = ".*_report.html"))$mtime %>% gsub(" .*","",.) %>% gsub("-","/",.)
      sn<-unlist(readyaml("classfile.yaml")) %>% unique() %>% length
      samtype<-oin[oin[,1]=="物种及样本",2]
      protype<-gsub("分析","",oin[1,2]) %>% gsub("蛋白质组","",.) %>% gsub("定量","",.) %>% gsub("标记","",.)%>% gsub("纯","",.)%>% gsub("蛋白组","",.)
      tmp<-c(proid,samtype,protype,file_user,sn,comn,complete_time)
      savetxt(data = paste(tmp,collapse = "\t"),
              filename = paste0(savepath,"/项目执行记录.txt"),append = T)
    }
  }
}

#' 出报告函数
#'
#' @param rawpath 
#' @param reportpath 
#' @param mode 
#' @param cloud 
#' @param ... 
#' @export
Proreport<-function(rawpath="rawdata/",reportpath="report/",cloud="F",mode="LM",...){
  library(stringr,dplyr)
  bqplot(inputpath=paste0(rawpath,"pic/"))
  difdata <- readxlsx(paste0(reportpath,"/result/4.Different_expression/差异表达矩阵.xlsx"), sheet = 1)
  if(!"GO_term" %in% names(difdata)){
    diffanno(reportpath = reportpath)
  }
  
  init<-readxlsx(paste0(reportpath,"/result/1.Project_information/overview_information.xlsx"),sheet="init")
  projname=init[1,2]
  ptmname<-str_extract_all(projname,pattern = '磷酸化|泛素化|乙酰化|糖基化|乳酸化|宏蛋白')[[1]][1]
  if(!is.na(ptmname)){
    system(paste0("cp /data/hstore3/public/propip/report/",ptmname,"/* ",rawpath))
  }else{
    cgname=str_extract_all(projname,pattern = '超微量|单细胞|DIA|Label Free|Deep高深度')[[1]][1]
    bjname=str_extract_all(projname,pattern = 'TMT|iTRAQ')[[1]][1]
    if(!is.na(cgname)){
      if(cgname=="Deep高深度"){
        system(paste0("cp /data/hstore3/public/propip/report/DIA/* ",rawpath))
      }else if(cgname=="单细胞"){
        system(paste0("cp /data/hstore3/public/propip/report/超微量/* ",rawpath))
      }else system(paste0("cp '/data/hstore3/public/propip/report/",cgname,"/'* ",rawpath))
      
    }else if(!is.na(bjname)){
      system(paste0("cp /data/hstore3/public/propip/report/TMT/* ",rawpath))
    }else{
      system(paste0("cp /data/hstore3/public/propip/report/纯分析/* ",rawpath))
    }
  }
  confirmfile<-setdiff(dir(paste0(reportpath,"/result/1.Project_information/")),"overview_information.xlsx")[1]
  # 合并图片
  setwd(rawpath)
  system("cp -r /data/hstore3/public/propip/report/Mini-report/* ./")
  system("cp pic/* src/images/")
  if(mode=="auto"){
    ty<-str_extract_all(reportwd,pattern = 'LM|QD|OE|SG|HY')[[1]][1]
    if(is.na(ty)){
      mode<-"LM"
    }else mode<-ty
  }
  if(!is.na(confirmfile)){
    cof<-readxlsx(paste0("../report/result/1.Project_information/",confirmfile))
    if("客户经理" %in% cof[,1]){
      if(cof[cof[,1]=="客户经理",2]=="张志永01"){
        mode="qds"
      }
    }
    
  }
  # 压缩keggmap
  if(dir.exists("../report/result/5.Enrichment/KEGG_map/")){
    path0=getwd()
    setwd("../report/result/5.Enrichment/")
    system("zip -qrm KEGG_map.zip KEGG_map")
    setwd(path0)
  }
  
  # 正式报告
  # system(paste0("python Report.py -m ",mode," -c ",cloud," > ../project_report_stderr.log"))
  system(paste0("python Report.py -m ",mode," -c ",cloud," >/dev/null"))
  # 简版报告
  system(paste0("python simreport.py -m ",mode," >/dev/null"))

  if(cloud=="T"){
    setwd("..")
    loctime=format(Sys.Date(),format="%Y%m")
    prolog(savepath = paste0("/data/hstore3/localPro/1.Projects/",loctime,"/"))
    locpath=paste0("/data/hstore3/localPro/1.Projects/",loctime,"/",init[3,2])
    minipath=paste0("/data/hstore3/localPro/5.Mini-report/",format(Sys.Date(),format="%Y"),"/",loctime)
    reportname=dir("./",pattern = ".zip")
    mininame=dir("./",pattern = "_report.html",full.names = T)
    file.remove(reportname,locpath,overwrite = T)
    file.copy(mininame,minipath,overwrite = T)
    system(paste0("cd ",locpath,";rm -rf ",gsub(".zip","",reportname),"; unzip -q ",reportname," -d ", gsub(".zip","",reportname)))
  }
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  parser$add_argument("-m","--mode", default = "auto",help = "报告模板类型，包括LM、OE、QD、SG、HY、GY、YZ、unlogo，默认为自动判断")
  parser$add_argument("-c","--cloud", default = "F",help = "报告是否云交付，默认为F")
  parser$add_argument("-rd","--rawpath", default = "rawdata/",help = "存放出报告文件路径，默认为rawdata/")
  parser$add_argument("-rp","--reportpath", default = "report/",help = "报告出具路径，默认为./report/")
  args <- parser$parse_args()
  Proreport <- do.call(Proreport,args = args)
}
