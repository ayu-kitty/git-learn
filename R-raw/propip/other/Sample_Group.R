#!/opt/conda/bin/Rscript

#' 提取蛋白分析所需样本信息表
#'
#' @param inputpath 
#' @export
prosamgroup<-function(inputpath="./"){
  pacman::p_load(dplyr,openxlsx,stringr)
  ##建立分类表
  project_type<- data.frame(
    c('TMT','single_mark','标记定量蛋白质组分析'),
    c('iTRAQ','single_mark','标记定量蛋白质组分析'),
    c('4D-DIA','DIA','定量蛋白质组分析'),
    c('3D-DIA','DIA','定量蛋白质组分析'),
    c('4D-Label Free','label_free','定量蛋白质组分析'),
    c('3D-Label Free','label_free','定量蛋白质组分析'),
    c('超微量蛋白质组','DIA','超微量蛋白质组分析'),
    c('4D-Label Free宏蛋白质组','label_free','宏蛋白质组分析'),
    c('4D-DIA宏蛋白质组','DIA','宏蛋白质组分析'),
    c('3D-Label Free宏蛋白质组','label_free','宏蛋白质组分析'),
    c('3D-DIA宏蛋白质组','DIA','宏蛋白质组分析'),
    c('4D-Olink蛋白质组','DIA','蛋白质组分析'),
	c('Deep高深度血液蛋白质组','DIA','高深度血液蛋白质组分析')
  )%>%t()%>%as.data.frame()%>%dplyr::rename(project_type=V1,analyze_type=V2)
  rownames(project_type)<- project_type[,1]
  PTM_type<- data.frame(
    c('磷酸化','phos'),
    c('泛素化','ubiqu'),
    c('乙酰化','acety'),
    c('糖基化','deamid'),
    c('乳酸化','lactyl')
  )%>%t()%>%as.data.frame()%>%dplyr::rename(PTM_type=V1,analyze_type=V2)
  rownames(PTM_type)<- PTM_type[,1]
  
  if(!file.exists(paste0(inputpath,"Sample_Group.xlsx"))){
    filename <- list.files(inputpath,pattern="(DLM|ZLM|LM|DQD|DOE|DZLM|DZQD).*.xlsx")
    if(length(filename ) == 0){
     stop("未找到项目登记单") 
    }
    init <- readdata(paste0(inputpath,filename),sheet = 1,row.names = 1)
    #######生成init文件#####
    init_sample <- data.frame()
    init_sample[1:7,1] <- c("项目类型","任务单号","项目编号","客户名称","联系人","物种及样本","完成时间")
    init_sample[1,2] <- ifelse(init["修饰类型",1]=="无",
                               paste0(str_extract(project_type[init["项目类别",1],1],"[^\\p{Han}]+")%>%na.omit(),project_type[init["项目类别",1],3]),
                               paste0(str_extract(project_type[init["项目类别",1],1],"[^\\p{Han}]+")%>%na.omit(),init["修饰类型",1],project_type[init["项目类别",1],3]))
    if(is.na(grep("-b",init["项目编号",1])[1])){
		init_sample[2,2] <- paste0(init["项目编号",1],"-b1")
	}else init_sample[2,2] <- init["项目编号",1]
    init_sample[3,2] <- init["项目编号",1]
    init_sample[4,2] <- init["客户名称",1]
    init_sample[5,2] <- init["联系人",1]
    init_sample[6,2] <- ifelse(is.na(grep(init["样本物种",1],init["样本类型",1])[1]),paste0(init["样本物种",1],init["样本类型",1]),grep(init["样本物种",1],init["样本类型",1],value = T))
    init_sample[7,2] <- format(Sys.Date(),"%Y年%m月")
    if(init_sample[5,2]==init_sample[4,2]|is.na(init_sample[5,2])){
      init_sample<- init_sample[-5,]
    }
    savetxt(init_sample, paste0(inputpath,"/init.txt"),quote=FALSE,col.names=F)
    
    ##生成project_tittle文件
    rownames(init_sample)<- init_sample[,1]
    bb<- gsub("分析","结题报告",init_sample["项目类型",2])%>%
      gsub(" ","_",.)
    name<-ifelse('联系人'%in%init_sample[,1],paste0('-',init_sample['联系人',2]),'')
    PTM<- ifelse(init["修饰类型",1]=="无","",init["修饰类型",1])
    project_title <- paste0(init_sample["项目编号",2],"-",init_sample["客户名称",2],name,"-",bb)
    write(x = project_title,file = paste0(inputpath,"/project_title.txt"))
    
    sample1 <- readdata(paste0(inputpath,filename),sheet = '样本分组信息')
    sample2 <- readdata(paste0(inputpath,filename),sheet = '差异比较信息')
    ##Sample_Group sheet1
    if (nrow(sample1)==0&nrow(sample2)==0){
      samp = readdata(paste0(inputpath,filename),sheet = '样本信息')
      sample_group = as.data.frame(matrix(data = NA ,nrow = nrow(samp),ncol = 2))
      sample_group[,1]<- samp$样本分析名称
      sample_group[,2]<- samp$样本分析名称
      colnames(sample_group) <- c("Sample","Group")
    }else{
      sample_group <- data.frame()
      for (i in 1:length(sample1[,1])) {
        sample_data <- data.frame()
        n <- length(strsplit(sample1[i,2],",")[[1]])
        if (n==1) {
          sample_data[1,1] <- sample1[i,2]
          sample_data[1,2] <- sample1[i,1]
          sample_group <- rbind(sample_group,sample_data)
        }else{
          sample_data[1:n,1] <- strsplit(sample1[i,2],",")[[1]][1:n]
          sample_data[1:n,2] <- sample1[i,1]
          sample_group <- rbind(sample_group,sample_data)
        }
      }
      colnames(sample_group) <- c("Sample","Group")
    }
    
    ##Sample_Group sheet2
    diff_group <- data.frame()
    for (j in 1:length(sample2[,1])) {
      diff_group[j,1] <- paste0(sample2[j,1],"/",sample2[j,2])
      diff_group[j,2] <- sample2$paired[j]
    }
    diff_group<- na.omit(diff_group)
    colnames(diff_group) <- c("比较组","配对")
    list_group <- list(sample_group,diff_group)
    names(list_group) <- c("样品","比较组")
    list_group=list_group[lapply(list_group,nrow)!=0]
    savexlsx(list_group, paste0(inputpath,"Sample_Group.xlsx"))
  }
}
if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-ip","--inputpath",type="character", default="./", help="下机数据文件存放路径，默认当前路径", metavar="character")
  
  args <- parser$parse_args()
  sgfetch_result <- do.call(prosamgroup,args = args)
}

