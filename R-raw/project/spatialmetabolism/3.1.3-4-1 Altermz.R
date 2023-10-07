#' Altermz
#' 修改定性结果添删mz

#' @param samplepath sample路径
#' @param peakpath peak路径
#' @param qualitativepath qualitative路径
#' @param addmz 添加的mz
#' @param delmz 删除的mz
#' @param addadducts 添加的加和离子
#' @param addformula 添加的分子式
#' @param addMeta 添加的代谢物名
#' @param version 定性版本
#' @param mode 离子模式
#' 
#' @export
Altermz<-function(samplepath="./sample/",
                peakpath=paste0(samplepath,"peak/"),
                qualitativepath=paste0(samplepath,"qualitative/"),
                addmz=NULL,
                delmz=NULL,
                addadducts=NULL,
                addformula=NULL,
                addMeta=NULL,
                version="2023new",
                mode="neg"
){
  #删mz
  if(length(delmz)>0){
    if(version=="2023new"){
      filename=paste0(qualitativepath,"adducts.xlsx")
    }else if(version=="2023old"){
      filename=paste0(peakpath,"样本离子统计.xlsx")
    }else{
      stop(paste0(version,"版本不对"))
    }
    
    qdata<-readdata(filename=filename,sheet=mode)
    allmz<-format(qdata$mz, nsmall = 5, trim = T)
    delmz<-format(delmz, nsmall = 5, trim = T)
    for(mz in delmz){
      if(mz %in% allmz){
        qdata$sampleselect[which(allmz==delmz)]<-FALSE
      }else{
        warning(paste0(mz,"不存在于样本离子统计中"), immediate. = T)
      }
    }
    savexlsx1(qdata,filename=filename,sheet=mode)
  }
  

  
  if(length(addmz)>0){
     
    #添加mz
    if(version=="2023new"){
      filename=paste0(qualitativepath,"adducts.xlsx")
    }else if(version=="2023old"){
      filename=paste0(peakpath,"样本离子统计.xlsx")
    }else{
      stop(paste0(version,"版本不对"))
    }
    
    qdata<-readdata(filename=filename,sheet=mode)
    allmz<-format(qdata$mz, nsmall = 5, trim = T)
    addmz<-format(addmz, nsmall = 5, trim = T)
    for(mz in addmz){
      if(mz %in% allmz){
        qdata$sampleselect[which(allmz==addmz)]<-TRUE
      }else{
        warning(paste0(mz,"不存在于样本离子统计中"), immediate. = T)
      }
    }
    savexlsx1(qdata,filename=filename,sheet=mode)
      
    
      #添加代谢物信息
      qdata<-readdata(filename=paste0(qualitativepath,"adducts.xlsx"),sheet=mode)
      annodata<-readdata(filename=paste0(qualitativepath,"anno.xlsx"),sheet=mode)
      annoalldata<-readdata(filename=paste0(qualitativepath,"anno-all.xlsx"),sheet=mode)
      if(version=="2023new"){
        annoalldata<-annoalldata[,-2]
      }
      allmz<-format(qdata$mz, nsmall = 5, trim = T)
      addmz<-format(addmz, nsmall = 5, trim = T)
      for(mz in addmz){
        if(mz %in% allmz){
          if(qdata$samplemeandata[which(allmz==mz)]!=0){
            adducts=addadducts[which(addmz==mz)]
            formula=addformula[which(addmz==mz)]
            meta=addMeta[which(addmz==mz)]
            qdata$adducts[which(allmz==mz)]<-adducts
            qdata$formula[which(allmz==mz)]<-formula
          }else{
            warning(paste0(mz,"表达量为0，故不添加"), immediate. = T)
            next
          }
        
        
        if(formula %in% annodata$Formula){
          formuladata<-annodata[annodata$Formula==formula,]
          if(!(meta %in% formuladata$Metabolites)){
            formulaalldata<-annoalldata[annoalldata$Formula==formula,]
            if(meta %in% formulaalldata$Metabolites){
              metanum=which(formulaalldata$Metabolites==meta)
              formulaalldata$num[metanum]=unique(formuladata$num)
              annodata<-rbind(annodata,formulaalldata[metanum,])
              annodata$num[which(annodata$Formula==formula)]=annodata$num[which(annodata$Formula==formula)]+1
            }else{
              warning(mz,"对应的代谢物信息",meta,"需要手动添加到anno.xlsx中")
            }
          }
        }else{
          warning(mz,"对应的分子式信息",formula,"需要手动添加到anno.xlsx中", immediate. = T)
        }
          

        }else{
          warning(paste0(mz,"不存在于样本离子统计中"), immediate. = T)
        }
        
      }
      savexlsx1(annodata,filename=paste0(qualitativepath,"anno.xlsx"),sheet=mode)
      savexlsx1(qdata,filename=paste0(qualitativepath,"adducts.xlsx"),sheet=mode)
    
  }
}
