#!/opt/conda/bin/Rscript

#' movemap
#'
#' 售后色谱图转移
#'
#' @export

remap_gclc<-function(oldpath_lc=paste0(path,name,"/",dataall["项目报告",2],"/2.色谱图/LCMS/"),
                     oldpath_gc=paste0(path,name,"/",dataall["项目报告",2],"/2.色谱图/GCMS/"),
                    path_lc=paste0("./",dataall["项目报告",2],"/2.色谱图/LCMS/"),
                    path_gc=paste0("./",dataall["项目报告",2],"/2.色谱图/GCMS/"),
                    rename=rename,
                     ...){
  for (i in 1:dim(rename)[1]) { 
  if(is.na(rename$deal[i])){
  if(rename$oldname[i]!=rename$newname[i]){
    lmbio::createdir(filename = "rename",force = force)
    oldmap_gc<-list.files(pattern =rename$oldname[i] ,path = oldpath_gc,full.names = T,all.files = T,recursive = T)
    file.copy(from = oldmap_gc ,to="./rename/GCMS/",recursive=T)
    oldname_gc<-list.files(pattern =rename$oldname[i] ,path = "./rename/GCMS/",full.names = T)
    newname_gc = gsub(pattern = rename$oldname[i],replacement =rename$newname[i] ,x = oldname_gc )
    file.rename(oldname_gc,newname_gc)
    file.copy(from = list.files(path = "./rename/GCMS/",full.names = T) ,to=path_gc,recursive=T)
    unlink(x = "rename",recursive = T)
    #lc
    oldmap_lc<-list.files(pattern =rename$oldname[i] ,path = oldpath_lc,full.names = T,all.files = T,recursive = T)
    file.copy(from = oldmap_lc ,to="./rename/LCMS/",recursive=T)
    oldname_lc<-list.files(pattern =rename$oldname[i] ,path = "./rename/LCMS/",full.names = T)
    newname_lc = gsub(pattern = rename$oldname[i],replacement =rename$newname[i] ,x = oldname_lc )
    file.rename(oldname_lc,newname_lc)
    file.copy(from = list.files(path = "./rename/LCMS/",full.names = T) ,to=path_lc,recursive=T)
    unlink(x = "rename",recursive = T)
  }else{
    map_gc<-list.files(pattern =rename$oldname[i] ,path = oldpath_gc,full.names = T,recursive = T)
    file.copy(from = map_gc ,to=path_gc,recursive=T)
    map_lc<-list.files(pattern =rename$oldname[i] ,path = oldpath_lc,full.names = T,recursive = T)
    file.copy(from = map_lc ,to=path_lc,recursive=T)
   }
  }
}
}

#' @export
move_map<-function(mappath=paste0(path,name,"/",dataall["项目报告",2]),
                  movepath=paste0("./",dataall["项目报告",2],"/2.色谱图/"),
                  ...){
  map<-list.files(dir(mappath,pattern = "色谱图$",full.names = T),full.names = T)
  file.copy(from = map ,to=movepath,recursive=T)
}

#' @export
remap<-function(oldpath=paste0(path,name,"/",dataall["项目报告",2],"/2.色谱图/"),
                newpath=paste0("./",dataall["项目报告",2],"/2.色谱图"),
                rename=rename,
                ...){
  for (i in 1:dim(rename)[1]) { 
    if(is.na(rename$deal[i])){
    if(rename$oldname[i]!=rename$newname[i]){
    lmbio::createdir(filename = "rename",force = force)
    oldmap<-list.files(pattern =rename$oldname[i] ,path = oldpath,full.names = T,all.files = T,recursive = T)
    file.copy(from = oldmap ,to="./rename/",recursive=T)
    oldname<-list.files(pattern =rename$oldname[i] ,path = "./rename/",full.names = T)
    newname = gsub(pattern = rename$oldname[i],replacement =rename$newname[i] ,x = oldname )
    file.rename(oldname,newname)
    file.copy(from = list.files(path = "./rename/",full.names = T) ,to=newpath,recursive=T)
    unlink(x = "rename",recursive = T)
  }else{
    map<-list.files(pattern =rename$oldname[i] ,path = oldpath,full.names = T)
    file.copy(from = map ,to=newpath,recursive=T)
    }
   }
  }
}
