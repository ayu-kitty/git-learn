#' @export
getspecies_tissuedata <- function(species_tissue = "小鼠",
                                  copypath = "./sample/qualitative"){
  
  species <- unlist(strsplit(x = species_tissue,split = "_"))[1]
  species_c <- switch (species,
                       "人"="hsa",
                       "小鼠"="mmu",
                       "大鼠"="rno",
                       "none")
  tissue <- unlist(strsplit(x = species_tissue,split = "_"))[2]
  tissue_c <- switch (tissue,
                      "心脏" = "heart",
                      "脑" = "brain",
                      "肾脏" = "kidney",
                      "肝脏" = "liver",
                      "肿瘤" = "tumor",
                      "肺" = "lung",
                      NULL)
  
  tissue_hmdb <- switch (tissue,
                         "心脏" = "heart",
                         "脑" = "Brain",
                         "肾脏" = "Kidney",
                         "肝脏" = "Liver",
                         "肿瘤" = "tumor",
                         "肺" = "Lung",
                         "脾脏" = "Spleen",
                         "胰" = "Intestine",
                         "骨骼肌" = "Skeletal Muscle",
                         "平滑肌" ="Smooth Muscle",
                         "表皮" = "Epidermis",
                         "脂肪组织" = "Adipose Tissue",
                         NULL)
  lcqualitativefrom <- "none"
  manqualitativefrom <- "none"
  hmdblocation <- "none"
  
  if(length(tissue_hmdb) != 0){
    if(!is.na(tissue_hmdb)){
      hmdblocation <- getmysqldata(table = "hmdblocation",
                                   tablelist = "hmdbid",
                                   wherename = "location",
                                   wheredata = tissue_hmdb)[,1]
    }
  }
  
  if(species_c != "none"){
    path <- list.files(path = databasepath(database = "database/database/qualitative/Space/LCMS",path = paste0(species,"/",tissue),exist = F),
                       pattern = "多项目非靶定性统计.xlsx",
                       full.names = T)
    if(length(path) == 0){
      path <- list.files(path = databasepath(database = "database/database/qualitative/Space/LCMS",path = species),
                         pattern = "多项目非靶定性统计.xlsx",
                         full.names = T)
    }
    if(length(path) >= 1){
      file.copy(from = path[1],to = copypath)
      lcqualitativefrom <- paste0(copypath,"/",basename(path[1]))
    }
  }
  
  if(length(tissue_c) != 0){
    if(!is.na(tissue_c)){
      path <- list.files(path = databasepath(database = "database/database/qualitative/Space/人工"),
                         pattern = paste0("人工数据库-",tissue_c,".xlsx"),
                         full.names = T)
      if(length(path) >= 1){
        file.copy(from = path[1],to = copypath)
        manqualitativefrom <- paste0(copypath,"/",basename(path[1]))
      }
    }
  }else if(species_c != "none"){
    path <- list.files(path = databasepath(database = "database/database/qualitative/Space/人工"),
                       pattern = "人工数据库.xlsx",
                       full.names = T)
    if(length(path) >= 1){
      file.copy(from = path[1],to = copypath) 
      manqualitativefrom <- paste0(copypath,"/",basename(path[1]))
    }
  }
  
  return(list(hmdbloaction = hmdblocation,
              lcqualitativefrom = lcqualitativefrom,
              manqualitativefrom = manqualitativefrom))
}