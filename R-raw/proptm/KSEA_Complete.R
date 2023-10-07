#!/opt/conda/bin/Rscript
#' 激酶活性计算
#' @export
KSEA_Complete<-function (KSData, phos, m.cutoff,p.cutoff,path)
{
  new = phos
  colnames(new)<- c("protein","SUB_MOD_RSD","SUB_GENE","log2FC",'p')
  new = new[complete.cases(new$log2FC), ]
  
  KSData.filtered = KSData
  KSData.dataset = merge(KSData.filtered, new)
  KSData.dataset = KSData.dataset[order(KSData.dataset$GENE),]
  KSData.dataset.abbrev = KSData.dataset[, c("GENE","KINASE","SUB_GENE", "SUB_MOD_RSD","log2FC","p" )]
  colnames(KSData.dataset.abbrev) = c("Kinase.Gene","Kinase.Name","Substrate.Gene","Substrate.Mod","log2FC","p")
  KSData.dataset.abbrev = KSData.dataset.abbrev[order(KSData.dataset.abbrev$Kinase.Gene,
                                                      KSData.dataset.abbrev$Substrate.Gene, KSData.dataset.abbrev$Substrate.Mod,
                                                      KSData.dataset.abbrev$p), ]
  if (nrow(KSData.dataset.abbrev)>1) {
    KSData.dataset.abbrev = aggregate(log2FC ~ Kinase.Gene + Kinase.Name +
                                        Substrate.Gene + Substrate.Mod , data = KSData.dataset.abbrev,
                                      FUN = mean)
  }else{
    KSData.dataset.abbrev = KSData.dataset.abbrev
  }
  
  KSData.dataset.abbrev = KSData.dataset.abbrev[order(KSData.dataset.abbrev$Kinase.Gene),]
  
  kinase.list = as.vector(KSData.dataset.abbrev$Kinase.Gene)
  kinase.list = as.matrix(table(kinase.list))
  if (nrow(kinase.list)>1) {
    Mean.FC = aggregate(log2FC ~ Kinase.Gene, data = KSData.dataset.abbrev,FUN = mean)##log2FC的计算方式为所有激酶对应底物的log2FC均值
  }else{
    Mean.FC = KSData.dataset.abbrev[,c('Kinase.Gene','log2FC')]
  }
  
  Mean.FC = Mean.FC[order(Mean.FC[, 1]), ]
  Mean.FC$mS = Mean.FC[, 2]##mS为log2FC
  Mean.FC$Enrichment = Mean.FC$mS/abs(mean(new$log2FC, na.rm = T))%>%as.numeric()##Enrichment的计算方法为mS除以所有激酶的log2FC均值的绝对值
  Mean.FC$m = kinase.list %>%as.numeric()#统计激酶对应的底物数量
  Mean.FC$z.score = (((Mean.FC$mS - mean(new$log2FC, na.rm = T)) *
                        sqrt(Mean.FC$m))/sd(new$log2FC, na.rm = T))%>%as.numeric()##z.score打分的计算方式为：激酶log2FC减所有激酶的log2FC均值乘以激酶底物数量平方根，再除以所有激酶的log2FC的标准差
  Mean.FC$p.value = pnorm(-abs(Mean.FC$z.score))%>%as.numeric()##p.value是计算z.score打分绝对值负值的正态分布
  Mean.FC$FDR = p.adjust(Mean.FC$p.value, method = "fdr")%>%as.numeric()
  savexlsx(KSData.dataset.abbrev, file = paste0(path,"/Kinase-Substrate Links.xlsx"))
  savexlsx(Mean.FC[order(Mean.FC$Kinase.Gene), -which(colnames(Mean.FC)=="mS")],file = paste0(path,"/KSEA Kinase Scores.xlsx"))
}
