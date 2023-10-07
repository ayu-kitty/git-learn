#!/opt/conda/bin/Rscript

#' Mulgetpseudotime
#'
#' 批量对拟时序数据进行数据获取
#'
#' @param moderange 正负离子模式
#' @param rdspath rds数据路径
#' @param annosample 是否根据样本出图
#' @param savepath 保存路径
#' @param asp 成像图比例
#' @param ... 见[Mulgetrdsdata()]
#'
#' @export
Mulgetpseudotime <- function(rdspath = "./sample/cluster/pseudotime/",
                             savepath = "./sample/map/pseudotime/",
                             moderange = c("neg", "pos"),
                             annosample = F,
                             asp = 1,
                             ...){
  # 获取拟时序数据
  Mulgetrdsdata(rdspath = rdspath,
                moderange = moderange,
                anamoudle = getpseudotime_data,
                savepath = savepath,
                ...)

  # mz筛选图
  Mulgetrdsdata(rdspath = rdspath,
                moderange = moderange,
                anamoudle = getpseudotime_plotordering,
                savepath = savepath,
                ...)

  # sscc聚类的拟时序图
  Mulgetrdsdata(rdspath = rdspath,
                moderange = moderange,
                anamoudle = getpseudotime_Cluster,
                savepath = savepath,
                ...)

  # # sscc聚类及State分面的拟时序图
  # Mulgetrdsdata(rdspath = rdspath,
  #               moderange = moderange,
  #               anamoudle = getpseudotime_ClustertoState,
  #               savepath = savepath,
  #               ...)

  # 拟时序图
  Mulgetrdsdata(rdspath = rdspath,
                moderange = moderange,
                anamoudle = getpseudotime_Pseudotime,
                savepath = savepath,
                ...)

  # Cluster分面的拟时序图
  Mulgetrdsdata(rdspath = rdspath,
                moderange = moderange,
                anamoudle = getpseudotime_PseudotimetoCluster,
                savepath = savepath,
                ...)

  # # State的拟时序图
  # Mulgetrdsdata(rdspath = rdspath,
  #               moderange = moderange,
  #               anamoudle = getpseudotime_State,
  #               savepath = savepath,
  #               ...)

  # 拟时序密度图
  Mulgetrdsdata(rdspath = rdspath,
                moderange = moderange,
                anamoudle = getpseudotime_Clusterdensity,
                savepath = savepath,
                ...)

  # Cluster分面的拟时序密度图
  Mulgetrdsdata(rdspath = rdspath,
                moderange = moderange,
                anamoudle = getpseudotime_Clusterdensity2,
                savepath = savepath,
                ...)

  # # State分面的拟时序密度图
  # Mulgetrdsdata(rdspath = rdspath,
  #               moderange = moderange,
  #               anamoudle = getpseudotime_ClusterdensitytoState,
  #               savepath = savepath,
  #               ...)

  # 强度分布图
  Mulgetrdsdata(rdspath = rdspath,
                moderange = moderange,
                anamoudle = getpseudotime_Clusterexpress,
                savepath = savepath,
                ...)

  if(annosample){
    # sscc聚类及Sample分面的拟时序图
    Mulgetrdsdata(rdspath = rdspath,
                  moderange = moderange,
                  anamoudle = getpseudotime_ClustertoSample,
                  savepath = savepath,
                  ...)

    # Sample分面的拟时序图
    Mulgetrdsdata(rdspath = rdspath,
                  moderange = moderange,
                  anamoudle = getpseudotime_PseudotimetoSample,
                  savepath = savepath,
                  ...)

    # Sample分面的拟时序密度图
    Mulgetrdsdata(rdspath = rdspath,
                  moderange = moderange,
                  anamoudle = getpseudotime_ClusterdensitytoSample,
                  savepath = savepath,
                  ...)

    # Sample及Cluster分面的拟时序密度图
    Mulgetrdsdata(rdspath = rdspath,
                  moderange = moderange,
                  anamoudle = getpseudotime_ClusterdensitytoSample2,
                  savepath = savepath,
                  ...)

  }else{
    # 拟时序图成像图
    Mulgetrdsdata(rdspath = rdspath,
                  moderange = moderange,
                  anamoudle = getpseudotime_Pseudotimeheatmap,
                  savepath = savepath,
                  asp = asp,
                  ...)
  }

}

