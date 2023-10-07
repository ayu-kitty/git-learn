#!/opt/conda/bin/Rscript

#' Mulgetrdsdata
#'
#' 批量对rds数据进行
#'
#' @param rdspath 数据路径
#' @param moderange 正负离子模式
#' @param anamoudle 分析模块
#' @param savepath 保存路径
#' @param ananame 分析名称
#' @param wantsample 想处理的样本名
#' @param ... 见anamoudle函数
#'
#' @export
Mulgetrdsdata <- function(rdspath = "./sample/cluster/sscc/",
                          moderange = c("neg", "pos"),
                          anamoudle = getssccfeaturetocor,
                          savepath = "./sample/map/sscc/",
                          ananame = NULL,
                          wantsample = NULL,
                          delsample = NULL,
                          ...) {
  for (mode in moderange) {

    # 获取.imzML文件绝对路径
    filename <- list.files(
      path = rdspath,
      pattern = paste0("-", mode, "-data.rds$"),
      full.names = F,
      recursive = T
    )

    if (length(filename) == 0) {
      print(paste0("在", rdspath, "目录下未找到", mode, "模式的rds文件"))
      next
    }

    # 获取样品名
    samplename <- gsub(
      pattern = paste0("-", mode, "-data.rds"),
      replacement = "",
      x = filename
    )

    if (!is.null(wantsample)) {
      samplename <- samplename[samplename %in% wantsample]
    }

    if (!is.null(delsample)) {
      samplename <- samplename[!(samplename %in% delsample)]
    }

    # 判断是否找到imzML文件
    if (length(samplename) > 0) {
      for (i in 1:length(samplename)) {
        anamoudle(
          samplename = samplename[i],
          savepath = ifelse(is.null(ananame),
            paste0(savepath, samplename[i], "/"),
            paste0(savepath, samplename[i], "/", ananame, "/")
          ),
          mode = mode,
          datapath = rdspath,
          ...
        )

        gc(verbose = T)
      }
    } else {
      print(paste0("在", rdspath, "目录下未找到", mode, "模式rds文件"))
    }
  }
}

#' Mulgetallrdsdata
#'
#' 批量对rds数据进行数据获取
#'
#' @param moderange 正负离子模式
#' @param rdspath rds数据路径
#' @param ssccrdspath sscc数据路径
#' @param tsnerdspath tsne数据路径
#' @param umaprdspath umap数据路径
#' @param savepath 保存路径
#' @param topFeaturesn 特征物质前几
#' @param asp 成像图长宽比
#'
#' @export
Mulgetallrdsdata <- function(rdspath = "./sample/cluster/",
                             ssccrdspath = paste0(rdspath, "sscc/"),
                             tsnerdspath = paste0(rdspath, "tsne/"),
                             umaprdspath = paste0(rdspath, "umap/"),
                             savepath = "./sample/map/sscc/",
                             moderange = c("neg", "pos"),
                             topFeaturesn = 10,
                             asp = 1,
                             ...) {
  # 相关性
  Mulgetrdsdata(
    rdspath = ssccrdspath,
    moderange = moderange,
    anamoudle = getssccfeaturetocor,
    savepath = savepath,
    ...
  )

  # 特征热图
  Mulgetrdsdata(
    rdspath = ssccrdspath,
    moderange = moderange,
    topFeaturesn = topFeaturesn,
    anamoudle = getsscctopfeaturetoheatmap,
    savepath = savepath,
    ...
  )

  # 聚类成像图
  Mulgetrdsdata(
    rdspath = ssccrdspath,
    moderange = moderange,
    anamoudle = getsscctoimage,
    savepath = savepath,
    asp = asp,
    ...
  )

  # 聚类质谱图
  Mulgetrdsdata(
    rdspath = ssccrdspath,
    moderange = moderange,
    anamoudle = getsscctoplot,
    savepath = savepath,
    ...
  )

  # 特征小提琴图
  Mulgetrdsdata(
    rdspath = ssccrdspath,
    moderange = moderange,
    topFeaturesn = topFeaturesn,
    anamoudle = getsscctopfeaturetoviolin,
    savepath = savepath,
    ananame = "特征物质-violin",
    ...
  )
  
  # tsne
  Mulgetrdsdata(
    rdspath = tsnerdspath,
    moderange = moderange,
    anamoudle = gettsnemap,
    savepath = savepath,
    ananame = "tsne",
    ssccpath = ssccrdspath,
    ...
  )

  # umap
  Mulgetrdsdata(
    rdspath = umaprdspath,
    moderange = moderange,
    anamoudle = getumapmap,
    savepath = savepath,
    ananame = "umap",
    ssccpath = ssccrdspath,
    ...
  )
}