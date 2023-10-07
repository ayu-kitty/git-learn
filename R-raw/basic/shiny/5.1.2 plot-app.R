#!/opt/conda/bin/Rscript

#' 图片下载
#' 
#' @param id ui编号 
#' @param downloadlabel 下载图标名称 
#' @param choices 选择项
#'
#' @export
plotdownloadui <- function(id = "plot",
                           downloadlabel = "图片保存",
                           choices = c("pdf", "png", "jpg")) {
  ns <- NS(id)
  tagList(
    radioButtons(inputId = ns("imagetype"), label = "图片格式", choices = choices, selected = "pdf"),
    downloadButton(ns("downloadplot"), downloadlabel)
  )
}

#' 图片下载
#' 
#' @param id ui编号 
#' @param downloadlabel 下载图标名称 
#' @param choices 选择项
#' @param right 图标位置
#'
#' @export
plotdownloadui2 <- function(id = "plot", 
                            downloadlabel = "图片保存",
                            right = T,
                            choices = c("pdf", "png", "jpg")) {
  ns <- NS(id)
  tagList(
    useShinyjs(),
    dropdownButton(
      inputId = ns("sumchartdrop"),
      label = "下载",
      icon = icon("gear"),
      status = "primary",
      circle = FALSE,
      right = right,
      radioButtons(inputId = ns("imagetype"), label = "图片格式", choices = choices, selected = "pdf"),
      downloadButton(ns("downloadplot"), downloadlabel)
    )
  )
}


#' 图片批量下载
#' 
#' @param id ui编号 
#' @param downloadlabel 下载图标名称 
#' @param choices 选择项
#'
#' @export
allplotdownloadui <- function(id = "plot", 
                              downloadlabel = "图片保存",
                              alldownloadlabel = "批量保存",
                              choices = c("pdf", "png", "jpg")) {
  ns <- NS(id)
  tagList(
    radioButtons(inputId = ns("imagetype"), label = "图片格式", choices = choices, selected = "pdf"),
    downloadButton(ns("downloadplot"), downloadlabel),
    downloadButton(ns("alldownloadData"), alldownloadlabel)
  )
}


#' 下载
#' 
#' @param id ui编号 
#' @param downloadlabel 下载图标名称 
#'
#' @export
downloadui <- function(id = "plot", downloadlabel = "结果保存") {
  ns <- NS(id)
  tagList(
    downloadButton(ns("downloadplot"), downloadlabel)
  )
}

#' 图片展示
#' 
#' @param id ui编号
#' 
#' @export
plotui <- function(id = "plot") {
  ns <- NS(id)
  tagList(
    shinycssloaders::withSpinner(
      plotOutput(outputId = ns("plot"), inline = T),
      type = 4
    )
  )
}
