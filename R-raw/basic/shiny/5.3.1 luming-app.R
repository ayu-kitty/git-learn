#' @export
lumingui <- function() {
  navbarPage(
    title = p("LumingOS", style = "font-size:40px;font-weight:bold;color:#5FA082"),
    navbarMenu(
      title = HTML('<font style="font-size:20px;font-weight:bold">项目分析</font>'),
      tabPanel(title = HTML('<font style="font-size:18px;font-weight:bold">代谢项目</font>'), metaanalystui(id = "meta"))
    ),
    navbarMenu(
      title = HTML('<font style="font-size:20px;font-weight:bold">分析流程</font>'),
      tabPanel(title = HTML('<font style="font-size:18px;font-weight:bold">趋势/时间序列分析</font>'), stemui(id = "stem")),
      tabPanel(title = HTML('<font style="font-size:18px;font-weight:bold">随机森林分析</font>'), randomForestui(id = "randomForest"))
    ),
    navbarMenu(
      title = HTML('<font style="font-size:20px;font-weight:bold">可视化</font>'),
      tabPanel(title = HTML('<font style="font-size:18px;font-weight:bold">热图</font>'), heatmapui(id = "heatmap")),
      tabPanel(title = HTML('<font style="font-size:18px;font-weight:bold">火山图</font>'), volcanoui(id = "volcano")),
      tabPanel(title = HTML('<font style="font-size:18px;font-weight:bold">相关性图</font>'), corrplotui(id = "corrplot")),
      tabPanel(title = HTML('<font style="font-size:18px;font-weight:bold">相关性网络图</font>'), corrnetworkui(id = "corrnetwork")),
      tabPanel(title = HTML('<font style="font-size:18px;font-weight:bold">韦恩图</font>'), vennui(id = "venn")),
      tabPanel(title = HTML('<font style="font-size:18px;font-weight:bold">ROC曲线</font>'), rocui(id = "roc")),
      tabPanel(title = HTML('<font style="font-size:18px;font-weight:bold">逻辑回归之ROC曲线</font>'), logrocui(id = "logroc"))
    ),
    navbarMenu(
      title = HTML('<font style="font-size:20px;font-weight:bold">小工具</font>'),
      tabPanel(title = HTML('<font style="font-size:18px;font-weight:bold">代谢物信息查询</font>'), compoundidui(id = "compoundid")) # ,
      # tabPanel(title = HTML('<font style="font-size:18px;font-weight:bold">中药数据库查询</font>'), tcmui(id = "tcm"))
    ),
    tabPanel(title = HTML('<font style="font-size:20px;font-weight:bold">帮助</font>'), helpui(id = "help")),
    tabPanel(title = HTML('<font style="font-size:20px;font-weight:bold">关于</font>'), infoui(id = "info")),
    inverse = F,
    collapsible = T,
    windowTitle = "LumingOS",
    id = "all"
  )
}

#' @export
lumingserver <- function(input, output, session) {
  
  # 项目分析
  callModule(metaanalystserver,"meta")
  
  # 分析流程
  callModule(stemserver, "stem")
  callModule(randomForestserver, "randomForest")
  
  # 可视化
  callModule(heatmapserver, "heatmap")
  callModule(volcanoserver, "volcano")
  callModule(corrplotserver, "corrplot")
  callModule(corrnetworkserver, "corrnetwork")
  callModule(vennserver, "venn")
  callModule(rocserver, "roc")
  callModule(logrocserver, "logroc")

  # 小工具
  callModule(compoundidserver, "compoundid")
  # callModule(tcmserver, "tcm")
  
  # 帮助
  callModule(helpserver, "help")
  
  # 公司信息
  callModule(infoserver, "info")
  
}


#' 云端分析
#'
#' @export
app_luming <- function(){
  suppressMessages(library(shiny))
  suppressMessages(library(shinydashboard))
  suppressMessages(library(shinyWidgets))
  suppressMessages(library(shinyjs))
  suppressMessages(library(DT))
  suppressMessages(library(colourpicker))
  addResourcePath("lmbio","/opt/conda/lib/R/library/lmbio/doc/")
  shinyApp(
    ui = lumingui(),
    server = lumingserver,
  )
}