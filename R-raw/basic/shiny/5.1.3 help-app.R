#!/opt/conda/bin/Rscript

#' @export
helpdoc <- list("项目流程"=list("空间代谢组" = "project_SpatialMetabolomeAnalysis.html",
                            "质谱流式" = "project_CytofAnalysis.html",
                            "中药成分鉴定分析" = "TcmReport.html",
                            "中药入血成分分析" = "BloodReport.html"),
                "分析流程"=list("蛋白与磷酸化联合分析" = "flow_capp_Analyze.html",
                            "蛋白数据预处理" = "profetch.html",
                            "多元统计分析" = "flow_mulstatistics_Analyze.html",
                            "热敏蛋白分析" = "flow_TPP_Analyze.html",
                            "代谢数据整理分析" = "flow_datafetch_meta.html",
                            "预处理分析" = "flow_predealdata_Analyze.html",
                            "批次校正流程"= "Batch_correct.html",
                            "分子对接分析" = "flow_autodock_Analyze.html",
                            "空代数据提取" = "flow_spatialmetabolism_getdata.html",
                            "空代共定位分析" = "spatial_colocalization.html",
                            "空代拟时序分析" = "flow_spatialmetabolism_pseudotime.html"),
                "绘图工具"=list("热图" = "map_common_heatmap.html",
                            "相关性图" = "map_common_corrplot.html",
                            "相关性网络图" = "map_common_corrnetwork.html",
                            "ROC及逻辑回归ROC" = "map_common_roc.html",
                            "Zscore" = "map_common_zscore.html",
                            "脂质相关绘图" = "Lipid_unique_plot.html",
                            "空代成像图及质谱图" = "map_spatialmetabolism_imzmlimage.html"),
                "数据库搭建"=list("保留时间预测" = "tool_gnn_rt.html",
                             "EMDB2.0数据库搭建流程" = "database_EMDB2.0.html"),
                "分析工具"=list("开发模式" = "tool_developmode.html",
                            "颜色选择" = "tool_colormap.html",
                            "调用路径" = "tool_path.html",
                            "数据读取与保存" = "tool_ReadWriteData.html",
                            "数据库版本选择" = "tool_selectfile.html",
                            "保留时间预测" = "tool_gnn_rt.html"))

#' @export
helpui <- function(id = "help") {
  ns <- NS(id)
  tagList(
    useShinyjs(),
    fixedRow(column(width = 1, offset = 10,    
                    dropdownButton(
                      inputId = ns("drop"),
                      label = "帮助文档",
                      icon = icon("magnifying-glass"),
                      status = "primary",
                      circle = FALSE,
                      right = T,
                      selectInput(inputId = ns("htmlpath"), 
                                  label = "选择查看的帮助文档:", 
                                  choices = helpdoc)
                    ))),
    uiOutput(outputId = ns("htmlshow"),inline = T)
  )
}

#' @export
helpserver <- function(input, output, session) {
  
  observe({
    output$htmlshow <- renderUI({
      # includeHTML(paste0("/opt/conda/lib/R/library/lmbio/doc/",input$htmlpath))
      tags$iframe(src = paste0("lmbio/",input$htmlpath),width = "100%",height = 800,style="border:none")
    })
  })
  
}

#' @export
app_help <- function() {
  suppressMessages(library(shiny))
  suppressMessages(library(shinyWidgets))
  suppressMessages(library(shinyjs))
  addResourcePath("lmbio","/opt/conda/lib/R/library/lmbio/doc/")
  
  shinyApp(
    ui = bootstrapPage(
      helpui(id = "help")
    ),
    server = function(input, output, session) {
      callModule(helpserver, "help")
    }
  )
  
}
