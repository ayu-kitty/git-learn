#!/opt/conda/bin/Rscript

#' @export
heatmapsetInput <- function(id = "heatmap",
                            right = T) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    dropdownButton(
      inputId = ns("heatmapdrop"),
      label = "设置",
      icon = icon("gear"),
      status = "primary",
      circle = FALSE,
      right = right,
      tabsetPanel(
        type = "tabs",
        tabPanel(
          "展示",
          br(),
          textInput(inputId = ns("main"), label = "标题名称:", value = "Heatmap"),
          colourpicker::colourInput(ns("col3"), "选择上调颜色:", "firebrick3", returnName = T),
          colourpicker::colourInput(ns("col2"), "选择过度颜色:", "white", returnName = T),
          colourpicker::colourInput(ns("col1"), "选择下调颜色:", "navy", returnName = T),
          selectInput(inputId = ns("angle_col"), label = "列名展示角度:", choices = c("270", "0", "45", "90", "315")),
          materialSwitch(inputId = ns("show_rownames"), label = "是否显示行名", value = TRUE, status = "primary"),
          materialSwitch(inputId = ns("show_colnames"), label = "是否显示列名", value = TRUE, status = "primary"),
          materialSwitch(inputId = ns("legend"), label = "是否显示图例", value = TRUE, status = "primary"),
          materialSwitch(inputId = ns("border_color"), label = "是否显示边框", value = FALSE, status = "primary"),
          colourpicker::colourInput(ns("col4"), "边框颜色:", "black", returnName = T)
        ),
        tabPanel(
          "运算",
          br(),
          materialSwitch(inputId = ns("log"), label = "是否log归一化", value = FALSE, status = "primary"),
          selectInput(inputId = ns("scale"), label = "标准化模式:", choices = c("row", "column", "none")),
          hr(),
          materialSwitch(inputId = ns("cluster_rows"), label = "是否横向聚类", value = TRUE, status = "primary"),
          sliderInput(inputId = ns("cutree_rows"), label = "横聚类分割:", value = 1, min = 1, max = 10, step = 1),
          hr(),
          materialSwitch(inputId = ns("cluster_cols"), label = "是否纵向聚类", value = FALSE, status = "primary"),
          sliderInput(inputId = ns("cutree_cols"), label = "纵聚类分割:", value = 1, min = 1, max = 10, step = 1),
          hr(),
          sliderInput(inputId = ns("rowgroup"), label = "横向注释标签数量:", value = 0, min = 0, max = 3, step = 1),
          sliderInput(inputId = ns("colgroup"), label = "纵向注释标签数量:", value = 0, min = 0, max = 3, step = 1)
        ),
        tabPanel(
          "图形",
          br(),
          p("此页面参数不设置将根据数据自动变化", style = "color:red"),
          numericInput(inputId = ns("cellheight"), label = "单元格高度:", value = NA),
          numericInput(inputId = ns("cellwidth"), label = "单元格宽度:", value = NA),
          numericInput(inputId = ns("height"), label = "图片高度:", value = NA),
          numericInput(inputId = ns("width"), label = "图片宽度:", value = NA),
          numericInput(inputId = ns("treeheight_col"), label = "纵聚类树高度:", value = NA),
          numericInput(inputId = ns("treeheight_row"), label = "横聚类树宽度:", value = NA),
          numericInput(inputId = ns("fontsize_row"), label = "行字体大小:", value = NA),
          numericInput(inputId = ns("fontsize_col"), label = "列字体大小:", value = NA)
        )
      )
    )
  )
}



#' @export
heatmapsetserver <- function(input, output, session, data) {
  parameter <- reactiveValues()
  
  observe({
    if (input$border_color) {
      shinyjs::show(id = "col4")
    } else {
      shinyjs::hide(id = "col4")
    }
    
    if (input$cluster_rows) {
      shinyjs::show(id = "cutree_rows")
    } else {
      shinyjs::hide(id = "cutree_rows")
    }
    
    if (input$cluster_cols) {
      shinyjs::show(id = "cutree_cols")
    } else {
      shinyjs::hide(id = "cutree_cols")
    }
  })
  
  observe({
    req(data$data)
    parameter$main <- input$main
    parameter$color <- colorRampPalette(c(input$col1, input$col2, input$col3))(1000)
    parameter$angle_col <- input$angle_col
    
    parameter$show_rownames <- input$show_rownames
    parameter$show_colnames <- input$show_colnames
    parameter$legend <- input$legend
    parameter$border_color <- ifelse(input$border_color, input$col4, NA)
    parameter$log <- input$log
    parameter$scale <- input$scale
    parameter$cluster_rows <- input$cluster_rows
    parameter$cutree_rows <- input$cutree_rows
    parameter$cluster_cols <- input$cluster_cols
    parameter$cutree_cols <- input$cutree_cols
    parameter$rowgroup <- input$rowgroup
    parameter$colgroup <- input$colgroup
    
    parameter$treeheight_row <- input$treeheight_row
    parameter$treeheight_col <- input$treeheight_col
    parameter$cellheight <- input$cellheight
    parameter$cellwidth <- input$cellwidth
    parameter$fontsize_row <- input$fontsize_row
    parameter$fontsize_col <- input$fontsize_col
    parameter$height <- input$height
    parameter$width <- input$width
  })
  
  return(parameter)
}

#' @export
heatmapplotserver <- function(input, output, session, data, parameter) {
  output$plot <- renderPlot(
    {},
    width = 6 * 72,
    height = 6 * 72
  )
  
  maprange <- reactiveValues()
  
  observe({
    req(data$data)
    req(parameter)
    
    output$plot <- renderPlot(
      {
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(data = data$data,imagetype = NA),args)
        
        returndata <- do.call(what = auto_heatmap,args = args)
        
        maprange$height <- returndata$height
        maprange$width <- returndata$width
      },
      height = if(is.null(maprange$height)){10 * 72}else{maprange$height * 72},
      width = if(is.null(maprange$width)){10 * 72}else{maprange$width * 72}
    )
  })
  
  observe({
    output$downloadplot <- downloadHandler(
      
      filename = function() {paste("heatmap", input$imagetype, sep = ".")},
      content = function(file) {
        req(data$data)
        req(parameter)
        
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(data = data$data,
                       imagetype = input$imagetype,
                       mapname = "heatmap"),args)
        
        returndata <- do.call(what = auto_heatmap,args = args)
      }
    )
  })
  
  observe({
    output$alldownloadData <- downloadHandler(
      filename = "热图.zip",
      content = function(file) {
        req(data$data)
        req(parameter)
        
        wd <- getwd()
        setwd(tempdir())
        
        name <- paste0("热图-", RandomCode(5))
        
        setwddir(filename = name, force = T)
        
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(filename = data$path,
                       imagetype = input$imagetype),args)
        
        returndata <- do.call(what = map_common_heatmap,args = args)
        
        setwd("../")
        zip::zip(zipfile = file, files = name)
        
        setwd(wd)
      }
    )
  })
}

#' @export
heatmapui <- function(id = "heatmap") {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        p("热图", style = "color:black;font-size:30px;font-weight:bold;"),
        hr(), readdataInput(id = ns("xlsx")),
        hr(), allplotdownloadui(id = ns("plot")),
        width = 3
      ),
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "热图展示",
            br(),
            fixedRow(column(width = 1, offset = 11, heatmapsetInput(id = ns("set")))),
            plotui(id = ns("plot"))
          ),
          tabPanel(
            "数据展示",
            br(),
            readdataOutput(id = ns("xlsx"))
          )
        )
      )
    )
  )
}

#' @export
heatmapserver <- function(input, output, session) {
  
  data <- callModule(readdataserver, "xlsx",row.names = T)
  parameter <- callModule(heatmapsetserver, "set", data = data)
  callModule(heatmapplotserver, "plot", data = data, parameter = parameter)
  
}

#' @export
app_heatmap <- function() {
  suppressMessages(library(shiny))
  suppressMessages(library(shinydashboard))
  suppressMessages(library(shinyWidgets))
  suppressMessages(library(shinyjs))
  suppressMessages(library(DT))
  suppressMessages(library(colourpicker))
  
  shinyApp(
    ui = bootstrapPage(
      heatmapui(id = "heatmap")
    ),
    server = function(input, output, session) {
      callModule(heatmapserver, "heatmap")
    }
  )
  
}
