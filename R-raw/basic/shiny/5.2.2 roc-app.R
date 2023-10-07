#!/opt/conda/bin/Rscript

#' @export
rocsetInput <- function(id = "roc",
                        right = T) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    dropdownButton(
      inputId = ns("rocdrop"),
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
          checkboxInput(ns("AUC"), "是否展示AUC", TRUE),
          sliderInput(ns("cex"), "AUC标签大小:", min = 1, max = 20, value = 4, step = 1),
          materialSwitch(inputId = ns("hide1"), label = "是否颜色修改", value = FALSE, status = "primary"),
          colourpicker::colourInput(ns("col1"), "选择物质1颜色:", "green4", returnName = T),
          colourpicker::colourInput(ns("col2"), "选择物质2颜色:", "blue3", returnName = T),
          colourpicker::colourInput(ns("col3"), "选择物质3颜色:", "firebrick", returnName = T),
          colourpicker::colourInput(ns("col4"), "选择物质4颜色:", "gold", returnName = T),
          colourpicker::colourInput(ns("col5"), "选择物质5颜色:", "darkviolet", returnName = T)
        ),
        tabPanel(
          "运算",
          br(),
          checkboxInput(ns("right"), "自动AUC>0.5", TRUE)
        ),
        tabPanel(
          "图形",
          br(),
          sliderInput(ns("nudge_x"), "AUC标签x轴偏移:", min = -1, max = 1, value = 0, step = 0.1),
          sliderInput(ns("nudge_y"), "AUC标签y轴偏移:", min = -1, max = 1, value = 0, step = 0.1),
          sliderInput(ns("width"), "图像宽度:", min = 5, max = 20, value = 5, step = 1),
          sliderInput(ns("height"), "图像高度:", min = 5, max = 20, value = 5, step = 1),
          sliderTextInput(inputId = ns("aspect.ratio"), label = "图像长宽比:", choices = c(1:10 / 10, 2:10), grid = TRUE, selected = 1)
        )
      )
    )
  )
}


#' @export
rocsetserver <- function(input, output, session, data) {
  
  parameter <- reactiveValues()
  
  observe({
    if (input$hide1) {
      shinyjs::show(id = "col1")
      shinyjs::show(id = "col2")
      shinyjs::show(id = "col3")
      shinyjs::show(id = "col4")
      shinyjs::show(id = "col5")
    } else {
      shinyjs::hide(id = "col1")
      shinyjs::hide(id = "col2")
      shinyjs::hide(id = "col3")
      shinyjs::hide(id = "col4")
      shinyjs::hide(id = "col5")
    }
  })
  
  observe({
    req(data$data)
    parameter$AUC <- input$AUC
    parameter$cex <- input$cex
    parameter$collist <- c(input$col1, input$col2, input$col3, input$col4, input$col5)
    parameter$right <- input$right
    parameter$nudge_x <- input$nudge_x
    parameter$nudge_y <- input$nudge_y
    parameter$width <- input$width
    parameter$height <- input$height
    parameter$aspect.ratio <- input$aspect.ratio
  })
  
  return(parameter)
}

#' @export
rocplotserver <- function(input, output, session, data, parameter) {
  output$plot <- renderPlot(
    {},
    width = 6 * 72,
    height = 6 * 72
  )
  
  observe({
    req(data$data)
    req(parameter)
    
    output$plot <- renderPlot(
      {
        number <- dim(data$data)[1]-1
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(data = data$data,
                       imagetype = NA,
                       number = number,
                       saveroc = F),args)
        
        returndata <- do.call(what = auto_roc,args = args)
      },
      width = parameter$width * 72,
      height = parameter$height * 72
    )
  })
  
  observe({
    output$downloadplot <- downloadHandler(
      filename = "roc.zip",
      content = function(file) {
        req(data$data)
        req(parameter)
        
        number <- dim(data$data)[1]-1
        
        wd <- getwd()
        setwd(tempdir())
        
        name <- paste0("ROC-", RandomCode(5))
        
        setwddir(filename = name, force = T)
        
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(data = data$data,
                       imagetype = input$imagetype,
                       number = number),args)
        
        returndata <- do.call(what = auto_roc,args = args)
        
        setwd("../")
        
        zip::zip(zipfile = file, files = name)
        setwd(wd)
      }
    )
  })
}

#' @export
rocui <- function(id = "roc") {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        p("ROC曲线", style = "color:black;font-size:30px;font-weight:bold;"),
        hr(), readdataInput2(id = ns("xlsx")),
        hr(), plotdownloadui(id = ns("plot")),
        width = 3
      ),
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "ROC图展示",
            br(),
            fixedRow(column(width = 1, offset = 11, rocsetInput(id = ns("set")))),
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
rocserver <- function(input, output, session) {
  data <- callModule(readdataserver2, "xlsx",)
  parameter <- callModule(rocsetserver, "set", data = data)
  callModule(rocplotserver, "plot", data = data, parameter = parameter)
}


#' @export
app_roc <- function() {
  suppressMessages(library(shiny))
  suppressMessages(library(shinydashboard))
  suppressMessages(library(shinyWidgets))
  suppressMessages(library(shinyjs))
  suppressMessages(library(DT))
  suppressMessages(library(colourpicker))
  
  shinyApp(
    ui = bootstrapPage(
      rocui(id = "roc")
    ),
    server = function(input, output, session) {
      callModule(rocserver, "roc")
    }
  )
}
