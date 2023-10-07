#' @export
corrnetworksetInput <- function(id = "corrnetwork",
                               right = T) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    dropdownButton(
      inputId = ns("cornetworkdrop"),
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
          selectInput(inputId = ns("algorithm"), label = "点排布模式:", 
                      choices = c("circle", "star",
                                  "nicely", "dh", "gem",
                                  "graphopt", "grid", "mds",
                                  "sphere", "randomly", "fr",
                                  "kk", "drl", "lgl"),
                      selected = "circle"),
          colourpicker::colourInput(ns("col1"), "选择正相关颜色:", "firebrick3", returnName = T),
          colourpicker::colourInput(ns("col2"), "选择负相关颜色:", "navy", returnName = T),
          sliderInput(inputId = ns("linewidth"), label = "连线粗细:", min = 0.1, max = 5, step = 0.1, value = c(0.2, 1)),
          sliderInput(inputId = ns("linealpha"), label = "连线透明度:", min = 0.1, max = 1, step = 0.1, value = 0.6)
        ),
        tabPanel(
          "运算",
          br(),
          sliderInput(inputId = ns("corfilter"), label = "相关性筛选:", min = 0.5, max = 0.99, step = 0.01, value = 0.95),
          sliderTextInput(inputId = ns("pfilter"), label = "P-value筛选:", choices = c(0.05, 0.01, 0.001), selected = 0.05)
        ),
        tabPanel(
          "图形",
          br(),
          p("此页面参数不设置将根据数据自动变化", style = "color:red"),
          numericInput(inputId = ns("pointsize"), label = "节点大小:", value = NA),
          numericInput(inputId = ns("textsize"), label = "节点标签大小:", value = NA),
          sliderInput(ns("width"), "图像宽度:", min = 5, max = 20, value = 10, step = 1),
          sliderInput(ns("height"), "图像高度:", min = 5, max = 20, value = 10, step = 1)
        )
      )
    )
  )
}


#' @export
corrnetworksetserver <- function(input, output, session, data) {
  parameter <- reactiveValues()
  
  observe({
    req(data$data)
    parameter$algorithm <- input$algorithm
    parameter$color <- c(
      "positive correlation" = input$col1,
      "negative correlation" = input$col2
    )
    parameter$linewidth <- input$linewidth
    parameter$linealpha <- input$linealpha
    parameter$corfilter <- input$corfilter
    parameter$pfilter <- input$pfilter
    parameter$pointsize <- input$pointsize
    parameter$textsize <- input$textsize
    
    parameter$width <- input$width
    parameter$height <- input$height
  })
  
  return(parameter)
}

#' @export
corrnetworkplotserver <- function(input, output, session, data, parameter) {
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
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(filename = data$data,imagetype = NA,savecorr = F),args)
        
        returndata <- do.call(what = map_common_corrnetwork,args = args)
      },
      width = parameter$width * 72,
      height = parameter$height * 72
    )
  })
  
  observe({
    output$downloadplot <- downloadHandler(
      filename = "相关性网络图.zip",
      content = function(file) {
        req(data$data)
        req(parameter)
        
        wd <- getwd()
        setwd(tempdir())
        
        name <- paste0("相关性网络图-", RandomCode(5))
        
        setwddir(filename = name, force = T)
        
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(filename = data$data,
                       imagetype = input$imagetype,
                       mapname = "corrnetwork"),args)
        
        returndata <- do.call(what = map_common_corrnetwork,args = args)
        
        setwd("../")
        zip::zip(zipfile = file, files = name)
        setwd(wd)
      }
    )
  })
  
  observe({
    output$alldownloadData <- downloadHandler(
      filename = "相关性网络图.zip",
      content = function(file) {
        req(data$data)
        req(parameter)
        
        wd <- getwd()
        setwd(tempdir())
        
        name <- paste0("相关性网络图-", RandomCode(5))
        
        setwddir(filename = name, force = T)
        
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(filename = data$path,
                       imagetype = input$imagetype),args)
        
        returndata <- do.call(what = map_common_corrnetwork,args = args)
        
        setwd("../")
        
        zip::zip(zipfile = file, files = name)
        
        setwd(wd)
      }
    )
  })
}

#' @export
corrnetworkui <- function(id = "corrnetwork") {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        p("相关性网络图", style = "color:black;font-size:30px;font-weight:bold;"),
        hr(), readdataInput(id = ns("xlsx")),
        hr(), allplotdownloadui(id = ns("plot")),
        width = 3
      ),
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "相关性网络图展示",
            br(),
            fixedRow(column(width = 1, offset = 11, corrnetworksetInput(id = ns("set")))),
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
corrnetworkserver <- function(input, output, session) {
  data <- callModule(readdataserver, "xlsx",row.names = T)
  parameter <- callModule(corrnetworksetserver, "set", data = data)
  callModule(corrnetworkplotserver, "plot", data = data, parameter = parameter)
}


#' @export
app_corrnetwork <- function() {
  suppressMessages(library(shiny))
  suppressMessages(library(shinydashboard))
  suppressMessages(library(shinyWidgets))
  suppressMessages(library(shinyjs))
  suppressMessages(library(DT))
  suppressMessages(library(colourpicker))
  
  shinyApp(
    ui = bootstrapPage(
      corrnetworkui(id = "corrnetwork")
    ),
    server = function(input, output, session) {
      callModule(corrnetworkserver, "corrnetwork")
    }
  )
}
