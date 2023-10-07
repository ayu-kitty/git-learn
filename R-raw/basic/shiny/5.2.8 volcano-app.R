#' @export
volcanosetInput <- function(id = "volcano",
                            right = T) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    dropdownButton(
      inputId = ns("volcanodrop"),
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
          numericInput(inputId = ns("cex"), label = "点大小:", value = 1),
          sliderInput(inputId = ns("size"), label = "点大小范围:", min = 0.1, max = 5, step = 0.1, value = c(0.5, 3)),
          colourpicker::colourInput(ns("col1"), "选择下调颜色:", "blue", returnName = T),
          colourpicker::colourInput(ns("col2"), "选择非显著颜色:", "grey", returnName = T),
          colourpicker::colourInput(ns("col3"), "选择上调颜色:", "red", returnName = T),
          materialSwitch(inputId = ns("adjustoutx"), label = "超出X轴范围的点是否在轴上显示", value = FALSE, status = "primary"),
          materialSwitch(inputId = ns("adjustouty"), label = "超出Y轴范围的点是否在轴上显示", value = FALSE, status = "primary")
        ),
        tabPanel(
          "运算",
          br(),
          numericInput(inputId = ns("drawFC"), label = "FC筛选标准:", value = NA),
          numericInput(inputId = ns("VIP"), label = "VIP筛选标准:", value = NA),
          numericInput(inputId = ns("p"), label = "P-value筛选标准:", value = 0.05)
        ),
        tabPanel(
          "图形",
          br(),
          p("此页面参数不设置将根据数据自动变化", style = "color:red"),
          numericInput(inputId = ns("x"), label = "x轴范围:", value = NA),
          numericInput(inputId = ns("y"), label = "y轴范围:", value = NA),
          sliderInput(ns("width"), "图像宽度:", min = 5, max = 20, value = 7, step = 1),
          sliderInput(ns("height"), "图像高度:", min = 5, max = 20, value = 7, step = 1),
          sliderTextInput(inputId = ns("aspect.ratio"), label = "图像长宽比:", choices = c(1:10 / 10, 2:10), grid = TRUE, selected = 1)
        )
      )
    )
  )
}


#' @export
volcanosetserver <- function(input, output, session, data) {
  parameter <- reactiveValues()
  
  observe({
    if (is.na(input$VIP)) {
      shinyjs::hide(id = "size")
      shinyjs::show(id = "cex")
    } else {
      shinyjs::show(id = "size")
      shinyjs::hide(id = "cex")
    }
  })
  
  observe({
    req(data$data)
    parameter$cex <- input$cex
    parameter$size <- input$size
    parameter$color <- c(input$col1, input$col2, input$col3)
    parameter$drawFC <- input$drawFC
    parameter$VIP <- input$VIP
    parameter$p <- input$p
    parameter$x <- input$x
    parameter$y <- input$y
    parameter$width <- input$width
    parameter$height <- input$height
    parameter$aspect.ratio <- input$aspect.ratio
    parameter$adjustoutx <- input$adjustoutx
    parameter$adjustouty <- input$adjustouty
  })
  
  return(parameter)
}

#' @export
volcanoplotserver <- function(input, output, session, data, parameter) {
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
        args <- c(list(data = data$data,imagetype = NA),args)
        
        returndata <- do.call(what = auto_volcano,args = args)
      },
      width = parameter$width * 72,
      height = parameter$height * 72
    )
  })
  
  observe({
    output$downloadplot <- downloadHandler(
      filename = function() {
        paste("volcano", input$imagetype, sep = ".")
      },
      content = function(file) {
        req(data$data)
        req(parameter)
        
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(data = data$data,
                       imagetype = input$imagetype,
                       mapname = "volcano"),args)
        
        returndata <- do.call(what = auto_volcano,args = args)
        
      }
    )
  })
  
  observe({
    output$alldownloadData <- downloadHandler(
      filename = "火山图.zip",
      content = function(file) {
        req(data$data)
        req(parameter)
        
        wd <- getwd()
        setwd(tempdir())
        
        name <- paste0("火山图-", RandomCode(5))
        
        setwddir(filename = name, force = T)
        
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(filename = data$path,
                       imagetype = input$imagetype),args)
        
        returndata <- do.call(what = map_common_volcano,args = args)
        
        setwd("../")
        zip::zip(zipfile = file, files = name)
        
        setwd(wd)
      }
    )
  })
}

#' @export
volcanoui <- function(id = "volcano") {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        p("火山图", style = "color:black;font-size:30px;font-weight:bold;"),
        hr(), readdataInput(id = ns("xlsx")),
        hr(), allplotdownloadui(id = ns("plot")),
        width = 3
      ),
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "火山图展示",
            br(),
            fixedRow(column(width = 1, offset = 11, volcanosetInput(id = ns("set")))),
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
volcanoserver <- function(input, output, session) {
  data <- callModule(readdataserver, "xlsx")
  parameter <- callModule(volcanosetserver, "set", data = data)
  callModule(volcanoplotserver, "plot", data = data, parameter = parameter)
}


#' @export
app_volcano <- function() {
  suppressMessages(library(shiny))
  suppressMessages(library(shinydashboard))
  suppressMessages(library(shinyWidgets))
  suppressMessages(library(shinyjs))
  suppressMessages(library(DT))
  suppressMessages(library(colourpicker))
  
  shinyApp(
    ui = bootstrapPage(
      volcanoui(id = "volcano")
    ),
    server = function(input, output, session) {
      callModule(volcanoserver, "volcano")
    }
  )
}
