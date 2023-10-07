#' @export
vennsetInput <- function(id = "venn",
                         right = T) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    dropdownButton(
      inputId = ns("venndrop"),
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
          selectInput(inputId = ns("fillcolor"), label = "选择呈现颜色:", choices = names(lmbio::palette_color)),
          checkboxInput(inputId = ns("col"), label = "是否呈现外框", value = FALSE),
          sliderInput(inputId = ns("cex"), label = "圈内数字大小:", min = 1, max = 5, step = 1, value = 2),
          sliderInput(inputId = ns("alpha"), label = "圈颜色透明度:", min = 0.1, max = 1, step = 0.1, value = 0.5),
          checkboxInput(inputId = ns("euler.d"), label = "是否启用欧拉图", value = TRUE),
          checkboxInput(inputId = ns("upset"), label = "是否Upset模式绘制", value = FALSE)
        ),
        tabPanel(
          "图形",
          br(),
          p("此页面参数不设置将根据数据自动变化", style = "color:red"),
          numericInput(inputId = ns("cat.cex"), label = "组别字体大小:", value = NA),
          sliderInput(ns("width"), "图像宽度:", min = 5, max = 20, value = 10, step = 1),
          sliderInput(ns("height"), "图像高度:", min = 5, max = 20, value = 10, step = 1)
        )
      )
    )
  )
}

#' @export
vennsetserver <- function(input, output, session, data) {
  parameter <- reactiveValues()
  
  observe({
    req(data$data)
    parameter$filllist <- SelectColors(palette = input$fillcolor,n = 5)
    parameter$collist <- if (input$col) {parameter$filllist} else {NA}
    parameter$cex <- input$cex
    parameter$alpha <- input$alpha
    parameter$euler.d <- input$euler.d
    parameter$scale <- input$euler.d
    parameter$upset <- input$upset
    parameter$cat.cex <- input$cat.cex
    parameter$width <- input$width
    parameter$height <- input$height
  })
  
  return(parameter)
}

#' @export
vennplotserver <- function(input, output, session, data, parameter) {
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
        args <- c(list(filename = data$data,
                       imagetype = NA,
                       savedata = F),args)
        
        returndata <- do.call(what = auto_vennmap,args = args)
      },
      width = parameter$width * 72,
      height = parameter$height * 72
    )
  })
  
  observe({
    output$downloadplot <- downloadHandler(
      filename = "韦恩图.zip",
      content = function(file) {
        req(data$data)
        req(parameter)
        
        wd <- getwd()
        setwd(tempdir())
        
        name <- paste0("韦恩图-", RandomCode(5))
        
        setwddir(filename = name, force = T)
        
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(filename = data$data,
                       imagetype = input$imagetype),args)
        
        returndata <- do.call(what = auto_vennmap,args = args)
        
        setwd("../")
        
        zip::zip(zipfile = file, files = name)
        setwd(wd)
      }
    )
  })
}

#' @export
vennui <- function(id = "venn") {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        p("韦恩图", style = "color:black;font-size:30px;font-weight:bold;"),
        hr(), readdataInput3(id = ns("xlsx")),
        hr(), plotdownloadui(id = ns("plot")),
        width = 3
      ),
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "韦恩图展示",
            br(),
            fixedRow(column(width = 1, offset = 11, vennsetInput(id = ns("set")))),
            plotui(id = ns("plot"))
          )
        )
      )
    )
  )
}

#' @export
vennserver <- function(input, output, session) {
  data <- callModule(readdataserver3, "xlsx")
  parameter <- callModule(vennsetserver, "set", data = data)
  callModule(vennplotserver, "plot", data = data, parameter = parameter)
}


#' @export
app_venn <- function() {
  suppressMessages(library(shiny))
  suppressMessages(library(shinydashboard))
  suppressMessages(library(shinyWidgets))
  suppressMessages(library(shinyjs))
  suppressMessages(library(DT))
  suppressMessages(library(colourpicker))
  
  shinyApp(
    ui = bootstrapPage(
      vennui(id = "venn")
    ),
    server = function(input, output, session) {
      callModule(vennserver, "venn")
    }
  )
}
