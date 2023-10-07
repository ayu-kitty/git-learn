#' @export
stemsetInput <- function(id = "stem",
                         right = T) {
  ns <- NS(id)
  tagList(
    useShinyjs(),
    dropdownButton(
      inputId = ns("stemdrop"),
      label = "设置",
      icon = icon("gear"),
      status = "primary",
      circle = FALSE,
      right = right,
      sliderInput(inputId = ns("cluster"), label = "聚类数量:", min = 3, max = 30, step = 1, value = 10),
      selectInput(inputId = ns("standardise"), label = "数据标准化:", choices = c("Z-score", "Log2"), selected = "Z-score"),
      materialSwitch(inputId = ns("filter_std"), label = "是否标准差筛选", value = FALSE, status = "primary"),
      sliderInput(inputId = ns("min.std"), label = "标准差筛选标准:", min = 0, max = 0.5, step = 0.01, value = 0)
    )
  )
}

#' @export
stemsetserver <- function(input, output, session, data) {
  parameter <- reactiveValues()
  
  observe({
    if (input$filter_std) {
      shinyjs::show(id = "min.std")
    } else {
      shinyjs::hide(id = "min.std")
    }
  })
  
  observe({
    req(data$data)
    parameter$cluster <- input$cluster
    parameter$standardise <- switch(input$standardise,
                                    "Z-score" = T,
                                    "Log2" = F)
    parameter$filter_std <- input$filter_std
    parameter$min.std <- input$min.std
  })
  
  return(parameter)
}

#' @export
stemplotserver <- function(input, output, session, data, parameter) {
  observe({
    output$downloadplot <- downloadHandler(
      filename = "时间序列分析.zip",
      content = function(file) {
        req(data$data)
        req(parameter)
        
        wd <- getwd()
        setwd(tempdir())
        
        name <- paste0("时间序列分析-", RandomCode(5))
        
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(filename = data$data,
                       resultPath = name),args)
        
        stemresult <- do.call(what = auto_stem, args = args)
        
        setwd("../")
        
        zip::zip(zipfile = file, files = name)
        setwd(wd)
      }
    )
  })
}

#' @export
stemui <- function(id = "stem") {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        p("趋势/时间序列分析", style = "color:black;font-size:30px;font-weight:bold;"),
        hr(), readdataInput(id = ns("xlsx")),
        hr(), downloadui(id = ns("plot"), downloadlabel = "结果保存"),
        width = 3
      ),
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "趋势/时间序列分析",
            br(),
            fixedRow(column(width = 1, offset = 11, stemsetInput(id = ns("set")))),
            h3("趋势/时间序列分析图较多，不进行展示，请直接下载结果")
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
stemserver <- function(input, output, session) {
  data <- callModule(readdataserver, "xlsx", row.names = T)
  parameter <- callModule(stemsetserver, "set", data = data)
  callModule(stemplotserver, "plot", data = data, parameter = parameter)
}


#' @export
app_stem <- function() {
  suppressMessages(library(shiny))
  suppressMessages(library(shinydashboard))
  suppressMessages(library(shinyWidgets))
  suppressMessages(library(shinyjs))
  suppressMessages(library(DT))
  suppressMessages(library(colourpicker))
  
  shinyApp(
    ui = bootstrapPage(
      stemui(id = "stem")
    ),
    server = function(input, output, session) {
      callModule(stemserver, "stem")
    }
  )
}