#' @export
randomForestsetInput <- function(id = "randomForest",
                                 right = T) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    dropdownButton(
      inputId = ns("randomForestdrop"),
      label = "设置",
      icon = icon("gear"),
      status = "primary",
      circle = FALSE,
      right = right,
      sliderInput(inputId = ns("trainRatio"), label = "训练/测试比例:", min = 0.1, max = 1, step = 0.1, value = 1),
      sliderInput(inputId = ns("importanceTop"), label = "呈现TOP数量:", min = 10, max = 50, step = 10, value = 30)
    )
  )
}


#' @export
randomForestsetserver <- function(input, output, session, data) {
  parameter <- reactiveValues()
  
  observe({
    req(data$data)
    parameter$trainRatio <- input$trainRatio
    parameter$importanceTop <- input$importanceTop
  })
  
  return(parameter)
}

#' @export
randomForestplotserver <- function(input, output, session, data, parameter) {
  observe({
    output$downloadplot <- downloadHandler(
      filename = "随机森林分析.zip",
      content = function(file) {
        req(data$data)
        req(parameter)
        
        wd <- getwd()
        setwd(tempdir())
        
        name <- paste0("随机森林分析-", RandomCode(5))
        
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(filename = data$data,
                       resultPath = name),args)
        
        randomforestresults <- do.call(what = randomForestAnalyst, args = args)
        
        setwd("../")
        
        zip::zip(zipfile = file, files = name)
        setwd(wd)
      }
    )
  })
}

#' @export
randomForestui <- function(id = "randomForest") {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        p("随机森林分析", style = "color:black;font-size:30px;font-weight:bold;"),
        hr(), readdataInput(id = ns("xlsx")),
        hr(), downloadui(id = ns("plot"), downloadlabel = "结果保存"),
        width = 3
      ),
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "随机森林分析",
            br(),
            fixedRow(column(width = 1, offset = 11, randomForestsetInput(id = ns("set")))),
            h3("随机森林分析图较多，不进行展示，请直接下载结果")
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
randomForestserver <- function(input, output, session) {
  data <- callModule(readdataserver, "xlsx", row.names = T)
  parameter <- callModule(randomForestsetserver, "set", data = data)
  callModule(randomForestplotserver, "plot", data = data, parameter = parameter)
}


#' @export
app_randomForest <- function() {
  suppressMessages(library(shiny))
  suppressMessages(library(shinydashboard))
  suppressMessages(library(shinyWidgets))
  suppressMessages(library(shinyjs))
  suppressMessages(library(DT))
  
  shinyApp(
    ui = bootstrapPage(
      randomForestui(id = "randomForest")
    ),
    server = function(input, output, session) {
      callModule(randomForestserver, "randomForest")
    }
  )
}