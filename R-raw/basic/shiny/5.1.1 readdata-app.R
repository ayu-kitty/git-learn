#!/opt/conda/bin/Rscript

#' 文件读取
#' 
#' @param id ui编号 
#'
#' @export
readdataInput <- function(id = "xlsx") {
  ns <- NS(id)
  tagList(
    fileInput(ns("file"), "上传文件", multiple = F, accept = c(".xlsx", ".txt", ".csv")),
    selectInput(ns("sheet"), "选择sheet文件:", choices = NULL)
  )
}

#' 文件读取后，选择物质
#' 
#' @param id ui编号 
#'
#' @export
readdataInput2 <- function(id = "xlsx", maxnumber = 5) {
  ns <- NS(id)
  tagList(
    fileInput(ns("file"), "上传文件", multiple = F, accept = c(".xlsx", ".txt", ".csv")),
    selectInput(ns("sheet"), "选择sheet文件:", choices = NULL),
    if(is.null(maxnumber)){
      selectInput(ns("feature"), "选择物质:", choices = NULL, multiple = F)
    }else{
      selectInput(ns("feature"), paste0("选择物质(最多", maxnumber, "个):"), choices = NULL, multiple = T)
    }
  )
}

#' 文件读取后，选择列名
#' 
#' @param id ui编号 
#'
#' @export
readdataInput3 <- function(id = "xlsx") {
  ns <- NS(id)
  tagList(
    fileInput(ns("file"), "上传文件", multiple = F, accept = c(".xlsx")),
    selectInput(ns("sheet"), "选择sheet文件:", choices = NULL, multiple = T),
    selectInput(ns("colname"), "选择共同列名:", choices = NULL)
  )
}

#' 文件读取后，选择列名
#' 
#' @param id ui编号 
#'
#' @export
readdataOutput <- function(id = "xlsx") {
  ns <- NS(id)
  tagList(
    shinycssloaders::withSpinner(
      DTOutput(ns("data")),
      type = 4
    )
  )
}


#' @export
readdataserver <- function(input, output, session, row.names = F) {
  
  data <- reactiveValues(data = NULL,
                         path = NULL,
                         sheet = NULL)
  
  observe({
    req(input$file)
    data$path <- input$file$datapath
    updateSelectInput(session,
                      inputId = "sheet",
                      label = "选择sheet文件:",
                      choices = getsheetname(input$file$datapath))
  })
  
  observe({
    req(input$file, input$sheet)
    
    data$data <- readdata(filename = input$file$datapath,
                          sheet = input$sheet,
                          row.names = row.names)
    
    data$sheet <- input$sheet
  })
  
  output$data <- renderDT({
    req(data$data)
    format(data$data, scientific = TRUE, digit = 2)
  },
  options = list(scrollX = TRUE),
  rownames = row.names
  )
  
  return(data)
}

#' @export
readdataserver2 <- function(input, output, session, row.names = T, maxnumber = 5) {
  
  data <- reactiveValues(data = NULL,
                         path = NULL,
                         sheet = NULL)
  
  observe({
    req(input$file)
    data$path <- input$file$datapath
    updateSelectInput(session,
                      inputId = "sheet",
                      label = "选择sheet文件:",
                      choices = getsheetname(input$file$datapath))
  })
  
  observe({
    req(input$file, input$sheet)
    
    data$rawdata <- readdata(filename = input$file$datapath,
                             sheet = input$sheet,
                             row.names = row.names)
    
    updateSelectInput(session,
                      inputId = "feature",
                      choices = row.names(data$rawdata)[-1])
  })
  
  observe({
    req(input$file, input$sheet, input$feature)
    
    range <- input$feature
    if(!is.null(maxnumber)){
      if (length(range) > maxnumber) {
        range <- range[1:maxnumber]
      }
    }
    
    data$data <- data$rawdata[c("Group", range),]
    
  })
  
  output$data <- renderDT({
    req(data$data)
    format(data$data, scientific = TRUE, digit = 2)
  },
  options = list(scrollX = TRUE),
  rownames = row.names
  )
  
  return(data)
}


#' @export
readdataserver3 <- function(input, output, session, row.names = F) {
  data <- reactiveValues(data = NULL,
                         path = NULL,
                         sheet = NULL,
                         colname = NULL)
  
  observe({
    req(input$file)
    data$path <- input$file$datapath
    updateSelectInput(session,
                      inputId = "sheet",
                      label = "选择sheet文件:",
                      choices = getsheetname(input$file$datapath))
  })
  
  observe({
    req(input$file, input$sheet)
    
    data$data <- readdata(filename = input$file$datapath,
                          sheet = input$sheet,
                          row.names = row.names)
    
    if (is.data.frame(data$data)) {
      choicesname <- NULL
    } else {
      choicesname <- colnames(data$data[[1]])
      if (length(data$data) > 1) {
        for (i in 2:length(data$data)) {
          choicesname <- intersect(choicesname, colnames(data$data[[i]]))
        }
      }
    }
    
    if (length(choicesname) != 0) {
      updateSelectInput(session,
                        inputId = "colname",
                        label = "选择共同列名:",
                        choices = choicesname)
    } else {
      updateSelectInput(session,
                        inputId = "colname",
                        label = "选择共同列名:",
                        choices = NULL)
    }
  })
  
  observe({
    req(input$file, input$sheet, input$colname)
    data$colname <- input$colname
  })
  
  return(data)
}

#' @export
app_readdata <- function() {
  suppressMessages(library(shiny))
  suppressMessages(library(shinydashboard))
  suppressMessages(library(DT))
  
  shinyApp(
    ui = bootstrapPage(readdataInput(id = "xlsx"),
                       readdataOutput(id = "xlsx")),
    server = function(input, output, session) {
      data <- callModule(readdataserver, "xlsx")
    }
  )
}

#' @export
app_readdata2 <- function() {
  suppressMessages(library(shiny))
  suppressMessages(library(shinydashboard))
  suppressMessages(library(DT))
  
  shinyApp(
    ui = bootstrapPage(readdataInput2(id = "xlsx"),
                       readdataOutput(id = "xlsx")),
    server = function(input, output, session) {
      data <- callModule(readdataserver2, "xlsx")
    }
  )
}

#' @export
app_readdata3 <- function() {
  suppressMessages(library(shiny))
  suppressMessages(library(shinydashboard))
  suppressMessages(library(DT))
  
  shinyApp(
    ui = bootstrapPage(readdataInput3(id = "xlsx"),
                       verbatimTextOutput("output")),
    server = function(input, output, session) {
      data <- callModule(readdataserver3, "xlsx")
      output$output <- renderPrint({reactiveValuesToList(data)})
    }
  )
}
