#!/opt/conda/bin/Rscript

#' @export
tcmui <- function(id = "tcm") {
  ns <- NS(id)
  
  tagList(
    sidebarLayout(
      sidebarPanel(
        p("中药数据库查询", style = "color:black;font-size:30px;font-weight:bold;"),
        hr(),
        textAreaInput(inputId = ns("Metabolites"), label = "代谢物列表:", height = "279px", resize = "none"),
        selectInput(inputId = ns("listfrom"), label = "代谢物方式:", choices = c("中文名", "英文名")),
        selectInput(inputId = ns("match"), label = "匹配方式:", choices = c("精确匹配","模糊匹配")),
        actionButton(inputId = ns("updata"), label = "提交"),
        width = 3
      ),
      mainPanel(
        shinycssloaders::withSpinner(
          DTOutput(outputId = ns("compoundid")),
          type = 4
        )
      )
    )
  )
}

#' @export
tcmserver <- function(input, output, session) {
  
  output$compoundid <- renderDT({})
  
  # hmdb转换
  observeEvent(input$updata, {
    req(input$Metabolites)
    
    output$compoundid <- renderDT(
      {
        if(isolate(input$listfrom) == "中文名"){
          wherename <- "cn_name"
        }else{
          wherename <- "compound_name"
        }
        
        data <- strsplit(isolate(input$Metabolites), split = "\n")[[1]]
        data <- data.frame("Metabolites" = data)
        data <- data[!is.na(data[, 1]), , drop = F]
        data <- data[data[, 1] != "\"", , drop = F]
        data[,1] <- gsub(pattern = "\"",replacement = "",x = data[,1])
        data[,1] <- gsub(pattern = "^\\s+|\\s+$",replacement = "",x = data[,1])
        data <- data[data[, 1] != "", , drop = F]
        data <- data[!duplicated(data[, 1]), , drop = F]
        
        if(isolate(input$match) == "精确匹配"){
          info <- getmysqldata(dbname = "cosa",
                               table = "herb_annotation",
                               wherename = wherename,
                               wheredata = data[,1])
        }else{
          info <- getmysqldata_vague(dbname = "cosa",
                                     table = "herb_annotation",
                                     wherename = wherename,
                                     wheredata = data[,1])
        }

        info
      },
      extensions = "Buttons",
      rownames = F,
      options = list(
        scrollX = TRUE,
        scrollCollapse = TRUE,
        scrollY = "500px",
        paging = F,
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel")
      )
    )
  })
}

#' @export
app_tcm <- function() {
  suppressMessages(library(shiny, quietly = T))
  suppressMessages(library(DT, quietly = T))
  suppressMessages(library(shinydashboard, quietly = T))
  
  shinyApp(
    ui = bootstrapPage(tcmui(id = "tcm")),
    server = function(input, output, session) {
      callModule(tcmserver, "tcm")
    }
  )
}
