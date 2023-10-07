#!/opt/conda/bin/Rscript

#' @export
compoundidui <- function(id = "compoundid") {
  ns <- NS(id)
  
  tagList(
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "获取代谢物信息",
            br(),
            textAreaInput(inputId = ns("Metabolites1"), label = "代谢物列表:", height = "279px", resize = "none"),
            # selectInput(inputId = ns("listfrom1"), label = "代谢物方式:", choices = c("代谢物名", "ID")),
            pickerInput(inputId = ns("metafrom1"), label = "信息来源:",
                        selected =  c("inchikey","metlin","hmdb","lipidmaps","smiles","cas","kegg","英文名","中文名"),
                        choices = c("inchikey","metlin","hmdb","lipidmaps","smiles","cas","kegg","英文名","中文名","来源"),
                        multiple = T,options = pickerOptions(actionsBox = TRUE, 
                                                             size = 10,
                                                             selectedTextFormat = "count > 3")),
            materialSwitch(inputId = ns("vague1"), label = "是否模糊匹配", value = F, status = "primary"),
            materialSwitch(inputId = ns("distinct1"), label = "是否显示全部", value = F, status = "primary"),
            actionButton(inputId = ns("updata1"), label = "提交")
          ),
          tabPanel(
            "获取mz信息",
            br(),
            p("最多提取50个代谢物"),
            textAreaInput(inputId = ns("Metabolites2"), label = "代谢物列表:", height = "279px", resize = "none"),
            # selectInput(inputId = ns("listfrom2"), label = "代谢物方式:", choices = c("代谢物名", "ID")),
            pickerInput(inputId = ns("metafrom2"), label = "信息来源:",
                        selected =  c("inchikey","metlin","hmdb","lipidmaps","smiles","cas","kegg","英文名","中文名"),
                        choices = c("inchikey","metlin","hmdb","lipidmaps","smiles","cas","kegg","英文名","中文名","来源"),
                        multiple = T,options = pickerOptions(actionsBox = TRUE, 
                                                             size = 10,
                                                             selectedTextFormat = "count > 3")),
            actionButton(inputId = ns("updata2"), label = "提交")
          ),
          tabPanel(
            "通过mz获取代谢物",
            br(),
            textAreaInput(inputId = ns("mz1"), label = "mz列表:", height = "279px", resize = "none"),
            selectInput(inputId = ns("listfrom4"), label = "正负离子形式:", choices = c("neg", "pos")),
            numericInput(inputId = ns("ppm1"), label = "ppm值:", value = 5),
            actionButton(inputId = ns("updata4"), label = "提交")
          )
        ),
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
compoundidserver <- function(input, output, session) {
  
  output$compoundid <- renderDT({})
  
  # hmdb转换
  observeEvent(input$updata1, {
    req(input$Metabolites1)
    
    output$compoundid <- renderDT(
      {
        data <- strsplit(isolate(input$Metabolites1), split = "\n")[[1]]
        data <- data.frame("ID" = data)
        data <- data[!is.na(data[, 1]), , drop = F]
        data <- data[data[, 1] != "\"", , drop = F]
        data[,1] <- gsub(pattern = "\"",replacement = "",x = data[,1])
        data[,1] <- gsub(pattern = "^\\s+|\\s+$",replacement = "",x = data[,1])
        data <- data[data[, 1] != "", , drop = F]
        data <- data[!duplicated(data[, 1]), , drop = F]
        
        data <- getmetainfo(data = data,
                            databasefrom =isolate(input$metafrom1),
                            idlist = "ID",
                            vague = isolate(input$vague1),
                            distinct = isolate(input$distinct1))
        
        data
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
  
  
  observeEvent(input$updata2, {
    req(input$Metabolites2)
    
    output$compoundid <- renderDT(
      {
        data <- strsplit(isolate(input$Metabolites2), split = "\n")[[1]]
        data <- data.frame("ID" = data)
        data <- data[!is.na(data[, 1]), , drop = F]
        data <- data[data[, 1] != "\"", , drop = F]
        data[,1] <- gsub(pattern = "\"",replacement = "",x = data[,1])
        data[,1] <- gsub(pattern = "^\\s+|\\s+$",replacement = "",x = data[,1])
        data <- data[data[, 1] != "", , drop = F]
        data <- data[!duplicated(data[, 1]), , drop = F]
        if(dim(data)[1] > 50){
          data <- data[1:50,,drop =F]
        }
        
        data <- getmetainfo(data = data,
                            databasefrom =isolate(input$metafrom2),
                            idlist = "ID")
        data <- getmetamz(data)
        
        data
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
  
  
  # hmdb转换
  observeEvent(input$updata4, {
    req(input$mz1)
    
    output$compoundid <- renderDT(
      {
        mz <- strsplit(isolate(input$mz1), split = "\n")[[1]]
        mz <- as.numeric(mz)
        # print(mz)
        data <- AddQualitative(mzdata = mz,
                               mode = isolate(input$listfrom4),
                               ppm = isolate(input$ppm1),
                               addmzinfo = T)
        # print(data)
        data
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
app_compundid <- function() {
  suppressMessages(library(shiny, quietly = T))
  suppressMessages(library(DT, quietly = T))
  suppressMessages(library(shinydashboard, quietly = T))
  suppressMessages(library(shinyWidgets))
  
  shinyApp(
    ui = bootstrapPage(compoundidui(id = "compoundid")),
    server = function(input, output, session) {
      callModule(compoundidserver, "compoundid")
    }
  )
}
