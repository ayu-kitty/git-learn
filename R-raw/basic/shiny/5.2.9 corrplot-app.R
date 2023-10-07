#' @export
corrplotsetInput <- function(id = "corrplot",
                             right = T) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    dropdownButton(
      inputId = ns("corrplotdrop"),
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
          textInput(inputId = ns("title"), label = "图展示名称:", value = "Correlation"),
          selectInput(inputId = ns("method"), label = "图形呈现方式:", choices = c("circle", "square", "ellipse", "number", "pie", "shade", "color")),
          selectInput(inputId = ns("type"), label = "图形呈现区域:", choices = c("upper", "full", "lower"), selected = "upper"),
          selectInput(inputId = ns("mode"), label = "选择标注模式:", choices = c("mode1", "mode2", "mode3", "mode4")),
          colourpicker::colourInput(ns("col3"), "选择正相关颜色:", "firebrick3", returnName = T),
          colourpicker::colourInput(ns("col2"), "选择过度颜色:", "white", returnName = T),
          colourpicker::colourInput(ns("col1"), "选择负相关颜色:", "navy", returnName = T),
          checkboxInput(ns("tl.pos"), "是否显示名称", TRUE),
          checkboxInput(ns("diag"), "是否呈现对角线", FALSE),
          numericInput(inputId = ns("tl.srt"), label = "字体角度:", value = 90)
        ),
        tabPanel(
          "运算",
          br(),
          selectInput(
            inputId = ns("order"), label = "选择排序方式:",
            choices = c("hclust", "original", "AOE", "FPC", "alphabet")
          ),
          selectInput(
            inputId = ns("adjust"), label = "P值矫正方式:",
            choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
            selected = "BH"
          ),
          checkboxInput(ns("trans"), "是否数据转置", TRUE),
        ),
        tabPanel(
          "图形",
          br(),
          p("此页面参数不设置将根据数据自动变化", style = "color:red"),
          numericInput(inputId = ns("tl.cex"), label = "字体大小:", value = NA),
          numericInput(inputId = ns("pch.cex"), label = "数字/形状大小:", value = NA),
          sliderInput(ns("width"), "图像宽度:", min = 5, max = 20, value = 10, step = 1),
          sliderInput(ns("height"), "图像高度:", min = 5, max = 20, value = 10, step = 1)
        )
      )
    )
  )
}


#' @export
corrplotsetserver <- function(input, output, session, data) {
  parameter <- reactiveValues()
  
  observe({
    req(data$data)
    
    parameter$title <- paste0("\n", input$title)
    parameter$method <- input$method
    parameter$type <- input$type
    parameter$col <- colorRampPalette(c(input$col1, input$col2, input$col3))(1000)
    parameter$tl.pos <- input$tl.pos
    parameter$diag <- input$diag
    parameter$tl.srt <- input$tl.srt
    parameter$order <- input$order
    parameter$adjust <- input$adjust
    parameter$trans <- input$trans
    parameter$tl.cex <- input$tl.cex
    parameter$pch.cex <- input$pch.cex
    parameter$number.cex <- input$pch.cex
    parameter$width <- input$width
    parameter$height <- input$height
    
    parameter$insig <- {
      if (input$mode == "mode1") {
        "label_sig"
      } else if (input$mode == "mode2") {
        "label_sig"
      } else if (input$mode == "mode3") {
        "blank"
      } else if (input$mode == "mode4") {
        "label_sig"
      }
    }
    
    parameter$sig.level <- {
      if (input$mode == "mode1") {
        c("", "", "")
      } else if (input$mode == "mode2") {
        c(0.01, 0.05)
      } else if (input$mode == "mode3") {
        0.05
      } else if (input$mode == "mode4") {
        c(0.001, 0.01, 0.05)
      }
    }
  })
  
  return(parameter)
}

#' @export
corrplotplotserver <- function(input, output, session, data, parameter) {
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
        
        returndata <- do.call(what = map_common_corrplot,args = args)
        
      },
      width = parameter$width * 72,
      height = parameter$height * 72
    )
  })
  
  observe({
    output$downloadplot <- downloadHandler(
      filename = "相关性图.zip",
      content = function(file) {
        req(data$data)
        req(parameter)
        
        wd <- getwd()
        setwd(tempdir())
        
        name <- paste0("相关性图-", RandomCode(5))
        
        setwddir(filename = name, force = T)
        
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(filename = data$data,
                       imagetype = input$imagetype,
                       mapname = "correlation"),args)
        
        returndata <- do.call(what = map_common_corrplot,args = args)
        
        setwd("../")
        zip::zip(zipfile = file, files = name)
        setwd(wd)
      }
    )
  })
  
  observe({
    output$alldownloadData <- downloadHandler(
      filename = "相关性图.zip",
      content = function(file) {
        req(data$data)
        req(parameter)
        
        wd <- getwd()
        setwd(tempdir())
        
        name <- paste0("相关性图-", RandomCode(5))
        
        setwddir(filename = name, force = T)
        
        args <- reactiveValuesToList(parameter)
        args <- args[!is.na(args)]
        args <- c(list(filename = data$path,
                       imagetype = input$imagetype),args)
        
        returndata <- do.call(what = map_common_corrplot,args = args)
        
        setwd("../")
        
        zip::zip(zipfile = file, files = name)
        
        setwd(wd)
      }
    )
  })
}

#' @export
corrplotui <- function(id = "corrplot") {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        p("相关性图", style = "color:black;font-size:30px;font-weight:bold;"),
        hr(), readdataInput(id = ns("xlsx")),
        hr(), allplotdownloadui(id = ns("plot")),
        width = 3
      ),
      mainPanel(
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "相关性图展示",
            br(),
            fixedRow(column(width = 1, offset = 11, corrplotsetInput(id = ns("set")))),
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
corrplotserver <- function(input, output, session) {
  data <- callModule(readdataserver, "xlsx",row.names = T)
  parameter <- callModule(corrplotsetserver, "set", data = data)
  callModule(corrplotplotserver, "plot", data = data, parameter = parameter)
}


#' @export
app_corrplot <- function() {
  suppressMessages(library(shiny))
  suppressMessages(library(shinydashboard))
  suppressMessages(library(shinyWidgets))
  suppressMessages(library(shinyjs))
  suppressMessages(library(DT))
  suppressMessages(library(colourpicker))
  
  shinyApp(
    ui = bootstrapPage(
      corrplotui(id = "corrplot")
    ),
    server = function(input, output, session) {
      callModule(corrplotserver, "corrplot")
    }
  )
}
