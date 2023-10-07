#' @export
metaanalystui <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    titlePanel("代谢分析"),
    navlistPanel(
      tabPanel(
        "第一步：信息获取",
        h3("信息获取"),
        hr(),
        textInput(inputId = ns("project_id1"), label = "项目编号:"),
        uiOutput(outputId = ns("backinfo1")),
        actionButton(ns("button1"), "目录生成"),
        actionButton(ns("button7"), "数据转格式"),
        br(),
        br(),
        textInput(inputId = ns("analysis_id1"), label = "分析编号:"),
        uiOutput(outputId = ns("backinfo2")),
        actionButton(ns("button2"), "信息获取")
      ),
      tabPanel(
        "第二步：实验报告",
        h3("实验报告"),
        hr(),
        selectInput(
          inputId = ns("project_id2"), label = "项目编号:",
          choices = list.dirs(path = projectpath(), recursive = F, full.names = F),
          selected = NULL
        ),
        selectInput(
          inputId = ns("mode1"), label = "报告模板:",
          choices = gsub(pattern = ".yaml",replacement = "",
                         x = list.files(selectfile(path = "report",file = "yamlfiles"),
                                        recursive = F, full.names = F)),
          selected = "lm"
        ),
        selectInput(
          inputId = ns("ana_type1"), label = "项目类型:",
          choices = list.dirs("/data/hstore4/database/mould/report/Oebio",
                              recursive = F, full.names = F)[!grepl("非靶向代谢",x =list.dirs("/data/hstore4/database/mould/report/Oebio",recursive = F, full.names = F))],
          selected = "全谱代谢-LCMS"
        ),
        uiOutput(outputId = ns("backinfo3")),
        actionButton(ns("button3"), "模板生成"),
        actionButton(ns("button4"), "报告生成")
      ),
      tabPanel(
        "第三步：项目分析",
        h3("项目分析"),
        hr(),
        selectInput(
          inputId = ns("project_id3"), label = "项目编号:",
          choices = list.dirs(path = projectpath(), recursive = F, full.names = F),
          selected = NULL
        ),
        selectInput(
          inputId = ns("analysis_id3"), label = "分析编号:",
          choices = NULL,
          selected = NULL
        ),
        selectInput(
          inputId = ns("mode2"), label = "报告模板:",
          choices = gsub(pattern = ".yaml",replacement = "",
                         x = list.files(selectfile(path = "report",file = "yamlfiles"),
                                        recursive = F, full.names = F)),
          selected = "lm"
        ),
        uiOutput(outputId = ns("backinfo4")),
        actionButton(ns("button5"), "质控分析"),
        actionButton(ns("button6"), "正式分析")
      ),
      widths = c(2, 10)
    )
  )
}


#' @export
metaanalystserver <- function(input, output, session) {
  observeEvent(input$button1, {
    req(input$project_id1)
    runtry <- try({
      if (!grepl(pattern = "^[0-9A-Za-z-]*$", x = input$project_id1)) {
        stop("项目编号不符合规则")
      }
      
      createproject(project_id = input$project_id1)
      
      updateSelectInput(session = session,
                        inputId = "project_id2",
                        choices = list.dirs(path = projectpath(), recursive = F, full.names = F),
                        selected = input$project_id1)
      updateSelectInput(session = session,
                        inputId = "project_id3",
                        choices = list.dirs(path = projectpath(), recursive = F, full.names = F),
                        selected = input$project_id1)
      updateSelectInput(session = session,
                        inputId = "analysis_id3",
                        choices = list.files(path = paste0(projectpath(), input$project_id1),
                                             pattern = input$project_id1, recursive = F, full.names = F))
    })
    
    if ("try-error" %in% class(runtry)) {
      output$backinfo1 <- renderUI({p(paste0(input$project_id1, "运行错误:", runtry[[1]]))})
    } else {
      output$backinfo1 <- renderUI({
        tagList(p(paste0(input$project_id1, "目录生成完成,请将质谱及搜库数据放到指定位置")),
                p(paste0("链接地址:", projectpath2(input$project_id1))))
      })
    }
  })
  
  observeEvent(input$button7, {
    req(input$project_id1)
    runtry <- try({
      mulmsconvert(project_id = input$project_id1, wait = F)
      
      updateSelectInput(session = session,
                        inputId = "project_id2",
                        choices = list.dirs(path = projectpath(), recursive = F, full.names = F),
                        selected = input$project_id1)
      updateSelectInput(session = session,
                        inputId = "project_id3",
                        choices = list.dirs(path = projectpath(), recursive = F, full.names = F),
                        selected = input$project_id1)
      updateSelectInput(session = session,
                        inputId = "analysis_id3",
                        choices = list.files(path = paste0(projectpath(), input$project_id1),
                                             pattern = input$project_id1, recursive = F, full.names = F))
    })
    
    if ("try-error" %in% class(runtry)) {
      output$backinfo1 <- renderUI({p(paste0(input$project_id1, "运行错误:", runtry[[1]]))})
    } else {
      output$backinfo1 <- renderUI({
        tagList(p(paste0(input$project_id1, "转格式中请耐心等待结果")),
                p(paste0("链接地址:", projectpath2(input$project_id1))))
      })
    }
  })
  
  
  observeEvent(input$button2, {
    req(input$project_id1)
    req(input$analysis_id1)
    runtry <- try({
      
      createanalysis(project_id = input$project_id1,
                     analysis_id = input$analysis_id1)
      
      updateSelectInput(session = session,
                        inputId = "project_id2",
                        choices = list.dirs(path = projectpath(), recursive = F, full.names = F),
                        selected = input$project_id1)
      updateSelectInput(session = session,
                        inputId = "project_id3",
                        choices = list.dirs(path = projectpath(), recursive = F, full.names = F),
                        selected = input$project_id1)
      
      # updateSelectInput(session = session,
      #                   inputId = "analysis_id2",
      #                   choices = list.files(path = paste0(projectpath(), input$project_id1),
      #                                        pattern = input$project_id1, recursive = F, full.names = F),
      #                   selected = input$analysis_id1)
      updateSelectInput(session = session,
                        inputId = "analysis_id3",
                        choices = list.files(path = paste0(projectpath(), input$project_id1),
                                             pattern = input$project_id1, recursive = F, full.names = F))
    })
    
    if ("try-error" %in% class(runtry)) {
      output$backinfo2 <- renderUI({p(paste0(input$project_id1, "运行错误:", runtry[[1]]))})
    } else {
      output$backinfo2 <- renderUI({
        tagList(p(paste0(input$project_id1, "信息获取完成")),
                p(paste0("链接地址:", projectpath2(paste0(input$project_id1, "/", input$analysis_id1)))))
      })
    }
  })
  
  
  # observeEvent(input$project_id2, {
  #   req(input$project_id2)
  #   filename <- list.files(
  #     path = projectpath(input$project_id2),
  #     pattern = input$project_id2, recursive = F, full.names = F
  #   )
  #   updateSelectInput(
  #     session = session,
  #     inputId = "analysis_id2",
  #     choices = filename
  #   )
  #   print(input$project_id2)
  #   print(filename)
  #   if (length(filename) == 0) {
  #     output$backinfo3 <- renderUI({
  #       p(paste0(input$project_id2, "下无分析编号目录"))
  #     })
  #   } else {
  #     output$backinfo3 <- renderUI({
  #       p("")
  #     })
  #   }
  # })
  
  observeEvent(input$project_id3, {
    req(input$project_id3)
    filename <- list.files(
      path = projectpath(input$project_id3),
      pattern = input$project_id3, recursive = F, full.names = F
    )
    print(input$project_id3)
    print(filename)
    updateSelectInput(
      session = session,
      inputId = "analysis_id3",
      choices = filename
    )
    
    if (length(filename) == 0) {
      output$backinfo4 <- renderUI({
        p(paste0(input$project_id3, "下无分析编号目录"))
      })
    } else {
      output$backinfo4 <- renderUI({
        p("")
      })
    }
  })
  
  observeEvent(input$mode1, {
    req(input$mode1)
    updateSelectInput(
      session = session,
      inputId = "mode2",
      choices = gsub(pattern = ".yaml",replacement = "",
                     x = list.files(selectfile(path = "report",file = "yamlfiles"),
                                    recursive = F, full.names = F)),
      selected = input$mode1
    )
  })
  
  observeEvent(input$button3, {
    req(input$project_id2)
    req(input$mode1)
    
    runtry <- try({
      base::print(input$project_id2)
      base::print(input$ana_type1)
      mkexreport_oebio(
        project_id = input$project_id2,
        type = input$ana_type1,
        excopy = T,
        exmd = F,
        mode = input$mode1)
    })
    
    if ("try-error" %in% class(runtry)) {
      output$backinfo3 <- renderUI({
        p(paste0(input$project_id2, "运行错误:", runtry[[1]]))
      })
    } else {
      output$backinfo3 <- renderUI({
        tagList(
          p(paste0(input$project_id2, "模板生成完成")),
          p(paste0("链接地址:", projectpath2(paste0(input$project_id2, "/", input$analysis_id2))))
          # a(paste0("链接地址:",projectpath2(paste0(input$project_id2,"/",input$analysis_id2))),
          #   href = projectpath2(paste0(input$project_id2,"/",input$analysis_id2)))
        )
      })
    }
  })
  
  observeEvent(input$button4, {
    req(input$project_id2)
    req(input$mode1)
    
    runtry <- try({
      base::print(input$project_id2)
      base::print(input$ana_type1)
      mkexreport_oebio(
        project_id = input$project_id2,
        type = input$ana_type1,
        excopy = F,
        exmd = T,
        mode = input$mode1)
    })
    
    if ("try-error" %in% class(runtry)) {
      output$backinfo3 <- renderUI({
        p(paste0(input$project_id2, "运行错误:", runtry[[1]]))
      })
    } else {
      output$backinfo3 <- renderUI({
        tagList(
          p(paste0(input$project_id2, "报告生成完成")),
          p(paste0("链接地址:", projectpath2(paste0(input$project_id2, "/", input$analysis_id2))))
          # a(paste0("链接地址:",projectpath2(paste0(input$project_id2,"/",input$analysis_id2))),
          #   href = projectpath2(paste0(input$project_id2,"/",input$analysis_id2)))
        )
      })
    }
  })
  
  observeEvent(input$button5, {
    req(input$project_id3)
    req(input$analysis_id3)
    req(input$mode2)
    runtry <- try({
      print(input$project_id3)
      print(input$analysis_id3)
      metanalystflow(project_id = input$project_id3,
                     analysis_id = input$analysis_id3,
                     process = "pre")
    })
    
    if ("try-error" %in% class(runtry)) {
      output$backinfo4 <- renderUI({
        p(paste0(input$project_id3, "运行错误:", runtry[[1]]))
      })
    } else {
      output$backinfo4 <- renderUI({
        tagList(
          p(paste0(input$project_id3, "初步分析中请耐心等待结果")),
          p(paste0("链接地址:", projectpath2(paste0(input$project_id3, "/", input$analysis_id3))))
          # a(paste0("链接地址:",projectpath2(paste0(input$project_id3,"/",input$analysis_id3))),
          #   href = projectpath2(paste0(input$project_id3,"/",input$analysis_id3)))
        )
      })
    }
  })
  
  observeEvent(input$button6, {
    req(input$project_id3)
    req(input$analysis_id3)
    req(input$mode2)
    runtry <- try({
      print(input$project_id3)
      print(input$analysis_id3)
      metanalystflow(project_id = input$project_id3,
                     analysis_id = input$analysis_id3,
                     process = "ana",
                     mode = input$mode2)
    })
    
    if ("try-error" %in% class(runtry)) {
      output$backinfo4 <- renderUI({
        p(paste0(input$project_id1, "运行错误:", runtry[[1]]))
      })
    } else {
      output$backinfo4 <- renderUI({
        tagList(
          p(paste0(input$project_id3, "正式分析中请耐心等待结果")),
          p(paste0("链接地址:", projectpath2(paste0(input$project_id3, "/", input$analysis_id3))))
          # a(paste0("链接地址:",projectpath2(paste0(input$project_id3,"/",input$analysis_id3))),
          #   href = projectpath2(paste0(input$project_id3,"/",input$analysis_id3)))
        )
      })
    }
  })
}

#' @export
metaanalyst_app <- function() {
  library(shiny, quietly = T)
  shinyApp(
    ui = metaanalystui(id = "metaanalyst"),
    server = function(input, output, session) {
      callModule(metaanalystserver, "metaanalyst")
    }
  )
}
