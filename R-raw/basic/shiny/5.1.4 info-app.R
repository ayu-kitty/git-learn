#!/opt/conda/bin/Rscript

#' @export
infoui <- function(id = "info") {
  ns <- NS(id)
  tagList(
    mainPanel(
      img(src = "https://www.lumingbio.com/upload/202011/1604899503.jpg", width = "100%"),
      h2("公司简介"),
      p("上海鹿明生物科技有限公司多年来，一直专注于生命科学和生命技术领域，
      是国内早期开展以蛋白组和代谢组为基础的多层组学整合实验与分析的团队。
      经过多年的发展沉淀，公司建立起了4D-DIA/LFQ/PRM、iTRAQ/TMT、DIA、PRM、
      修饰蛋白组等蛋白组学技术平台和空间代谢组学、全谱代谢组、靶向代谢组、
      拟靶向代谢组、脂质组、精准靶向等代谢组学技术平台以及相应的数据整合分析平台，
      并建立了科学完整的服务流程和精细规范的操作标准。同时鹿明生物的蛋白组学，
      代谢组学技术广泛应用于疾病标志物发掘、分型诊断、精准用药、药代药动、药物表征等多个领域。"),
      p("2021年8月上海鹿明生物科技有限公司（欧易生物子公司）搬迁新办公场地
        与上海欧易生物医学科技有限公司（总公司）在一幢大楼。)"),
      p("公司地址位于：上海市闵行区联航路1188号25幢5号楼"),
      br(),
      h3("联系方式:"),
      h3(a("lujw@lumingbio.com")),
      img(src = "https://www.lumingbio.com/uploadfile/202102/643bb8535758c00.png", height = "200x"),
      h4("提醒:仅内部使用，请勿商用。", style = "color:red"),
      h4("受限于服务器,如卡顿请见谅。", style = "color:red")
    )
  )
}


#' @export
infoserver <- function(input, output, session) {
}

#' @export
app_info <- function() {
  suppressMessages(library(shiny))
  suppressMessages(library(shinydashboard))
  # suppressMessages(library(Biobase))
  
  shinyApp(
    ui = bootstrapPage(infoui(id = "info")),
    server = function(input, output, session) {
      callModule(infoserver, "info")
    }
  )
}
