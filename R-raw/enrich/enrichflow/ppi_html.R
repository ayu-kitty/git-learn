#!/opt/conda/bin/Rscript

#' PPI 动态网络图
#'
#' @param savepath 储存路径
#' @param conmpare 比较组名称
#' @param fontfamily 字体
#' @param ... 
#' @export
ppi_html<- function(ppi_node=ppi_node,ppi_fc=ppi_fc,savepath = "./ppi/", compare = conmpare , fontfamily = "sans", ...){
  pacman::p_load(networkD3,htmltools,dplyr,stringr)
  diff2<- diff1[,c("Accession","description")]
  ppi_fc<- ppi_fc[which(ppi_fc$Degree!=0),]
  ppi_fc$num<- seq(from=0, to=nrow(ppi_fc)-1)
  ppi_fc$group<- ifelse("Foldchange"%in%names(ppi_fc),ifelse(ppi_fc$Foldchange>=1,"Up","Down"),"All")
  ppi_node$source <- sapply(ppi_node$node1,function(x){
    ppi_fc$num[which(x==ppi_fc$Accession)]
  })
  ppi_node$target<- sapply(ppi_node$node2,function(x){
    ppi_fc$num[which(x==ppi_fc$Accession)]
  })
  ppi_node$value<- ppi_node$combined_score*0.005
  Links<- ppi_node[,c("source","target","value")]
  Nodes <- as.data.frame(unclass(ppi_fc), stringsAsFactors = TRUE)
  
  
  ColourScale <- 'd3.scaleOrdinal()
            .domain(["Up", "Down"])
           .range(["#E6A0C4", "#7294D4"]);'
  script <- 'd3.select(this).select("circle").transition().duration(750).attr("r", 30);
confirm("Acceccion  :  " +  d.name + "\\n" +
                 "Gene Name  :  " +  (d.GN) + "\\n" +
                 "Foldchange  :  " +  (d.FC) + "\\n" +
                 "Type  :  " +  (d.group) + "\\n" );'
  
  fp<- forceNetwork(Links = Links,#边数据集
                    Nodes = Nodes,# 节点数据集
                    Source = 'source',#边数据集中起点对应的列
                    Target = 'target',# 边数据集中终点对应的列
                    Value = 'value',# 边数据集中边的宽度对应的列
                    NodeID = 'Accession',# 节点数据集中节点名称对应的列
                    Group = "group",# 节点数据集中节点分组对应的列
                    width = 1200,# 图宽度
                    height = 1000,# 图高度
                    zoom = T,# 图是否允许缩放
                    bounded=T,# 图是否有边界
                    legend=T,# 图是否显示图例
                    linkDistance=50,
                    opacity = 1,
                    opacityNoHover = 0.5,# 鼠标没有停留时其他节点名称的透明度
                    charge=-240, #节点斥力大小(负值越大斥力越大)
                    colourScale= JS(ColourScale),# 节点颜色
                    Nodesize = "Degree",# 节点比例大小
                    fontFamily = fontfamily,# 节点名称的字体
                    fontSize = 16,# 节点名称的字号
                    arrows = F,# 边是否显示箭头
                    linkColour = "grey",# 边颜色，Cols可以是一个预先设置的列表
                    clickAction = script)
  fp$x$nodes$FC<- Nodes$Foldchange
  fp$x$nodes$GN<- Nodes$Gene.Name
  fp$x$nodes$group<- Nodes$group
  
  fp <- htmlwidgets::onRender(
    fp,'function(el,x){
debugger;
  var optArray = [];
  for (var i = 0; i < x.nodes.name.length - 1; i++) {
    optArray.push(x.nodes.name[i]);
  }
  optArray = optArray.sort();
  $(function () {
    $("#search").autocomplete({
      source: optArray
    });
  });
  d3.select(".ui-widget button").node().onclick=searchNode;
  function searchNode() {
    debugger;
    //find the node
    var selectedVal = document.getElementById("search").value;
    var svg = d3.select(el).select("svg");
    var node = d3.select(el).selectAll(".node");
    if (selectedVal == "none") {
      node.style("stroke", "white").style("stroke-width", "2");
    } else {
      var selected = node.filter(function (d, i) {
        return d.name != selectedVal;
      });
      selected.style("opacity", "0");
      var link = svg.selectAll(".link")
      link.style("opacity", "0");
      d3.selectAll(".node, .link").transition()
        .duration(5000)
        .style("opacity", 2);
    }
  }
}'
  )
  network<-
    browsable(
      attachDependencies(
        tagList(
          tags$style('
        body{background-color: #ffffff !important}
        .nodetext{fill: #000000}
        .legend text{fill: #000000}'),
          HTML('
  <div class="ui-widget">
      <input id="search">
      <button type="button">Search</button>
  </div>'),
          tags$h1(
            inputfile),fp),
        list(
          rmarkdown::html_dependency_jquery(),
          rmarkdown::html_dependency_jqueryui())))
  save_html(network, file = paste0(path,"/network.html"),libdir = "lib")
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn=-1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()
  # 基本参数
  parser$add_argument("-fa","--fontfamily",default = "sans", help = "字体,默认为Arial")
  parser$add_argument("-in","--inputfile",default = "group", help = "比较组名称")
  args <- parser$parse_args()
  ppi_network_anno <- do.call(ppi_network_anno,args = args)
}
