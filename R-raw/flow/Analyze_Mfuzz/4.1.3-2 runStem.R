#!/opt/conda/bin/Rscript

#' runStem
#'
#' 时序分析，只需要在路径下准备内部分析单,数据矩阵和差异代谢物
#'
#' @param dif_file 差异表达矩阵的表格路径，默认为差异表达矩阵.xlsx
#' @param mtrx_file 数据矩阵的路径，默认为数据矩阵.xlsx
#' @param grp_order 时序组的顺序（必须）
#' @param registr 内部分析单/项目登记单的路径，默认为内部分析单.xlsx
#' @param cluster 聚类数量,默认为2的（组数-1）次幂
#' @param sep 表单名（比较组）的分割符
#' @param protein 是否为蛋白项目
#' @param filter 是否对表单进行筛选
#' @param ...
#'
#' @export
runStem <- function ( dif_file = '差异表达矩阵.xlsx',
                      mtrx_file = '数据矩阵.xlsx',
                      grp_order = c(),
                      registr = "内部分析单.xlsx",
                      cluster = 2 ** (length(grp_order) - 1),
                      sep = ifelse(
                        grepl("内部分析单", basename(registr)),
                        "-vs-", "/"
                      ),
                      protein = F,
                      filter = T,
                      ...
          ){
    suppressMessages(library(tidyverse))
    suppressMessages(library(dplyr))
    suppressMessages(library(stringr))
  
  # if(!any(grp_order)){ stop("未获取到时间组顺序，请使用参数grp_order") }
  registr_type <-  basename(registr)
  if(grepl("内部分析单", registr_type)){
    # 获得分组信息
    group_df <- lmbio::readxlsx(registr, sheet = '样本基本信息') %>% 
                select("样本分析名称", "分组")
  } else if (grepl("项目登记单", registr_type)){
    group_df <- lmbio::readxlsx(registr, sheet = '分组信息') %>% 
                select("作图编号", "样品分组")
  } else if (grepl("Sample_Group", registr_type)){
    group_df <- lmbio::readxlsx(registr, sheet = '样品')
  } else {
    stop("错误的分析单类型")
  }
  colnames(group_df) <- c('samples','group')
  
  name_label <- ifelse(protein, "Accession", "Metabolites")
  usecols <- c(name_label, group_df$samples)  
  usecols <- usecols[!grepl("^QC", usecols, ignore.case=T)]
  
  sheets <- getsheetname(dif_file)

  union_data <- data.frame()
  
  # 遍历所有sheet, 并取并集
  for (sheet in sheets) {
    if(filter & (
       !str_split(sheet, sep)[[1]][1] %in% grp_order | 
       !str_split(sheet, sep)[[1]][2] %in% grp_order) 
       ) next
    sheet_data <- lmbio::readxlsx(dif_file, sheet = sheet)
    sheet_data <- sheet_data[colnames(sheet_data) %in% usecols]
    # 将当前sheet的数据添加到union_data数据框
    union_data <- bind_rows(union_data, sheet_data)
    union_data <- union_data[!duplicated(union_data[[name_label]]), ]
  }
  
  whole_data_matrix <- lmbio::readxlsx(mtrx_file, sheet = '数据矩阵') 
  whole_data_matrix <- whole_data_matrix[, colnames(whole_data_matrix) %in% usecols]
  all_data <- whole_data_matrix[whole_data_matrix[[name_label]] %in% union_data[[name_label]], ]
  
  # 转置
  all_data <- all_data %>%
              t() %>%
              as.data.frame()
  
  # 将第一行数据转化为列名
  colnames(all_data) <-  all_data[1, ]
  all_data <- all_data[-1, ]
  all_data$samples <- rownames(all_data)
  
  #合并数据框并转化为数值类型
  merged_df <- left_join(all_data, group_df, by = 'samples')
  columns_to_convert <- names(merged_df)[1 : (length(names(merged_df)) - 2)]
  merged_df[columns_to_convert] <- as.data.frame(lapply(merged_df[columns_to_convert],as.numeric))
  
  #将samples列转化为行名
  rownames(merged_df) <-  merged_df$samples
  merged_df <- merged_df %>% select(-"samples")
  
  # 计算均值（将空值替换为最小值）
  averages <- merged_df %>%
    group_by(group) %>%
    mutate_at(vars(-group), 
              ~ ifelse(is.na(.),  
                        merged_df %>% select(-group) %>% min(na.rm = TRUE),
                        .)) %>%
    summarise_all("mean") %>%
    t() %>%
    as.data.frame()
  
  colnames(averages) <- averages[1, ]
  averages <- averages[-1, ]
  averages[] <- lapply(averages, as.numeric)
  averages <- averages[, grp_order] 
  
  filename = '趋势分析矩阵.xlsx'
  openxlsx::write.xlsx(averages, file = filename, rowNames = T)
  
  #运行时序分析
  lmbio::auto_stem(filename = filename, cluster = cluster, ...)
}

if(this.path::is.main() & !is.null(whereami::thisfile_r())){
  options(warn = -1)
  suppressMessages(library("lmbio"))
  parser <- ArgumentParser()

  parser$add_argument("-d", "--dif_file", default = '差异表达矩阵.xlsx', help = "差异表达矩阵的表格路径，默认为差异表达矩阵.xlsx")
  parser$add_argument("-m", "--mtrx_file", default = "数据矩阵.xlsx", help = "数据矩阵的路径，默认为数据矩阵.xlsx")
  parser$add_argument("-g", "--grp_order", default = c(), help = "时序组的顺序（必须）", nargs="+")
  parser$add_argument("-r", "--registr", default = "内部分析单.xlsx", help = "内部分析单/项目登记单的路径，默认为内部分析单.xlsx")
  parser$add_argument("-c", "--cluster", default = 0, type = "integer", help = "聚类数量，默认为2的（组数-1）次幂")
  parser$add_argument("-s", "--sep", default = "-vs-", help = "表单名（比较组）的分割符，默认为-vs-")
  parser$add_argument("-p", "--protein", default = F, type = "logical", help = "是否为蛋白项目，默认为F")
  parser$add_argument("-f", "--filter", default = T, type = "logical", help = "是否对表单进行筛选，默认为T")
  
  args <- parser$parse_args()
  if(args$cluster == 0){
    args$cluster <- 2 ** (length(args$grp_order) - 1)
  }
  do.call(runStem, args = args)
}

