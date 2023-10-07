
#' @export
cleanenv <- function(pattern = "^arg"){
  pacman::p_load((.packages()), character.only = TRUE)
  pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
  value <- ls(pos = 1,all.names = T)
  value <- value[!grepl(pattern = pattern,x = value)]
  rm(list = c(value),pos = 1)
  gc(reset = TRUE)
  suppressMessages(library("lmbio"))
}

