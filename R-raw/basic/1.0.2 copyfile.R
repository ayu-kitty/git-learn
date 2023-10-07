#! /opt/conda/bin/Rscript

#' @export
file.copy <- function(from,to,...){
  createdir(filename = to)
  base::file.copy(from = from,to = to,...)
  Sys.chmod(paths = to,mode = "0777",use_umask = F)
}
