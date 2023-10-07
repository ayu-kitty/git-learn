#!/opt/conda/bin/Rscript

#' @export
getsuitrow_col <- function(n = 21,
                           num = 10){
  row <- ((n-1) %/% num +1)
  col <- ceiling(n / row)
  
  return(c(row,col))
}

#' @export
getsuitlab_range <- function(maxn,minn){
  
  if ((maxn - minn) == 0) {
    x <- maxn
    lab <- c(x)
    ran <- c(3, 3)
  } else if ((maxn - minn) == 1) {
    lab <- c(minn, minn + 1)
    ran <- c(3, 4)
  } else if ((maxn - minn) == 2) {
    lab <- c(minn, minn + 1, minn + 2)
    ran <- c(3, 5)
  } else if ((maxn - minn) %% 5 == 0 || (maxn - minn + 1) %% 5 == 0) {
    x <- round((maxn - minn) / 5)
    lab <- c(minn, minn + x, minn + 2 * x, minn + 3 * x, minn + 4 * x, minn + 5 * x)
    ran <- c(1, 5)
  } else if ((maxn - minn) %% 4 == 0 || (maxn - minn + 1) %% 4 == 0) {
    x <- round((maxn - minn) / 4)
    lab <- c(minn, minn + x, minn + 2 * x, minn + 3 * x, minn + 4 * x)
    ran <- c(1, 5)
  } else if ((maxn - minn) %% 3 == 0 || (maxn - minn + 1) %% 3 == 0) {
    x <- round((maxn - minn) / 3)
    lab <- c(minn, minn + x, minn + 2 * x, minn + 3 * x)
    ran <- c(2, 5)
  } else {
    x <- floor((maxn - minn) / 4)
    lab <- c(minn, minn + x, minn + 2 * x, minn + 3 * x, minn + 4 * x)
    ran <- c(1, 5)
  }
  
  return(list(lab = lab,
              range = ran))
}
