
#' @export
addblack_color <- function(color,num = 50,...){
  red <- col2rgb(color)[1]
  if(red > num){
    if((red - num) > 255){
      red <- 255
    }else{
      red <- red - num
    }
  }else{red <- 0}
  green <- col2rgb(color)[2]
  if(green > num){
    if((green - num) > 255){
      green <- 255
    }else{
      green <- green - num
    }
  }else{green <- 0}
  blue <- col2rgb(color)[3]
  if(blue > num){
    if((blue - num) > 255){
      blue <- 255
    }else{
      blue <- blue - num
    }
  }else{blue <- 0}
  
  color <- rgb(red = red,green = green,blue = blue,maxColorValue = 255,...)
  
  return(color)
}




