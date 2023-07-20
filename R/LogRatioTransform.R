#' @export png.iclr
png.iclr <- function(x){
  exp(x)/sum(exp(x))
}

#' @export png.clr
png.clr <- function(x){
  gx <- sum(log(x))/length(x)
  log(x)-gx
}
