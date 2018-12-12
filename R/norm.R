#' @export
gm_mean = function(x){
  if(any(x < 0, na.rm = TRUE)){
    stop("Negative values not allowed")
  }
  exp(mean(log(x)))
}

#' @export
norm_clr <- function(x){
  gm <- apply(x, 2, function(y) gm_mean(y))
  return(t(log(t(x)/gm)))
}

#' @export
norm_alr <- function(x, ref = nrow(x)){
  return(apply(x, 2, function(y) log(y/y[ref]))[-ref,])
}


