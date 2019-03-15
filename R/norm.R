gm_mean = function(x){
  if(any(x < 0, na.rm = TRUE)){
    stop("Negative values not allowed")
  }
  exp(mean(log(x)))
}

#' Centered log-ratio normalization
#'
#' @param x numeric matrix. Samples are columns
#' @return A CLR normalized count_table
#' @examples norm_clr(matrix(1:100, nrow = 10, ncol = 10))
#' @export
norm_clr <- function(x){
  gm <- apply(x, 2, function(y) gm_mean(y))
  return(t(log(t(x)/gm)))
}

#' Additive log-ratio normalization
#'
#' @param x numeric matrix. Samples are columns
#' @param ref reference feature
#' @return An ALR normalized count_table
#' @examples norm_alr(matrix(1:100, nrow = 10, ncol = 10))
#' @export
norm_alr <- function(x, ref = nrow(x)){
  return(apply(x, 2, function(y) log(y/y[ref]))[-ref,])
}


