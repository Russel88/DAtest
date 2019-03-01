#' Print results from \code{testDA}
#'
#' @param x The output from the \code{testDA} function
#' @param ... Additional arguments for \code{print}
#' @return Prints results from testDA
#' @export

print.DA <- function(x, ...){
  
  xx <- x$table
  print(xx, row.names = FALSE, ...)
  
}
