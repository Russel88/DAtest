#' Summary of results from \code{testDA}
#'
#' @param object The output from the \code{testDA} function
#' @param sort Sort methods by \code{c("AUC","FPR","Spike.detect.rate")}
#' @param ... Additional arguments for \code{print}
#' @import stats
#' @import methods
#' @import utils
#' @export

summary.DA <- function(object, sort = "AUC", ...){
  
  # Find medians
  output.summary.auc <- aggregate(AUC ~ Method, data = object$table, FUN = function(x) round(median(x),3))
  output.summary.fpr <- aggregate(FPR ~ Method, data = object$table, FUN = function(x) round(median(x),3))
  output.summary.sdr <- aggregate(Spike.detect.rate ~ Method, data = object$table, FUN = function(x) round(median(x),3))
  
  # Merge
  df <- merge(merge(output.summary.auc,output.summary.fpr, by = "Method"),output.summary.sdr, by = "Method")
  
  # Sort
  if(sort == "AUC") df <- df[order(df$AUC, decreasing = TRUE),]
  if(sort == "FPR") df <- df[order(df$FPR, decreasing = FALSE),]
  if(sort == "Spike.detect.rate") df <- df[order(df$Spike.detect.rate, decreasing = TRUE),]
  
  print(df, row.names = FALSE, ...)
  
}

#' Print results from \code{testDA}
#'
#' @param x The output from the \code{testDA} function
#' @param ... Additional arguments for \code{print}
#' @export

print.DA <- function(x, ...){
  
  xx <- x$table
  print(xx, row.names = FALSE, ...)
  
}
