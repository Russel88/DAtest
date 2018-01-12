#' Summary of results from \code{testDA}
#'
#' @param object The output from the \code{testDA} function
#' @param sort Sort methods by \code{c("AUC","FPR","Spike.detect.rate","Rank")}
#' @param ... Additional arguments for \code{print}
#' @import stats
#' @import methods
#' @import utils
#' @export

summary.DA <- function(object, sort = "Rank", ...){
  
  # Find medians
  output.summary.auc <- aggregate(AUC ~ Method, data = object$table, FUN = function(x) round(median(x),3))
  output.summary.fpr <- aggregate(FPR ~ Method, data = object$table, FUN = function(x) round(median(x),3))
  output.summary.sdr <- aggregate(Spike.detect.rate ~ Method, data = object$table, FUN = function(x) round(median(x),3))
  
  # Merge
  df <- merge(merge(output.summary.auc,output.summary.fpr, by = "Method"),output.summary.sdr, by = "Method")
  
  # Rank
  nobad <- df[df$FPR <= 0.05,]
  auc.r <- rank(nobad$AUC)
  sdr.r <- rank(nobad$Spike.detect.rate)
  mean.r <- (auc.r + sdr.r)/2
  nobad$Rank <- nrow(nobad) - mean.r
  bad <- df[df$FPR > 0.05,]
  bad <- bad[order(bad$AUC, decreasing = TRUE),]
  bad$Rank <- NA
  df <- rbind(nobad,bad)
  
  # Sort
  if(sort == "AUC") df <- df[order(df$AUC, decreasing = TRUE),]
  if(sort == "FPR") df <- df[order(df$FPR, decreasing = FALSE),]
  if(sort == "Spike.detect.rate") df <- df[order(df$Spike.detect.rate, decreasing = TRUE),]
  if(sort == "Rank") df <- df[order(df$Rank, decreasing = FALSE),]
  
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
