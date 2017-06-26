#' Summary and print of results from testDA
#'
#' @param DA The output from the testDA function
#' 
#' @export

summary.DA <- function(DA){
  
  output.summary.auc <- aggregate(AUC ~ Method, data = DA$table, FUN = function(x) round(median(x),3))
  output.summary.fpr <- aggregate(FPR ~ Method, data = DA$table, FUN = function(x) round(median(x),3))
  output.summary.sdr <- aggregate(Spike.detect.rate ~ Method, data = DA$table, FUN = function(x) round(median(x),3))
  
  df <- merge(merge(output.summary.auc,output.summary.fpr, by = "Method"),output.summary.sdr, by = "Method")
  
  df <- df[order(df$AUC, decreasing = TRUE),]
  
  print(df, row.names = FALSE)
  
}

#' @export

print.DA <- function(DA){
  
  print(DA$table, row.names = FALSE)
  
}
