#' Summary of results from testDA
#'
#' @param DA The output from the testDA function
#' @param sort Sort methods by c("AUC","FPR","Spike.detect.rate")
#' @export

summary.DA <- function(DA, sort = "AUC"){
  
  output.summary.auc <- aggregate(AUC ~ Method, data = DA$table, FUN = function(x) round(median(x),3))
  output.summary.fpr <- aggregate(FPR ~ Method, data = DA$table, FUN = function(x) round(median(x),3))
  output.summary.sdr <- aggregate(Spike.detect.rate ~ Method, data = DA$table, FUN = function(x) round(median(x),3))
  
  df <- merge(merge(output.summary.auc,output.summary.fpr, by = "Method"),output.summary.sdr, by = "Method")
  
  if(sort == "AUC") df <- df[order(df$AUC, decreasing = TRUE),]
  if(sort == "FPR") df <- df[order(df$FPR, decreasing = FALSE),]
  if(sort == "Spike.detect.rate") df <- df[order(df$Spike.detect.rate, decreasing = TRUE),]
  
  print(df, row.names = FALSE)
  
}

#' Print results from testDA
#'
#' @param DA The output from the testDA function
#' @export

print.DA <- function(DA){
  
  print(DA$table, row.names = FALSE)
  
}
