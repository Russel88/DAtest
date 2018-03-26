#' Summary of results from \code{testDA}
#'
#' @param object The output from the \code{testDA} function
#' @param sort Sort methods by \code{c("AUC","FPR","Spike.detect.rate","Score")}
#' @param boot If TRUE will use bootstrap for confidence limits of the Score, else will compute the limits from the original table. Recommended to be TRUE unless \code{R >= 100} in \code{testDA}
#' @param prob Confidence limits for Score. Default \code{90\%} = \code{c(0.05,0.095)}
#' @param N Number of bootstraps. Default 1000
#' @param boot.seed Random seed for reproducibility of bootstraps
#' @param ... Additional arguments for \code{print}
#' @import stats
#' @import methods
#' @import utils
#' @export

summary.DA <- function(object, sort = "Score", boot = TRUE, prob = c(0.05,0.95), N = 1000, boot.seed = 1, ...){
  
  # Find medians
  output.summary.auc <- aggregate(AUC ~ Method, data = object$table, FUN = function(x) round(median(x),3))
  output.summary.fpr <- aggregate(FPR ~ Method, data = object$table, FUN = function(x) round(median(x),3))
  output.summary.sdr <- aggregate(Spike.detect.rate ~ Method, data = object$table, FUN = function(x) round(median(x),3))
  output.summary.fdr <- aggregate(FDR ~ Method, data = object$table, FUN = function(x) round(median(x),3))
  
  # Merge
  df <- merge(merge(merge(output.summary.auc,output.summary.fpr, by = "Method", all = TRUE),output.summary.fdr, by = "Method"),output.summary.sdr, by = "Method")
  
  # Score
  df$Score <- round(df$AUC * df$Spike.detect.rate - df$FDR,3)
  
  # Interval
  object$table$Score <- object$table$AUC * object$table$Spike.detect.rate - object$table$FDR
  
  if(boot){
    set.seed(boot.seed)
    boots <- lapply(unique(object$table$Method), function(x) object$table[object$table$Method == x,][
      sample(rownames(object$table[object$table$Method == x,]),N,replace = TRUE),
      ])
    boot.score <- lapply(boots,function(y) aggregate(Score ~ Method, data = y, FUN = function(x) round(quantile(x,probs = prob),3)))
    score.cl <- do.call(rbind,boot.score)
  } else {
    score.cl <- aggregate(Score ~ Method, data = object$table, FUN = function(x) round(quantile(x,probs = prob),3))
  }
 
  df <- merge(df, as.matrix(score.cl), by = "Method")
  
  # Sort
  if(sort == "AUC") df <- df[order(df$AUC, decreasing = TRUE),]
  if(sort == "FPR") df <- df[order(df$FPR, decreasing = FALSE),]
  if(sort == "Spike.detect.rate") df <- df[order(df$Spike.detect.rate, decreasing = TRUE),]
  if(sort == "Score") df <- df[order(df$Score, decreasing = TRUE),]
  
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
