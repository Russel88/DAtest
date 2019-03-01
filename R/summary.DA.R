#' Summary of results from \code{testDA}
#'
#' @param object The output from the \code{testDA} function
#' @param sort Sort methods by \code{c("AUC","FPR","Power","Score")}
#' @param boot If TRUE will use bootstrap for confidence limits of the Score, else will compute the limits from the original table. Recommended to be TRUE unless \code{R >= 100} in \code{testDA}
#' @param prob Confidence limits for Score. Default \code{90\%} = \code{c(0.05,0.095)}
#' @param N Number of bootstraps. Default 1000
#' @param decimals Precision in output. Default 2
#' @param ... Additional arguments for \code{print}
#' @return Prints summary from the DAtest function. '*' next to method, means that it's median score is overlapping with the 90% confidence limits of the best method
#' @import stats
#' @import methods
#' @import utils
#' @export

summary.DA <- function(object, sort = "Score", boot = TRUE, prob = c(0.05,0.95), N = 1000, decimals = 2, ...){
  
  # Find medians
  output.summary.auc <- aggregate(AUC ~ Method, data = object$table, FUN = median)
  output.summary.fpr <- aggregate(FPR ~ Method, data = object$table, FUN = median)
  output.summary.sdr <- aggregate(Power ~ Method, data = object$table, FUN = median)
  output.summary.fdr <- aggregate(FDR ~ Method, data = object$table, FUN = median)
  
  # Merge
  df <- merge(merge(merge(output.summary.auc,output.summary.fpr, by = "Method", all = TRUE),output.summary.fdr, by = "Method"),output.summary.sdr, by = "Method")
  
  # Score
  df$Score <- (df$AUC-0.5) * df$Power - df$FDR
  
  # Interval
  object$table$Score <- (object$table$AUC-0.5) * object$table$Power - object$table$FDR
  
  if(boot){
    boots <- lapply(unique(object$table$Method), function(x) object$table[object$table$Method == x,][
      sample(rownames(object$table[object$table$Method == x,]),N,replace = TRUE),
      ])
    boot.score <- lapply(boots,function(y) aggregate(Score ~ Method, data = y, FUN = function(x) quantile(x,probs = prob)))
    score.cl <- do.call(rbind,boot.score)
  } else {
    score.cl <- aggregate(Score ~ Method, data = object$table, FUN = function(x) quantile(x,probs = prob))
  }
 
  df <- merge(df, as.matrix(score.cl), by = "Method")
  df <- df[order(df$Score,df[,7],df[,8], decreasing = TRUE),]
  if(df[1,"Score"] == 0) warning("Best Score is equal to zero!\nYou might want to re-run with a higher effectSize or pruned dataset (see preDA)")
  
  df <- cbind(data.frame(Method = df[,1]),apply(df[,2:ncol(df)],2,function(x) round(as.numeric(x), decimals)))
  
  mat <- c(FALSE)
  for(i in seq_len(nrow(df))){
    mat[i] <- df[i,6] >= df[1,7]
  }
  df$` ` <- " "
  df[mat,]$` ` <- "*" 
  
  # Sort
  if(sort == "AUC") df <- df[order(df$AUC, decreasing = TRUE),]
  if(sort == "FPR") df <- df[order(df$FPR, decreasing = FALSE),]
  if(sort == "Power") df <- df[order(df$Power, decreasing = TRUE),]

  print(df, row.names = FALSE, ...)
  
}

