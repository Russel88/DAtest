#' Summary of results from \code{powerDA}
#'
#' @param x The output from the \code{powerDA} function
#' @param ... Additional arguments for \code{print}
#' @export
summary.DAPower <- function(x, ...){
    
    x <- as.data.frame(unclass(x[[1]]))
    
    if(x$Method[1] == "SAMseq (sam)"){
      # Find medians
      output.summary.sdr <- aggregate(Power ~ EffectSize, data = x, FUN = function(y) round(median(y),3))
      output.summary.auc <- aggregate(AUC ~ EffectSize, data = x, FUN = function(y) round(median(y),3))
      output.summary.fdr <- aggregate(FDR ~ EffectSize, data = x, FUN = function(y) round(median(y),3))
      
      # Merge
      df <- merge(merge(output.summary.sdr,output.summary.auc, by = "EffectSize"),output.summary.fdr, by = "EffectSize")
      colnames(df) <- c("EffectSize","Power","AUC","FDR")
    } else {
      if(x$Method[1] %in% c("ALDEx2 t-test (adx)","ALDEx2 wilcox (adx)")){
        # Find medians
        output.summary.sdr <- aggregate(Power ~ Method + EffectSize, data = x, FUN = function(y) round(median(y),3))
        output.summary.auc <- aggregate(AUC ~ Method + EffectSize, data = x, FUN = function(y) round(median(y),3))
        output.summary.fpr <- aggregate(FPR ~ Method + EffectSize, data = x, FUN = function(y) round(median(y),3))
        output.summary.fdr <- aggregate(FDR ~ Method + EffectSize, data = x, FUN = function(y) round(median(y),3))
        
        # Merge
        df <- cbind(output.summary.sdr,output.summary.auc[,3],output.summary.fpr[,3],output.summary.fdr[,3])
        df <- as.data.frame(df[order(df$Method),])
        colnames(df) <- c("Method","EffectSize","Power","AUC","FPR","FDR")
      } else {
        # Find medians
        output.summary.sdr <- aggregate(Power ~ EffectSize, data = x, FUN = function(y) round(median(y),3))
        output.summary.auc <- aggregate(AUC ~ EffectSize, data = x, FUN = function(y) round(median(y),3))
        output.summary.fpr <- aggregate(FPR ~ EffectSize, data = x, FUN = function(y) round(median(y),3))
        output.summary.fdr <- aggregate(FDR ~ EffectSize, data = x, FUN = function(y) round(median(y),3))
        
        # Merge
        df <- merge(merge(merge(output.summary.sdr,output.summary.auc, by = "EffectSize"),output.summary.fpr, by = "EffectSize"),output.summary.fdr, by = "EffectSize")
        colnames(df) <- c("EffectSize","Power","AUC","FPR","FDR")
      }
    } 
  
    print(df, row.names = FALSE, ...)
    
}
