#' Summary of results from \code{powerDA}
#'
#' @param object The output from the \code{powerDA} function
#' @param ... Additional arguments for \code{print}
#' @param decimals Precision in output. Default 2
#' @return Prints summary results from the powerDA function
#' @export
summary.DAPower <- function(object, decimals = 2, ...){
    
    x <- as.data.frame(unclass(object[[1]]))
    
    if(x$Method[1] == "SAMseq (sam)"){
      # Find medians
      output.summary.sdr <- aggregate(Power ~ EffectSize, data = x, FUN = function(y) median)
      output.summary.auc <- aggregate(AUC ~ EffectSize, data = x, FUN = function(y) median)
      output.summary.fdr <- aggregate(FDR ~ EffectSize, data = x, FUN = function(y) median)
      
      # Merge
      df <- merge(merge(output.summary.sdr,output.summary.auc, by = "EffectSize"),output.summary.fdr, by = "EffectSize")
      colnames(df) <- c("EffectSize","Power","AUC","FDR")
    } else {
      if(x$Method[1] %in% c("ALDEx2 t-test (adx)","ALDEx2 wilcox (adx)")){
        # Find medians
        output.summary.sdr <- aggregate(Power ~ Method + EffectSize, data = x, FUN = median)
        output.summary.auc <- aggregate(AUC ~ Method + EffectSize, data = x, FUN = median)
        output.summary.fpr <- aggregate(FPR ~ Method + EffectSize, data = x, FUN = median)
        output.summary.fdr <- aggregate(FDR ~ Method + EffectSize, data = x, FUN = median)
        
        # Merge
        df <- cbind(output.summary.sdr,output.summary.auc[,3],output.summary.fpr[,3],output.summary.fdr[,3])
        df <- as.data.frame(df[order(df$Method),])
        colnames(df) <- c("Method","EffectSize","Power","AUC","FPR","FDR")
      } else {
        # Find medians
        output.summary.sdr <- aggregate(Power ~ EffectSize, data = x, FUN = median)
        output.summary.auc <- aggregate(AUC ~ EffectSize, data = x, FUN = median)
        output.summary.fpr <- aggregate(FPR ~ EffectSize, data = x, FUN = median)
        output.summary.fdr <- aggregate(FDR ~ EffectSize, data = x, FUN = median)
        
        # Merge
        df <- merge(merge(merge(output.summary.sdr,output.summary.auc, by = "EffectSize"),
                          output.summary.fpr, by = "EffectSize"),
                    output.summary.fdr, by = "EffectSize")
        colnames(df) <- c("EffectSize","Power","AUC","FPR","FDR")
      }
    } 
  
    df <- cbind(data.frame(EffectSize = df[,1]),
                apply(df[, 2:ncol(df)], 2, function(x) round(as.numeric(x), decimals)))

    print(df, row.names = FALSE, ...)
    
}
