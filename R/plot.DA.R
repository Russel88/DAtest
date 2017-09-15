#' Plotting results from testDA
#'
#' @param x The output from the testDA function
#' @param sort Sort methods by c("AUC","FPR")
#' @param p Logical. Should the p-value distribution be plotted (only p-values from non-spiked features)
#' @param bins Integer. Number of bins in p-value histograms
#' @param ... Additional plotting arguments
#' @import ggplot2
#' @importFrom cowplot ggdraw
#' @importFrom cowplot draw_plot
#' @export

plot.DA <- function(x, sort = "AUC", p = FALSE, bins = 20, ...){
  
  if(p){
    pval.all <- lapply(x$results, function(x) lapply(x, function(y) y[,c("pval","Method","Spiked")]))
    df.all <- do.call(rbind, do.call(rbind,pval.all))
    df.all <- df.all[df.all$Spiked == "No",]
    
    ggplot(df.all, aes(pval)) +
      theme_bw() +
      geom_histogram(bins=bins) +
      facet_wrap(~Method, scales = "free_y") +
      ylab("Density") +
      xlab("P-value")
    
  } else {
    
    if(sort == "AUC") {
      auc.median <- aggregate(AUC ~ Method, data = x$table, FUN = median)
      x$table$Method <- factor(x$table$Method, levels = auc.median[order(auc.median$AUC, decreasing = TRUE),"Method"])
    }
    if(sort == "FPR") {
      fpr.median <- aggregate(FPR ~ Method, data = x$table, FUN = median)
      x$table$Method <- factor(x$table$Method, levels = fpr.median[order(fpr.median$FPR, decreasing = FALSE),"Method"])
    }
    
    p1 <- ggplot(x$table, aes(Method, FPR)) +
      theme_bw() +
      coord_cartesian(ylim = c(0,1)) +
      geom_hline(yintercept = 0.05, colour = "red") +
      geom_point() +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",colour="red",width=0.75) +
      ylab("False Positive Rate") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      xlab(NULL) +
      theme(panel.grid.minor = element_blank())
    
    p2 <- ggplot(x$table, aes(Method, AUC)) +
      theme_bw() +
      coord_cartesian(ylim = c(min(x$table$AUC, 0.45),1)) +
      geom_hline(yintercept = 0.5, colour = "red") +
      geom_point() +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",colour="red",width=0.75) +
      ylab("Area Under the Curve") +
      theme(axis.text.x = element_blank(),
            panel.grid.minor = element_blank()) +
      xlab(NULL) +
      scale_y_continuous(labels=function(x) sprintf("%.2f", x))
    
    pp <- cowplot::ggdraw(...) +
      cowplot::draw_plot(p2, 0, .5, 1, .5) +
      cowplot::draw_plot(p1, 0, 0, 1, .5)
    
    suppressWarnings(pp)
    
  }
  
}