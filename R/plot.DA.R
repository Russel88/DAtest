#' Plotting results from testDA
#'
#' @param DA The output from the testDA function
#' @param sort Sort methods by c("AUC","FPR","Spike.detect.rate")
#' @param p Logical. Should the p-value distribution be plotted
#' @param bins Integer. Number of bins in p-value histograms
#' @param adj Logical. Whether the histograms should show adjusted p-values 
#' @import ggplot2 cowplot
#' 
#' @export

plot.DA <- function(DA, sort = "AUC", p = FALSE, bins = 50, adj = FALSE){
  
  require(ggplot2, quietly = TRUE)
  require(cowplot, quietly = TRUE)
  
  if(p){
    pval.all <- lapply(DA$results, function(x) lapply(x, function(y) y[,c("pval","pval.adj","Method")]))
    df.all <- do.call(rbind, do.call(rbind,pval.all))
    
    if(adj){
      ggplot(df.all, aes(pval.adj)) +
        theme_bw() +
        geom_histogram(bins=50) +
        facet_wrap(~Method, scales = "free_y") +
        ylab("Density") +
        xlab("Adjusted P-value")
      
    } else {
      ggplot(df.all, aes(pval)) +
        theme_bw() +
        geom_histogram(bins=50) +
        facet_wrap(~Method, scales = "free_y") +
        ylab("Density") +
        xlab("P-value")
    }
    
  } else {
    
    if(sort == "AUC") {
      auc.median <- aggregate(AUC ~ Method, data = DA$table, FUN = median)
      DA$table$Method <- factor(DA$table$Method, levels = auc.median[order(auc.median$AUC, decreasing = TRUE),"Method"])
    }
    if(sort == "FPR") {
      fpr.median <- aggregate(FPR ~ Method, data = DA$table, FUN = median)
      DA$table$Method <- factor(DA$table$Method, levels = fpr.median[order(fpr.median$FPR, decreasing = FALSE),"Method"])
    }
    if(sort == "Spike.detect.rate") {
      spr.median <- aggregate(Spike.detect.rate ~ Method, data = DA$table, FUN = median)
      DA$table$Method <- factor(DA$table$Method, levels = fpr.median[order(fpr.median$Spike.detect.rate, decreasing = TRUE),"Method"])
    }
    
    p1 <- ggplot(DA$table, aes(Method, FPR)) +
      theme_bw() +
      coord_cartesian(ylim = c(0,1)) +
      geom_hline(yintercept = 0.05, colour = "red") +
      geom_point() +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",colour="red",width=0.75) +
      ylab("False Positive Rate") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      xlab(NULL) +
      theme(axis.text.x = element_blank(),
            panel.grid.minor = element_blank())
    
    p2 <- ggplot(DA$table, aes(Method, AUC)) +
      theme_bw() +
      coord_cartesian(ylim = c(min(DA$table$AUC, 0.45),1)) +
      geom_hline(yintercept = 0.5, colour = "red") +
      geom_point() +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",colour="red",width=0.75) +
      ylab("Area Under the Curve") +
      theme(axis.text.x = element_blank(),
            panel.grid.minor = element_blank()) +
      xlab(NULL) +
      scale_y_continuous(labels=function(x) sprintf("%.2f", x))
    
    p3 <- ggplot(DA$table, aes(Method, Spike.detect.rate)) +
      theme_bw() +
      coord_cartesian(ylim = c(0,1)) +
      geom_point() +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",colour="red",width=0.75) +
      ylab("Spike-In Detection Rate") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            panel.grid.minor = element_blank()) +
      xlab(NULL)
    
    pp <- ggdraw() +
      draw_plot(p1, 0, .7, 1, .3) +
      draw_plot(p2, 0, .4, 1, .3) +
      draw_plot(p3, 0, 0, 1, .4)
    
    suppressWarnings(pp)
    
  }
  
}