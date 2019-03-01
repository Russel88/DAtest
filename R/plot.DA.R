#' Plotting results from \code{testDA}
#'
#' @param x The output from the \code{testDA} function
#' @param sort Sort methods by median \code{c("AUC","FDR","Power","Score")}
#' @param p Logical. Should the p-value distribution be plotted (only p-values from non-spiked features)
#' @param bins Integer. Number of bins in p-value histograms
#' @param ... Additional arguments for \code{ggdraw}
#' @return Plots the output from testDA
#' @import ggplot2
#' @importFrom cowplot ggdraw
#' @importFrom cowplot draw_plot
#' @importFrom graphics pairs plot
#' @export

plot.DA <- function(x, sort = "Score", p = FALSE, bins = 20, ...){
  
  if(p){
    # For plotting p-value histograms
    ## Extract p-values and remove all Spiked
    pval.all <- lapply(x$results, function(x) lapply(x, function(y) y[,c("pval","Method","Spiked")]))
    df.all <- suppressWarnings(do.call(rbind, do.call(rbind,pval.all)))
    df.all <- df.all[df.all$Spiked == "No",]
    df.all <- df.all[!df.all$Method %in% c("ANCOM (anc)","SAMseq (sam)"),]
    
    ## Plot
    ggplot(df.all, aes_string(x = "pval")) +
      theme_bw() +
      geom_histogram(bins=bins) +
      facet_wrap(~Method, scales = "free_y") +
      ylab("Density") +
      xlab("P-value")
    
  } else {
    
    # Score
    ## Find medians
    auc.median <- aggregate(AUC ~ Method, data = x$table, FUN = function(x) round(median(x),3))
    fdr.median <- aggregate(FDR ~ Method, data = x$table, FUN = function(x) round(median(x),3))
    sdr.median <- aggregate(Power ~ Method, data = x$table, FUN = function(x) round(median(x),3))
    
    ## Merge
    df <- merge(merge(auc.median,fdr.median, by = "Method"),sdr.median, by = "Method")
    
    # Score
    df$Score <- round((df$AUC-0.5) * df$Power - df$FDR,3)
    
    # Sort the reults
    if(sort == "AUC") {
      x$table$Method <- factor(x$table$Method, levels = df[order(df$AUC, decreasing = TRUE),"Method"])
    }
    if(sort == "FDR") {
      x$table$Method <- factor(x$table$Method, levels = df[order(df$FDR, decreasing = FALSE),"Method"])
    }
    if(sort == "Power") {
      x$table$Method <- factor(x$table$Method, levels = df[order(df$Power, decreasing = TRUE),"Method"])
    }
    if(sort == "Score") {
      x$table$Method <- factor(x$table$Method, levels = df[order(df$Score, decreasing = TRUE),"Method"])
    }
    
    # Define FDR and AUC plots
    p1 <- ggplot(x$table, aes_string(x = "Method", y = "FDR")) +
      theme_bw() +
      coord_cartesian(ylim = c(0,1)) +
      geom_point() +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",colour="red",width=0.75) +
      ylab("False Discovery Rate") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      xlab(NULL) +
      theme(panel.grid.minor = element_blank())
    
    p2 <- ggplot(x$table, aes_string(x = "Method", y = "AUC")) +
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
    
    p3 <- ggplot(x$table, aes_string(x = "Method", y = "Power")) +
      theme_bw() +
      geom_point() +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",colour="red",width=0.75) +
      ylab("Power") +
      theme(axis.text.x = element_blank(),
            panel.grid.minor = element_blank()) +
      xlab(NULL) +
      scale_y_continuous(labels=function(x) sprintf("%.2f", x))
    
    # Plot it
    pp <- cowplot::ggdraw(...) +
      cowplot::draw_plot(p2, 0, .70, 1, .30) +
      cowplot::draw_plot(p3, 0, .40, 1, .30) +
      cowplot::draw_plot(p1, 0, 0, 1, .40)
    
    suppressWarnings(pp)
    
  }
  
}
