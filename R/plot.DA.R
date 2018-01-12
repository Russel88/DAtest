#' Plotting results from \code{testDA}
#'
#' @param x The output from the \code{testDA} function
#' @param sort Sort methods by \code{c("AUC","FPR","Spike.detect.rate","Rank")}
#' @param p Logical. Should the p-value distribution be plotted (only p-values from non-spiked features)
#' @param bins Integer. Number of bins in p-value histograms
#' @param ... Additional arguments for \code{ggdraw}
#' @import ggplot2
#' @importFrom cowplot ggdraw
#' @importFrom cowplot draw_plot
#' @export

plot.DA <- function(x, sort = "Rank", p = FALSE, bins = 20, ...){
  
  if(p){
    # For plotting p-value histograms
    ## Extract p-values and remove all Spiked
    pval.all <- lapply(x$results, function(x) lapply(x, function(y) y[,c("pval","Method","Spiked")]))
    df.all <- suppressWarnings(do.call(rbind, do.call(rbind,pval.all)))
    df.all <- df.all[df.all$Spiked == "No",]
    df.all <- df.all[!df.all$Method %in% c("ANCOM (anc)","SAMseq (sam)"),]
    
    ## Plot
    ggplot(df.all, aes(pval)) +
      theme_bw() +
      geom_histogram(bins=bins) +
      facet_wrap(~Method, scales = "free_y") +
      ylab("Density") +
      xlab("P-value")
    
  } else {
    
    # Rank
    ## Find medians
    output.summary.auc <- aggregate(AUC ~ Method, data = x$table, FUN = function(x) round(median(x),3))
    output.summary.fpr <- aggregate(FPR ~ Method, data = x$table, FUN = function(x) round(median(x),3))
    output.summary.sdr <- aggregate(Spike.detect.rate ~ Method, data = x$table, FUN = function(x) round(median(x),3))
    
    ## Merge
    df <- merge(merge(output.summary.auc,output.summary.fpr, by = "Method"),output.summary.sdr, by = "Method")
    
    ## Rank
    nobad <- df[df$FPR <= 0.05,]
    auc.r <- rank(nobad$AUC)
    sdr.r <- rank(nobad$Spike.detect.rate)
    mean.r <- (auc.r + sdr.r)/2
    nobad$Rank <- nrow(nobad) - mean.r
    bad <- df[df$FPR > 0.05,]
    bad <- bad[order(bad$AUC, decreasing = TRUE),]
    bad$Rank <- NA
    df <- rbind(nobad,bad)
    df <- df[order(df$Rank, decreasing = FALSE),]
    
    # Sort the reults
    if(sort == "AUC") {
      auc.median <- aggregate(AUC ~ Method, data = x$table, FUN = median)
      x$table$Method <- factor(x$table$Method, levels = auc.median[order(auc.median$AUC, decreasing = TRUE),"Method"])
    }
    if(sort == "FPR") {
      fpr.median <- aggregate(FPR ~ Method, data = x$table, FUN = median)
      x$table$Method <- factor(x$table$Method, levels = fpr.median[order(fpr.median$FPR, decreasing = FALSE),"Method"])
    }
    if(sort == "Spike.detect.rate") {
      sdr.median <- aggregate(Spike.detect.rate ~ Method, data = x$table, FUN = median)
      x$table$Method <- factor(x$table$Method, levels = sdr.median[order(sdr.median$Spike.detect.rate, decreasing = TRUE),"Method"])
    }
    if(sort == "Rank") {
      x$table$Method <- factor(x$table$Method, levels = df$Method)
    }
    
    # Colour bad ones
    cob <- data.frame(x = rep(c((nrow(nobad)+0.5),(nrow(df)+1)),2),
                      ymin = rep(-Inf,4),
                      ymax = rep(Inf,4),
                      FPR = rep(0,4),
                      AUC = rep(0,4),
                      Spike.detect.rate = rep(0,4))
    
    # Define FPR and AUC plots
    p1 <- ggplot(x$table, aes(Method, FPR)) +
      theme_bw() +
      coord_cartesian(ylim = c(0,1)) +
      geom_hline(yintercept = 0.05, colour = "red") +
      geom_point() +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",colour="red",width=0.75) +
      ylab("False Positive Rate") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      xlab(NULL) +
      theme(panel.grid.minor = element_blank()) +
      geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax),data = cob, fill = "red", alpha = 0.2)
    
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
      scale_y_continuous(labels=function(x) sprintf("%.2f", x)) +
      geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax),data = cob, fill = "red", alpha = 0.2)
    
    p3 <- ggplot(x$table, aes(Method, Spike.detect.rate)) +
      theme_bw() +
      geom_point() +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",colour="red",width=0.75) +
      ylab("Spike detect rate") +
      theme(axis.text.x = element_blank(),
            panel.grid.minor = element_blank()) +
      xlab(NULL) +
      scale_y_continuous(labels=function(x) sprintf("%.2f", x)) +
      geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax),data = cob, fill = "red", alpha = 0.2)
    
    # Plot it
    pp <- cowplot::ggdraw(...) +
      cowplot::draw_plot(p2, 0, .70, 1, .30) +
      cowplot::draw_plot(p3, 0, .40, 1, .30) +
      cowplot::draw_plot(p1, 0, 0, 1, .40)
    
    suppressWarnings(pp)
    
  }
  
}