#' Plotting results from \code{powerDA}
#'
#' @param x The output from the \code{powerDA} function
#' @param ... Additional arguments for \code{ggdraw}
#' @import ggplot2
#' @importFrom cowplot ggdraw
#' @importFrom cowplot draw_plot
#' @export
plot.DAPower <- function(x, ...){
  
  z <- as.data.frame(unclass(x[[1]]))
  z$EffectSize <- log2(z$EffectSize)
    
  # Define plots
  p1 <- ggplot(z, aes_string(x = "EffectSize", y = "Power")) +
    theme_bw() +
    coord_cartesian(ylim = c(0,1)) +
    geom_point() +
    geom_smooth(method = "loess", colour = "red", alpha = 0.4) +
    ylab("Power") +
    xlab("Log2 Effect Size") +
    theme(panel.grid.minor = element_blank())
    
  p2 <- ggplot(z, aes_string(x = "EffectSize", y = "AUC")) +
    theme_bw() +
    coord_cartesian(ylim = c(min(0.45, min(z$AUC)),1)) +
    geom_point() +
    geom_hline(yintercept = 0.5) +
    geom_smooth(method = "loess", colour = "red", alpha = 0.4) +
    ylab("Area Under the Curve") +
    xlab("Log2 Effect Size") +
    theme(panel.grid.minor = element_blank())
    
  p3 <- ggplot(z, aes_string(x = "EffectSize", y = "FPR")) +
    theme_bw() +
    coord_cartesian(ylim = c(0,(max(z$FPR)+0.05))) +
    geom_point() +
    geom_hline(yintercept = x[[2]]) +
    geom_smooth(method = "loess", colour = "red", alpha = 0.4) +
    ylab("False Positive Rate") +
    xlab("Log2 Effect Size") +
    theme(panel.grid.minor = element_blank())
    
  p4 <- ggplot(z, aes_string(x = "EffectSize", y = "FDR")) +
    theme_bw() +
    coord_cartesian(ylim = c(0,(max(z$FDR)+0.05))) +
    geom_point() +
    geom_hline(yintercept = x[[3]]) +
    geom_smooth(method = "loess", colour = "red", alpha = 0.4) +
    ylab("False Discovery Rate") +
    xlab("Log2 Effect Size") +
    theme(panel.grid.minor = element_blank())
   
  if(z$Method[1] %in% c("ALDEx2 t-test (adx)","ALDEx2 wilcox (adx)")){
    p1 <- p1 + facet_grid(.~Method)
    p2 <- p2 + facet_grid(.~Method)
    p3 <- p3 + facet_grid(.~Method)
    p4 <- p4 + facet_grid(.~Method)
  }
   
  if(all(is.na(z$Power))) p1 <- NULL
  if(all(is.na(z$AUC))) p2 <- NULL
  if(all(is.na(z$FPR))) p3 <- NULL
  if(all(is.na(z$FDR))) p4 <- NULL
   
  # Plot it
  pp <- cowplot::ggdraw(...) 
  if(!is.null(p1)) pp <- pp + cowplot::draw_plot(p1, 0, 0.5, 0.48, 0.48)
  if(!is.null(p2)) pp <- pp + cowplot::draw_plot(p2, 0.5, 0.5, 0.48, 0.48) 
  if(!is.null(p3)) pp <- pp + cowplot::draw_plot(p3, 0, 0, 0.48, 0.48) 
  if(!is.null(p4)) pp <- pp + cowplot::draw_plot(p4, 0.5, 0, 0.48, 0.48)
  
  suppressWarnings(pp)
    
}
