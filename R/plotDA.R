#' Plotting results from testDA
#'
#' @param DA The output from the testDA function
#' @import ggplot2 cowplot
#' 
#' @export

plotDA <- function(DA){
  
  require(ggplot2, quietly = TRUE)
  require(cowplot, quietly = TRUE)
  
  DA$table$Method <- factor(DA$table$Method, levels = DA$summary[order(DA$summary$AUC, decreasing = TRUE),"Method"])
  
  p1 <- ggplot(DA$table, aes(Method, FPR)) +
    theme_bw() +
    coord_cartesian(ylim = c(0,1)) +
    geom_hline(yintercept = 0.05, colour = "red") +
    geom_point() +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",colour="red") +
    ylab("False Positive Rate") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab(NULL) +
    theme(axis.text.x = element_blank(),
          panel.grid.minor = element_blank())
  
  p2 <- ggplot(DA$table, aes(Method, AUC)) +
    theme_bw() +
    coord_cartesian(ylim = c(min(DA$table$AUC),1)) +
    geom_hline(yintercept = 0.5, colour = "red") +
    geom_point() +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,geom = "crossbar",colour="red") +
    ylab("Area Under the Curve") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.minor = element_blank()) +
    xlab(NULL)
  
  suppressWarnings(plot_grid(p1, p2, nrow=2))
  
}