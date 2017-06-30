#' LIMMA
#'
#' Some is borrowed from:
#' http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html

#' @export

DA.lim <- function(count_table, outcome, p.adj, relative, log = FALSE, delta = 1){
  
  library(limma, quietly = TRUE)
  
  if(log) count_table <- log(count_table + delta)
  
  if(relative){
    count.rel <- apply(count_table,2,function(x) x/sum(x))
  } else {
    count.rel <- count_table
  }

  design <- model.matrix(~outcome)
  n <- dim(count.rel)[1]
  fit <- lmFit(count.rel, design)
  fit.eb <- eBayes(fit)
  Estimate <- fit.eb$coefficients[, 2]
  df.residual <- fit.eb$df.residual
  df.prior <- rep(fit.eb$df.prior, n)
  s2.prior <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.stat <- fit.eb$t[, 2]
  pval <- fit.eb$p.value[, 2]
  pval.adj <- p.adjust(pval, method = p.adj)
  res <- data.frame(Estimate, t.stat, pval, pval.adj, df.residual, df.prior, s2.prior, s2, s2.post)
  res$Feature <- rownames(res)
  res$Method <- "LIMMA"

  return(res)  
}

