#' ANCOM-BC
#'
#' Implementation of ANCOM-BC for \code{DAtest}
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param covars Either a named list with covariables, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param out.all If TRUE, will run global test which will produce one p-value for the \code{predictor}. If FALSE will run standard test and will output p-value from one level of the predictor specified by \code{coeff}. If NULL (default) set as TRUE for multi-class \code{predictor} and FALSE otherwise
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param coeff Integer. The beta coefficient and p-value will be associated with this coefficient. This coefficient is by default compared to the intercept (1. level of \code{predictor}) Default 2, i.e. the 2. level of the \code{predictor}.
#' @param allResults If TRUE will return raw results from the \code{ancombc} function
#' @param ... Additional arguments for the \code{ancombc} function
#' @return A data.frame with with results.
#' @examples
#' # Creating random count_table and predictor
#' set.seed(4)
#' mat <- matrix(rnbinom(200, size = 0.1, mu = 500), nrow = 20, ncol = 10)
#' rownames(mat) <- 1:20
#' pred <- c(rep("Control", 5), rep("Treatment", 5))
#' 
#' # Running ANCOM-BC
#' res <- DA.abc(data = mat, predictor = pred)
#' @export
DA.abc <- function(data, predictor, covars = NULL, out.all = NULL, p.adj = "fdr", coeff = 2, allResults = FALSE, ...){
  
  ok1 <- tryCatch({
    loadNamespace("ANCOMBC")
    TRUE
  }, error=function(...) FALSE)
  ok2 <- tryCatch({
    loadNamespace("phyloseq")
    TRUE
  }, error=function(...) FALSE)
  ok <- ok1 && ok2
  
  if (ok){
    # Convert to phyloseq
    if(is(data, "phyloseq")){
      phy_data <- data
      org_pred <- unlist(phyloseq::sample_data(data)[, predictor])
      if(!is.null(covars)){
        form <- paste(predictor, paste(covars, collapse="+"), sep="+")
      } else {
        form <- predictor
      }
      
    } else {
      org_pred <- predictor
      predictor <- "predictor"
      otu_table <-  phyloseq::otu_table(data, taxa_are_rows = TRUE)
      if(!is.null(covars)){
        samp_table <- phyloseq::sample_data(data.frame(row.names = colnames(data), predictor = org_pred))
        form <- paste("predictor", paste(names(covars), collapse="+"), sep="+")
        for(covar_sub in names(covars)){
          samp_table[, covar_sub] <- covars[covar_sub]
        }
      } else {
        samp_table <-  phyloseq::sample_data(data.frame(row.names = colnames(data), predictor = org_pred))
        form <- "predictor"
      }
      phy_data <- phyloseq::phyloseq(otu_table, samp_table)
    }
    
    coeff.ref <- 1
    if(coeff == coeff.ref) stop("coeff and coeff.ref cannot be the same")
    
    # out.all
    if(is.null(out.all)){
      if(length(unique(org_pred)) == 2) out.all <- FALSE
      if(length(unique(org_pred)) > 2) out.all <- TRUE
    }
    
    # Run test
    if(out.all){
      abc_res <- ANCOMBC::ancombc(phy_data, formula = form, p_adj_method = p.adj, group = predictor, global = TRUE, ...)
      
      res <- data.frame(Feature = rownames(abc_res$res_global),
                        Beta = NA,
                        pval = abc_res$res_global$p_val,
                        pval.adj = abc_res$res_global$q_val)
      
    }
    if(!out.all){
      abc_res <- ANCOMBC::ancombc(phy_data, formula = form, p_adj_method = p.adj, global = FALSE, ...)
      
      res <- data.frame(Feature = rownames(abc_res$res$beta),
                        Beta = as.numeric(abc_res$res$beta[, coeff - 1]),
                        pval = as.numeric(abc_res$res$p_val[, coeff - 1]),
                        pval.adj = as.numeric(abc_res$res$q_val[, coeff - 1]))
     
    }
    
    if(!is.numeric(predictor)){
      res$ordering <- NA
      res[!is.na(res$Beta) & res$Beta > 0,"ordering"] <- paste0(levels(as.factor(org_pred))[coeff],">",levels(as.factor(org_pred))[coeff.ref])
      res[!is.na(res$Beta) & res$Beta < 0,"ordering"] <- paste0(levels(as.factor(org_pred))[coeff.ref],">",levels(as.factor(org_pred))[coeff])
    }
    
    res$Method <- "ANCOM-BC (abc)"
    
    if(is(data, "phyloseq")) res <- addTax(data, res)
    
    if(allResults) return(abc_res) else return(res)
  } else {
    stop("ANCOM-BC package required")
  }
    
}

