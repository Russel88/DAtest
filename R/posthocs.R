#' Run \code{drop1} on all features from \code{DAtest} results with \code{allResults = TRUE}
#'
#' Works on "zpo", "znb", "qpo", "neb", "poi". Non-paired "lrm", "llm", "llm2", "lma", "lmc"
#' @param results Output from a \code{DA."test"} function with \code{allResults = TRUE}
#' @param test Which test to use to calculate p-values. See \code{drop1} for details. Default "Chisq"
#' @param p.adj P-value adjustment method. See \code{p.adjust} for details. Default "fdr"
#' @param ... Additional arguments for \code{drop1} function
#' @export
DA.drop1 <- function(results, test = "Chisq", p.adj = "fdr", ...){
  
  # Check input
  if(is.data.frame(results) | !is.list(results)) stop("results should be the output from DA.zpo, DA.znb, DA.qpo, DA.neb, DA.poi, DA.lrm, DA.lma, DA.lmc, DA.llm or DA.llm2 with allResults=TRUE")
  
  # Class
  k <- 1
  while(class(results[[k]])[1] == "NULL"){
    k<- k+1
  } 
  xclass <- class(results[[k]])
  
  # Check class
  if(any("lme" %in% xclass)) stop("drop1 does not work on mixed-effect linear models. Use DA.anova")
  if(!any(c("lm","glm","zeroinfl","negbin","glmerMod") %in% xclass)) stop("results should be the output from DA.zpo, DA.znb, DA.qpo, DA.neb, DA.poi, DA.lrm, DA.llm, DA.lma, DA.lmc, or DA.llm2 with allResults=TRUE")
  
  # Run tests
  xres <- lapply(results, function(x) tryCatch(drop1(x, test = test, ...),error = function(e) NA))
  
  # Extract results
  if(all(xclass == "lm")){
    AIC <- lapply(xres, function(x) tryCatch(x[,4],error = function(e) NA))
    pv <-  lapply(xres, function(x) tryCatch(x[,5],error = function(e) NA))
    
    AIC <- do.call(rbind,AIC[lapply(AIC, length) > 1])
    pv <- do.call(rbind,pv[lapply(pv, length) > 1])
    pva <- apply(pv, 2, function(x) p.adjust(x, method=p.adj))
    
    colnames(AIC) <- paste0("AIC_",rownames(drop1(results[[k]], ...)))
    colnames(pv) <- paste0("pval_",rownames(drop1(results[[k]], ...)))
    colnames(pva) <- paste0("pval.adj_",rownames(drop1(results[[k]], ...)))
    
    res <- cbind(AIC,pv,pva)
    
  }
  
  if(xclass[1] == "glm"){
    
    if(is.na(results[[k]]$aic)){
      
      Dev <- lapply(xres, function(x) tryCatch(x[,2],error = function(e) NA))
      pv <-  lapply(xres, function(x) tryCatch(x[,4],error = function(e) NA))
      
      Dev <- do.call(rbind,Dev[lapply(Dev, length) > 1])
      pv <- do.call(rbind,pv[lapply(pv, length) > 1])
      pva <- apply(pv, 2, function(x) p.adjust(x, method=p.adj))
      
      colnames(Dev) <- paste0("Deviance_",rownames(drop1(results[[k]], ...)))
      colnames(pv) <- paste0("pval_",rownames(drop1(results[[k]], ...)))
      colnames(pva) <- paste0("pval.adj_",rownames(drop1(results[[k]], ...)))
      
      res <- cbind(Dev,pv,pva)
      
    } else {
      
      AIC <- lapply(xres, function(x) tryCatch(x[,3],error = function(e) NA))
      LRT <- lapply(xres, function(x) tryCatch(x[,4],error = function(e) NA))
      pv <-  lapply(xres, function(x) tryCatch(x[,5],error = function(e) NA))
      
      AIC <- do.call(rbind,AIC[lapply(AIC, length) > 1])
      LRT <- do.call(rbind,LRT[lapply(LRT, length) > 1])
      pv <- do.call(rbind,pv[lapply(pv, length) > 1])
      pva <- apply(pv, 2, function(x) p.adjust(x, method=p.adj))
      
      colnames(AIC) <- paste0("AIC_",rownames(drop1(results[[k]], ...)))
      colnames(LRT) <- paste0("LRT",rownames(drop1(results[[k]], ...)))
      colnames(pv) <- paste0("pval_",rownames(drop1(results[[k]], ...)))
      colnames(pva) <- paste0("pval.adj_",rownames(drop1(results[[k]], ...)))
      
      res <- cbind(AIC,LRT,pv,pva)
      
    }
    
  }
  
  if(all(xclass == "zeroinfl" | xclass == "glmerMod")){
    AIC <- lapply(xres, function(x) tryCatch(x[,2],error = function(e) NA))
    LRT <- lapply(xres, function(x) tryCatch(x[,3],error = function(e) NA))
    pv <-  lapply(xres, function(x) tryCatch(x[,4],error = function(e) NA))
    
    AIC <- do.call(rbind,AIC[lapply(AIC, length) > 1])
    LRT <- do.call(rbind,LRT[lapply(LRT, length) > 1])
    pv <- do.call(rbind,pv[lapply(pv, length) > 1])
    pva <- apply(pv, 2, function(x) p.adjust(x, method=p.adj))
    
    colnames(AIC) <- paste0("AIC_",rownames(drop1(results[[k]], ...)))
    colnames(LRT) <- paste0("LRT",rownames(drop1(results[[k]], ...)))
    colnames(pv) <- paste0("pval_",rownames(drop1(results[[k]], ...)))
    colnames(pva) <- paste0("pval.adj_",rownames(drop1(results[[k]], ...)))
    
    res <- cbind(AIC,LRT,pv,pva)
    
  }
  
  if(xclass[1] == "negbin"){
    AIC <- lapply(xres, function(x) tryCatch(x[,3],error = function(e) NA))
    LRT <- lapply(xres, function(x) tryCatch(x[,4],error = function(e) NA))
    pv <-  lapply(xres, function(x) tryCatch(x[,5],error = function(e) NA))
    
    AIC <- do.call(rbind,AIC[lapply(AIC, length) > 1])
    LRT <- do.call(rbind,LRT[lapply(LRT, length) > 1])
    pv <- do.call(rbind,pv[lapply(pv, length) > 1])
    pva <- apply(pv, 2, function(x) p.adjust(x, method=p.adj))
    
    colnames(AIC) <- paste0("AIC_",rownames(drop1(results[[k]], ...)))
    colnames(LRT) <- paste0("LRT",rownames(drop1(results[[k]], ...)))
    colnames(pv) <- paste0("pval_",rownames(drop1(results[[k]], ...)))
    colnames(pva) <- paste0("pval.adj_",rownames(drop1(results[[k]], ...)))
    
    res <- cbind(AIC,LRT,pv,pva)
    
  }
  
  res <- as.data.frame(res[,colSums(is.na(res)) != nrow(res)])
  return(res)
}


#' Run \code{anova} on all features from \code{DAtest} results with \code{allResults = TRUE}
#'
#' Works on "lrm", "llm", "llm2", "lma", "lmc". Non-paired "neb"
#' @param results Output from a \code{DA."test"} function with \code{allResults = TRUE}
#' @param p.adj P-value adjustment method. See \code{p.adjust for details}. Default "fdr"
#' @param ... Additional arguments for \code{anova} function
#' @export
DA.anova <- function(results, p.adj = "fdr", ...){
  
  # Check input
  if(is.data.frame(results) | !is.list(results)) stop("results should be the output from DA.lrm, DA.lma, DA.lmc, DA.llm, DA.llm2 or DA.neb with allResults=TRUE")
  
  # Class
  k <- 1
  while(class(results[[k]])[1] == "NULL"){
    k<- k+1
  } 
  xclass <- class(results[[k]])
  
  # Check class
  if(any("glmerMod" %in% xclass)) stop("anova does not work on mixed-effect negative binomial models. Use DA.drop1")
  if(!any(c("lm","nebgin","lme") %in% xclass)) stop("results should be the output from DA.lrm, DA.lma, DA.lmc, DA.llm, DA.llm2 or DA.neb with allResults=TRUE")
  
  # Run tests
  if(all(xclass == "lme")){
    
    pv <-  lapply(results, function(x) tryCatch(anova(x, ...)[,4],error = function(e) NA))
    
    pv <- do.call(rbind,pv[lapply(pv, length) > 1])
    pva <- apply(pv, 2, function(x) p.adjust(x, method=p.adj))
    
    colnames(pv) <- paste0("pval_",rownames(anova(results[[k]], ...)))
    colnames(pva) <- paste0("pval.adj_",rownames(anova(results[[k]], ...)))
    
    res <- cbind(pv,pva)
    
  }
  
  if(xclass[1] == "negbin" | xclass[1] == "lm"){
    
    pv <-  lapply(results, function(x) tryCatch(anova(x)[,5],error = function(e) NA))
    
    pv <- do.call(rbind,pv[lapply(pv, length) > 1])
    pva <- apply(pv, 2, function(x) p.adjust(x, method=p.adj))
    
    colnames(pv) <- paste0("pval_",rownames(anova(results[[k]])))
    colnames(pva) <- paste0("pval.adj_",rownames(anova(results[[k]])))
    
    res <- cbind(pv,pva)
    
  }
  
  res <- as.data.frame(res[,colSums(is.na(res)) != nrow(res)])
  return(res)
}


#' Run \code{TukeyHSD} on all features from \code{DAtest} results with \code{allResults = TRUE}
#'
#' Works on "aov", "lao", "lao2", "aoc", "aoa"
#' @param results Output from a \code{DA."test"} function with \code{allResults = TRUE}
#' @param variable Which variable to test. Default predictor. Alternatively, the name of a covar
#' @param p.adj P-value adjustment method. See \code{p.adjust} for details
#' @param ... Additional arguments for \code{TukeyHSD} function
#' @export
DA.TukeyHSD <- function(results, variable = "predictor", p.adj = "fdr", ...){
  
  # Check input
  if(is.data.frame(results) | !is.list(results)) stop("results should be the output from DA.aov, DA.aoa, DA.aoc, DA.lao or DA.lao2 with allResults=TRUE")
  
  # Class
  k <- 1
  while(class(results[[k]])[1] == "NULL"){
    k<- k+1
  } 
  xclass <- class(results[[k]])
  
  # Check class and results
  if(xclass[1] != "aov") stop("results should be the output from DA.aov, DA.aoa, DA.aoc, DA.lao or DA.lao2 with allResults=TRUE")
  if(!variable %in% attr(results[[k]]$terms,"term.labels")) stop(paste(variable,"not found in the models."))
  
  # Run test
  pv <-  lapply(results, function(x) tryCatch(as.data.frame(TukeyHSD(x, ...)[variable])[,4],error = function(e) NA))
  
  # Extract results
  pvs <- do.call(rbind,pv[lapply(pv, length) > 1])
  colnames(pvs) <- rownames(as.data.frame(TukeyHSD(results[[k]], ...)[variable]))
  pva <- apply(pvs, 2, function(x) p.adjust(x, method=p.adj))
  
  colnames(pvs) <- paste0("pval_",colnames(pvs))
  colnames(pva) <- paste0("pval.adj_",colnames(pva))
  
  res <- as.data.frame(cbind(pvs,pva))
  
  return(res)
}


#' Run \code{lsmeans} on all features from \code{DAtest} results with \code{allResults = TRUE}
#'
#' Pairwise testing on predictor and covars. Works on "poi", "neb", "lrm", "lma", "lmc", "llm", "llm2", "qpo", "znb", "zpo".
#' 
#' Require the \code{lsmeans} package
#' @param results Output from a \code{DA."test"} function with \code{allResults = TRUE}
#' @param variable Which variable to test. Default predictor. Alternatively, the name of a covar
#' @param predictor If results come from a paired "lrm", "llm", "lma", "lmc" or "llm2" supply the original predictor variable in the form of as a vector
#' @param covars If results come from a paired "lrm", "lma", "lmc", "llm" or "llm2" supply the original covars in the form of a named list
#' @param p.adj P-value adjustment method. See \code{p.adjust} for details
#' @param ... Additional arguments for \code{lsmeans} function
#' @export
DA.lsmeans <- function(results, variable = "predictor", predictor = NULL, covars = NULL, p.adj = "fdr", ...){

  # Check input
  if(is.data.frame(results) | !is.list(results)) stop("results should be the output from DA.poi, DA.neb, DA.lrm, DA.lma, DA.lmc, DA.llm, DA.llm2, DA.qpo, DA.znb or DA.zpo with allResults=TRUE")

  library(lsmeans)
  
  # Class
  k <- 1
  while(class(results[[k]])[1] == "NULL"){
    k<- k+1
  } 
  xclass <- class(results[[k]])
  
  # Check class and extract covars if necessary
  if(!any(c("lm","lme","glm","zeroinfl","negbin","glmerMod") %in% xclass)) stop("results should be the output from DA.zpo, DA.znb, DA.qpo, DA.neb, DA.poi, DA.lrm, DA.lma, DA.lmc, DA.llm or DA.llm2 with allResults=TRUE")
  
  if(class(results[[k]])[1] == "lme"){
    form <<- as.formula(paste("x ~",paste(attr(results[[1]]$terms,"term.labels"), collapse = "+")))
    if(is.null(predictor)) stop("predictor has to be supplied for a paired lrm, lma, lmc, llm and llm2")
    if(!is.null(covars)){
      for(i in 1:length(covars)){
        assign(names(covars)[i],covars[[i]])
      }
    }
  }

  # Run test and extract p-values and estimates
  mc <- lapply(results, function(x) tryCatch(summary(pairs(lsmeans(x, variable))),error = function(e) NA))
  pv <- lapply(mc, function(x) as.data.frame(x)$p.value)
  est <- lapply(mc, function(x) as.data.frame(x)$estimate)
  
  # Combine results
  pvs <- do.call(rbind,pv[lapply(pv, length) > 1])
  est <- do.call(rbind,est[lapply(est, length) > 1])
  colnames(pvs) <- summary(pairs(lsmeans(results[[k]], variable)))$contrast
  pva <- apply(pvs, 2, function(x) p.adjust(x, method=p.adj))
  
  colnames(est) <- paste0("estimate_",colnames(pvs))
  colnames(pvs) <- paste0("pval_",colnames(pvs))
  colnames(pva) <- paste0("pval.adj_",colnames(pva))
  
  res <- as.data.frame(cbind(est,pvs,pva))
  if(class(results[[k]])[1] == "lme") rm(form, envir = .GlobalEnv)
  return(res)
}
