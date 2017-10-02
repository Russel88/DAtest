#' Run drop1 on all features from DAtest results with allResults = TRUE
#'
#' Works on "zpo", "znb", "qpo", "neb", "poi". Non-paired "lrm", "llm", "llm2"
#' @param results Output from a DA."test" function with allResults = TRUE
#' @param test Which test to use to calculate p-values. See ?drop1 for details
#' @param p.adj P-value adjustment method. See ?p.adjust for details
#' @param ... Additional arguments for drop1 function
#' @export
DA.drop1 <- function(results, test = "Chisq", p.adj = "fdr", ...){
  
  # Class
  k <- 1
  while(class(results[[k]])[1] == "NULL"){
    k<- k+1
  } 
  xclass <- class(results[[k]])
  
  if(!any(c("lm","glm","zeroinfl","negbin","glmerMod") %in% xclass)) stop(paste("Class should be one of lm, glm, zeroinfl, negbin or glmerMod and not:",xclass))
  
  xres <- lapply(results, function(x) tryCatch(drop1(x, test = test, ...),error = function(e) NA))
  
  if(all(xclass == "lm")){
    AIC <- lapply(xres, function(x) x[,4])
    pv <-  lapply(xres, function(x) x[,5])
    
    AIC <- do.call(rbind,AIC[lapply(AIC, length) > 1])
    pv <- do.call(rbind,pv[lapply(pv, length) > 1])
    pva <- apply(pv, 2, function(x) p.adjust(x, method=p.adj))
    
    colnames(AIC) <- paste0("AIC_",rownames(drop1(results[[1]], ...)))
    colnames(pv) <- paste0("pval_",rownames(drop1(results[[1]], ...)))
    colnames(pva) <- paste0("pval.adj_",rownames(drop1(results[[1]], ...)))
    
    res <- cbind(AIC,pv,pva)
    
  }
  
  if(xclass[1] == "glm"){
    
    if(is.na(results[[k]]$aic)){
      
      Dev <- lapply(xres, function(x) x[,2])
      pv <-  lapply(xres, function(x) x[,4])
      
      Dev <- do.call(rbind,Dev[lapply(Dev, length) > 1])
      pv <- do.call(rbind,pv[lapply(pv, length) > 1])
      pva <- apply(pv, 2, function(x) p.adjust(x, method=p.adj))
      
      colnames(Dev) <- paste0("Deviance_",rownames(drop1(results[[1]], ...)))
      colnames(pv) <- paste0("pval_",rownames(drop1(results[[1]], ...)))
      colnames(pva) <- paste0("pval.adj_",rownames(drop1(results[[1]], ...)))
      
      res <- cbind(Dev,pv,pva)
      
    } else {
      
      AIC <- lapply(xres, function(x) x[,3])
      LRT <- lapply(xres, function(x) x[,4])
      pv <-  lapply(xres, function(x) x[,5])
      
      AIC <- do.call(rbind,AIC[lapply(AIC, length) > 1])
      LRT <- do.call(rbind,LRT[lapply(LRT, length) > 1])
      pv <- do.call(rbind,pv[lapply(pv, length) > 1])
      pva <- apply(pv, 2, function(x) p.adjust(x, method=p.adj))
      
      colnames(AIC) <- paste0("AIC_",rownames(drop1(results[[1]], ...)))
      colnames(LRT) <- paste0("LRT",rownames(drop1(results[[1]], ...)))
      colnames(pv) <- paste0("pval_",rownames(drop1(results[[1]], ...)))
      colnames(pva) <- paste0("pval.adj_",rownames(drop1(results[[1]], ...)))
      
      res <- cbind(AIC,LRT,pv,pva)
      
    }
    
  }
  
  if(all(xclass == "zeroinfl" | xclass == "glmerMod")){
    AIC <- lapply(xres, function(x) x[,2])
    LRT <- lapply(xres, function(x) x[,3])
    pv <-  lapply(xres, function(x) x[,4])
    
    AIC <- do.call(rbind,AIC[lapply(AIC, length) > 1])
    LRT <- do.call(rbind,LRT[lapply(LRT, length) > 1])
    pv <- do.call(rbind,pv[lapply(pv, length) > 1])
    pva <- apply(pv, 2, function(x) p.adjust(x, method=p.adj))
    
    colnames(AIC) <- paste0("AIC_",rownames(drop1(results[[1]], ...)))
    colnames(LRT) <- paste0("LRT",rownames(drop1(results[[1]], ...)))
    colnames(pv) <- paste0("pval_",rownames(drop1(results[[1]], ...)))
    colnames(pva) <- paste0("pval.adj_",rownames(drop1(results[[1]], ...)))
    
    res <- cbind(AIC,LRT,pv,pva)
    
  }
  
  if(xclass[1] == "negbin"){
    AIC <- lapply(xres, function(x) x[,3])
    LRT <- lapply(xres, function(x) x[,4])
    pv <-  lapply(xres, function(x) x[,5])
    
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
  
  return(res)
}


#' Run anova on all features from DAtest results with allResults = TRUE
#'
#' Works on "lrm", "llm", "llm2". Non-paired "neb"
#' @param results Output from a DA."test" function with allResults = TRUE
#' @param p.adj P-value adjustment method. See ?p.adjust for details
#' @param ... Additional arguments for anova function
#' @export
DA.anova <- function(results, p.adj = "fdr", ...){
  
  # Class
  k <- 1
  while(class(results[[k]])[1] == "NULL"){
    k<- k+1
  } 
  xclass <- class(results[[k]])
  
  if(!any(c("lm","nebgin","lme") %in% xclass)) stop(paste("Class should be one of lm, lme or negbin and not:",xclass))
  
  if(all(xclass == "lme")){
    
    pv <-  lapply(results, function(x) tryCatch(anova(x, ...)[,4],error = function(e) NA))
    
    pv <- do.call(rbind,pv[lapply(pv, length) > 1])
    pva <- apply(pv, 2, function(x) p.adjust(x, method=p.adj))
    
    colnames(pv) <- paste0("pval_",rownames(anova(results[[1]], ...)))
    colnames(pva) <- paste0("pval.adj_",rownames(anova(results[[1]], ...)))
    
    res <- cbind(pv,pva)
    
  }
  
  if(xclass[1] == "negbin" | xclass[1] == "lm"){
    
    pv <-  lapply(results, function(x) tryCatch(anova(x)[,5],error = function(e) NA))
    
    pv <- do.call(rbind,pv[lapply(pv, length) > 1])
    pva <- apply(pv, 2, function(x) p.adjust(x, method=p.adj))
    
    colnames(pv) <- paste0("pval_",rownames(anova(results[[1]])))
    colnames(pva) <- paste0("pval.adj_",rownames(anova(results[[1]])))
    
    res <- cbind(pv,pva)
    
  }
  
  
  return(res)
}


#' Run TukeyHSD on all features from DAtest results with allResults = TRUE
#'
#' Works on "aov", "lao", "lao2"
#' @param results Output from a DA."test" function with allResults = TRUE
#' @param variable Which variable to test. Default predictor. Alternatively, the name of a covar
#' @param p.adj P-value adjustment method. See ?p.adjust for details
#' @param ... Additional arguments for TukeyHSD function
#' @export
DA.TukeyHSD <- function(results, variable = "predictor", p.adj = "fdr", ...){
  
  # Class
  k <- 1
  while(class(results[[k]])[1] == "NULL"){
    k<- k+1
  } 
  xclass <- class(results[[k]])
  
  if(xclass[1] != "aov") stop(paste("Class should be aov and not:",xclass))
  if(!variable %in% attr(results[[k]]$terms,"term.labels")) stop(paste(variable,"not found in the models."))
  
  pv <-  lapply(results, function(x) tryCatch(as.data.frame(TukeyHSD(x, ...)[variable])[,4],error = function(e) NA))
  
  pvs <- do.call(rbind,pv[lapply(pv, length) > 1])
  colnames(pvs) <- rownames(as.data.frame(TukeyHSD(results[[k]], ...)[variable]))
  pva <- apply(pvs, 2, function(x) p.adjust(x, method=p.adj))
  
  colnames(pvs) <- paste0("pval_",colnames(pvs))
  colnames(pva) <- paste0("pval.adj_",colnames(pva))
  
  res <- cbind(pvs,pva)
  
  return(res)
}


#' Run lsmeans on all features from DAtest results with allResults = TRUE
#'
#' Pairwise testing on predictor and covars. Works on "poi", "neb", "lrm", "llm", "llm2", "qpo", "znb", "zpo".
#' 
#' Require the lsmeans package
#' @param results Output from a DA."test" function with allResults = TRUE
#' @param variable Which variable to test. Default predictor. Alternatively, the name of a covar
#' @param predictor If results come from a paired "lrm", "llm" or "llm2" supply the original predictor variable in the form of as a vector
#' @param covars If results come from a paired "lrm", "llm" or "llm2" supply the original covars in the form of a named list
#' @param p.adj P-value adjustment method. See ?p.adjust for details
#' @param ... Additional arguments for lsmeans function
#' @export
DA.lsmeans <- function(results, variable = "predictor", predictor = NULL, covars = NULL, p.adj = "fdr", ...){
  
  library(lsmeans)
  
  # Class
  k <- 1
  while(class(results[[k]])[1] == "NULL"){
    k<- k+1
  } 
  
  if(class(results[[k]])[1] == "lme"){
    form <<- as.formula(paste("x ~",paste(attr(results[[1]]$terms,"term.labels"), collapse = "+")))
    if(is.null(predictor)) stop("predictor has to be supplied for a paired lrm, llm and llm2")
    if(!is.null(covars)){
      for(i in 1:length(covars)){
        assign(names(covars)[i],covars[[i]])
      }
    }
  }
  
  mc <- lapply(results, function(x) tryCatch(summary(pairs(lsmeans(x, variable))),error = function(e) NA))
  pv <- lapply(mc, function(x) as.data.frame(x)$p.value)
  est <- lapply(mc, function(x) as.data.frame(x)$estimate)
  
  pvs <- do.call(rbind,pv[lapply(pv, length) > 1])
  est <- do.call(rbind,est[lapply(est, length) > 1])
  colnames(pvs) <- summary(pairs(lsmeans(results[[k]], variable)))$contrast
  pva <- apply(pvs, 2, function(x) p.adjust(x, method=p.adj))
  
  colnames(est) <- paste0("estimate_",colnames(pvs))
  colnames(pvs) <- paste0("pval_",colnames(pvs))
  colnames(pva) <- paste0("pval.adj_",colnames(pva))
  
  res <- cbind(est,pvs,pva)
  if(class(results[[k]])[1] == "lme") rm(form, envir = .GlobalEnv)
  return(res)
}

