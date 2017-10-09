#' Estimate runtime of testDA on large datasets
#' 
#' Estimate the runtime of testDA from running on a subset of the features. Intended for datasets with at least 5000 features.
#' 
#' Runtime of all methods are expected to scale linearly with the number of features, except "anc" and "bay" which are modelled with a 2. order polynomial.
#' @param data Either a matrix with counts/abundances, OR a phyloseq object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples, and there should be rownames
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if data is a phyloseq object the name of the variable in sample_data in quotation. If the predictor is numeric it will be treated as such in the analyses
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if data is a phyloseq object the name of the variable in sample_data in quotation. Only for "anc", "poi", "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq", "ds2", "lrm", "llm", "llm2", "lim", "lli", "lli2" and "zig"
#' @param covars Either a named list with covariates, OR if data is a phyloseq object a character vector with names of the variables in sample_data(data)
#' @param subsamples Vector with numbers of features to subsample to estimate runtime for fast methods
#' @param subsamples.slow Vector with numbers of features to subsample to estimate runtime for slow methods
#' @param tests Fast methods to include
#' @param tests.slow Slow methods to include
#' @param R Intended number of repeats for the testDA function
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available. Set to 1 for sequential computing.
#' @param print.res If TRUE will print the results, alternatively will return a data.frame with the results.
#' @param ... Additional arguments for the testDA function
#' @return A data.frame if print.res is FALSE
#' @importFrom parallel detectCores
#' @export
runtimeDA <- function(data, predictor, paired = NULL, covars = NULL, subsamples = c(500,1000,1500,2000,2500), subsamples.slow = c(100,200,300,400,500), 
                      tests =  c("ds2","sam", "qua", "fri", "vli", "qpo", "poi", "pea", "wil", "ttt", "ltt", "ltt2", "erq", "erq2","ere", "ere2", "msf", "zig", "lim", "lli", "lli2", "aov", "lao", "lao2", "kru", "lrm", "llm", "llm2", "spe"), 
                      tests.slow = c("neb", "bay", "per", "zpo", "znb", "rai", "adx"), R = 10, cores = (detectCores()-1), print.res = TRUE, ...){
  
  stopifnot(exists("data"),exists("predictor"))

  if(cores > 10){
    ANSWER <- readline(paste("You are about to run runtimeDA using",cores,"cores. Enter y to proceed "))
    if(ANSWER != "y") stop("Process aborted")
  }
  
  # Extract from phyloseq
  if(class(data) == "phyloseq"){
    if(length(predictor) > 1 | length(paired) > 1) stop("When data is a phyloseq object predictor and paired should only contain the name of the variables in sample_data")
    if(!predictor %in% sample_variables(data)) stop(paste(predictor,"is not present in sample_data(data)"))
    if(!is.null(paired)){
      if(!paired %in% sample_variables(data)) stop(paste(paired,"is not present in sample_data(data)"))
    }
    count_table <- otu_table(data)
    if(!taxa_are_rows(data)) count_table <- t(count_table)
    predictor <- suppressWarnings(as.matrix(sample_data(data)[,predictor]))
    if(!is.null(paired)) paired <- suppressWarnings(as.factor(as.matrix(sample_data(data)[,paired])))
    if(!is.null(covars)){
      covars.n <- covars
      covars <- list()
      for(i in 1:length(covars.n)){
        covars[[i]] <- suppressWarnings(as.matrix(sample_data(data)[,covars.n[i]]))
      }
      names(covars) <- covars.n
    } 
  } else {
    count_table <- data
  }
  
  # Remove Features not present in any samples
  if(sum(rowSums(count_table) == 0) != 0) message(paste(sum(rowSums(count_table) == 0),"empty features removed"))
  count_table <- count_table[rowSums(count_table) > 0,]
  
  # Run subsets fast
  message("Running fast methods")
  test.list <- list()
  for(i in 1:length(subsamples)){
    
    # Subset and ensure that no samples are empty
    j <- 0
    while(j == 0){
      sub <- count_table[sample(rownames(count_table),subsamples[i]),]
      if(any(colSums(sub) == 0)) j <- 0 else j <- 1
    }
    
    # Run test
    sub.test <- testDA(sub, predictor, paired, covars, R = 1, tests = tests, cores = cores, core.check = FALSE, ...)
    test.list[[i]] <- sub.test
  }

  # Run subsets slow
  message("Running slow methods")
  test.slow.list <- list()
  for(i in 1:length(subsamples.slow)){
    
    # Subset and ensure that no samples are empty
    j <- 0
    while(j == 0){
      sub <- count_table[sample(rownames(count_table),subsamples.slow[i]),]
      if(any(colSums(sub) == 0)) j <- 0 else j <- 1
    }
    
    # Run test
    sub.test <- testDA(sub, predictor, paired, covars, R = 1, tests = tests.slow, cores = cores, core.check = FALSE, ...)
    test.slow.list[[i]] <- sub.test
  }
  
  tests.list <- append(test.list,test.slow.list) 
  subsamps <- c(subsamples,subsamples.slow)
  
  # Extrapolate
  runtimes <- lapply(tests.list, function(x) x$run.times)
  runtimes <- lapply(1:length(runtimes), function(x) cbind(runtimes[[x]],subsamps[[x]]))
  runtimes <- do.call(rbind,runtimes)
  
  # Which tests have been run
  all.tests <- unique(rownames(runtimes))
  all.tests <- names(table(rownames(runtimes))[table(rownames(runtimes))>1])
  
  # Collect data
  extra <- data.frame(Test = all.tests,
                      Minutes = NA,
                      Minutes. = NA)
  colnames(extra)[3] <- paste0(colnames(extra[3]),"R=",R)
  
  for(i in all.tests){
    extra.sub <- runtimes[rownames(runtimes) == i,]
    if(i %in% c("anc","bay")){
      fit <- lm(Minutes ~ poly(V2,2), data = as.data.frame(extra.sub))
    } else {
      fit <- lm(Minutes ~ V2, data = as.data.frame(extra.sub))
    }
    extra[extra$Test == i,2] <- round(predict(fit, newdata = data.frame(V2 = nrow(count_table))),2)
    extra[extra$Test == i,3] <- round(predict(fit, newdata = data.frame(V2 = nrow(count_table)))*R,2)
  }
  
  extra[extra[,2] < 0,"Minutes"] <- 0
  extra[extra[,3] < 0,"Minutes"] <- 0
  
  # Order extra
  extra <- extra[order(extra$Minutes, decreasing = TRUE),]
  
  if(print.res){
    # Print the results
    message("Estimated run times.\nWith cores=1 the runtime will be the sum of them all.\nWith more cores the actual runtime will decrease asymptotically towards the slowest test")
    print(extra, row.names = FALSE)
  } else {
    return(extra)
  }
}