#' Estimate runtime of \code{testDA} on large datasets
#' 
#' Estimate the runtime of \code{testDA} from running on a subset of the features. Intended for datasets with at least 2000 features.
#' 
#' Outputs the estimated times for running each method 1 time. With cores=1 the runtime will be the sum of them all. With more cores the actual runtime will decrease asymptotically towards the slowest test
#' 
#' Runtime of all methods are expected to scale linearly with the number of features, except "anc" and "bay" which are modelled with a 2. order polynomial.
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples, and there should be rownames
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation. If the \code{predictor} is numeric it will be treated as such in the analyses
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation. Only for "poi", "per", "ttt", "ltt", "ltt2", "neb", "wil", "erq", "ds2", "lrm", "llm", "llm2", "lim", "lli", "lli2" and "zig"
#' @param covars Either a named list with covariates, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param subsamples Vector with numbers of features to subsample to estimate runtime for fast methods
#' @param subsamples.slow Vector with numbers of features to subsample to estimate runtime for slow methods
#' @param tests Fast methods to include
#' @param tests.slow Slow methods to include
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available. Set to 1 for sequential computing.
#' @param ... Additional arguments for the \code{testDA} function
#' @return A data.frame with estimated runtimes for 1 run
#' @importFrom parallel detectCores
#' @examples 
#' # Creating large random count_table and predictor
#' set.seed(5)
#' mat <- matrix(rnbinom(150000, size = 0.5, mu = 500), nrow = 10000, ncol = 10)
#' rownames(mat) <- 1:10000
#' pred <- c(rep("A", 5), rep("B", 5))
#' 
#' # Use runtimeDA to predict total runtime for all features
#' # This example uses 1 core (cores = 1). 
#' # Remove the cores argument to get it as high (and thereby fast) as possible.
#' # Also, in this example only a subset of tests are run.
#' runtimeDA(mat, pred, cores = 1, tests = c("ttt","wil"), tests.slow = c("neb"))
#' @export
runtimeDA <- function(data, predictor, paired = NULL, covars = NULL, subsamples = c(500,1000,1500,2000), subsamples.slow = c(100,150,200,250), 
                      tests =  c("sam", "qua", "fri", "vli", "qpo", "pea", "wil", "ttt", "ltt", "ltt2","ere", "ere2", "msf", "zig", "lim", "lli", "lli2", "aov", "lao", "lao2", "kru", "lrm", "llm", "llm2", "spe", "aoa", "aoc", "tta", "ttc", "lma", "lmc", "lia", "lic"), 
                      tests.slow = c("mva", "neb", "bay", "per", "ds2", "ds2x", "zpo", "znb", "adx", "poi", "erq", "erq2"), cores = (detectCores()-1), ...){
  
  stopifnot(exists("data"),exists("predictor"))

  if(cores > 20){
    ANSWER <- readline(paste("You are about to run runtimeDA using",cores,"cores. Enter y to proceed "))
    if(ANSWER != "y") stop("Process aborted")
  }
  
  # Extract from phyloseq
  if(is(data, "phyloseq")){
    DAdata <- DA.phyloseq(data, predictor, paired, covars)
    count_table <- DAdata$count_table
    predictor <- DAdata$predictor
    paired <- DAdata$paired
    covars <- DAdata$covars
  } else {
    count_table <- data
  }
  if(!is.null(covars)){
    for(i in seq_along(covars)){
      assign(names(covars)[i], covars[[i]])
    }
  }
  
  # Remove Features not present in any samples
  if(sum(rowSums(count_table) == 0) != 0) message(paste(sum(rowSums(count_table) == 0),"empty features removed"))
  count_table <- count_table[rowSums(count_table) > 0,]
  
  # Trim for small datasets
  subsamples <- subsamples[subsamples <= nrow(count_table)]
  subsamples.slow <- subsamples.slow[subsamples.slow <= nrow(count_table)]
  if(length(subsamples) < 2) stop("At least two subsamples are needed to estimate runtime")
  if(length(subsamples.slow) < 2) stop("At least two subsamples are needed to estimate runtime")
  
  # predictor
  if(is.numeric(predictor[1])){
    message("predictor is assumed to be a quantitative variable")
  } else {
    message(paste("predictor is assumed to be a categorical variable with",length(unique(predictor)),"levels:",paste(levels(as.factor(predictor)),collapse = ", ")))
  }
  
  # Covars
  if(!is.null(covars)){
    for(i in seq_along(covars)){
      if(is.numeric(covars[[i]][1])){
        message(paste(names(covars)[i],"is assumed to be a quantitative variable"))
      } else {
        message(paste(names(covars)[i],"is assumed to be a categorical variable with",length(unique(covars[[i]])),"levels:",paste(levels(as.factor(covars[[i]])),collapse = ", ")))
      }
    }
  }
  
  # Run subsets fast
  count_table <- count_table[order(rowSums(count_table), decreasing = TRUE), ]
  
  message("Running fast methods")
  test.list <- list()
  for(i in seq_along(subsamples)){
    cat(paste("\n",subsamples[i],"features"),fill = TRUE)
    # Subset
    sub <- count_table[seq_len(subsamples[i]), ]
    if(sum(colSums(sub) == 0) != 0) warning(paste(sum(colSums(sub) == 0),"empty samples removed"))
    sub <- sub[, colSums(sub) > 0]
    
    # Run test
    sub.test <- testDA(sub, predictor, paired, covars, R = 1, tests = tests, cores = cores, core.check = FALSE, verbose = FALSE, ...)
    test.list[[i]] <- sub.test
  }

  # Run subsets slow
  message("Running slow methods")
  test.slow.list <- list()
  for(i in seq_along(subsamples.slow)){
    cat(paste("\n",subsamples.slow[i],"features"),fill = TRUE)
    # Subset
    sub <- count_table[seq_len(subsamples.slow[i]),]
    if(sum(colSums(sub) == 0) != 0) warning(paste(sum(colSums(sub) == 0),"empty samples removed"))
    sub <- sub[, colSums(sub) > 0]
    
    # Run test
    sub.test <- testDA(sub, predictor, paired, covars, R = 1, tests = tests.slow, cores = cores, core.check = FALSE, verbose = FALSE)
    test.slow.list[[i]] <- sub.test
  }
  
  tests.list <- append(test.list,test.slow.list) 
  subsamps <- c(subsamples,subsamples.slow)
  
  # Extrapolate
  runtimes <- lapply(tests.list, function(x) x$run.times)
  runtimes <- lapply(seq_along(runtimes), function(x) cbind(runtimes[[x]],subsamps[[x]],rownames(runtimes[[x]])))
  runtimes <- do.call(rbind,runtimes)
  colnames(runtimes) <- c("Minutes","SubSamp","Test")
  
  # Which tests have been run
  all.tests <- names(table(runtimes$Test)[table(runtimes$Test)>1])
  
  # Collect data
  extra <- data.frame(Test = all.tests,
                      Minutes = NA)
  
  for(i in all.tests){
    extra.sub <- runtimes[runtimes$Test == i,]
    if(i %in% c("bay")){
      fit <- lm(Minutes ~ poly(SubSamp,2), data = as.data.frame(extra.sub))
    } else {
      fit <- lm(Minutes ~ SubSamp, data = as.data.frame(extra.sub))
    }
    extra[extra$Test == i,2] <- round(predict(fit, newdata = data.frame(SubSamp = nrow(count_table))),2)
  }
  extra[extra[,2] < 0,"Minutes"] <- 0

  # Order extra
  extra <- extra[order(extra$Minutes, decreasing = TRUE),]
  
  return(extra)

}
