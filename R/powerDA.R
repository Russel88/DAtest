#' Estimating (empirical) statistical power
#'
#' Estimating (empirical) statistical power for a specific differential abundance and expression method on a specific dataset
#' @param data Either a data.frame with counts/abundances, OR a \code{phyloseq} object. If a data.frame is provided rows should be taxa/genes/proteins and columns samples, and there should be rownames
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation. If the \code{predictor} is numeric it will be treated as such in the analyses
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation.
#' @param covars Either a named list with covariates, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param test Character. Which test to include. See \code{testDA} for details on the implemented tests. 
#' @param effectSizes Numeric. The effect sizes for the spike-ins. Default \code{c(2,4,8,16,32)}
#' @param alpha.p p-value threshold for false positive rates. Default 0.05
#' @param alpha.q q-value threshold for determining significance for \code{empirical power}. Default 0.1. This will change \code{fdr.output} for "sam". 
#' @param p.adj Character. Method for p-value adjustment. See \code{p.adjust} for details. Default "fdr"
#' @param R Integer. Number of times to run the tests. Default 5
#' @param relative Logical. TRUE (default) for compositional data, FALSE for absolute abundances or pre-normalized data.
#' @param k Vector of length 3. Number of Features to spike in each tertile (lower, mid, upper). E.g. \code{k=c(5,10,15)}: 5 features spiked in low abundance tertile, 10 features spiked in mid abundance tertile and 15 features spiked in high abundance tertile. Default NULL, which will spike 2 percent of the total amount of features in each tertile (a total of 6 percent), but minimum c(5,5,5)
#' @param cores Integer. Number of cores to use for parallel computing. Default one less than available. Set to 1 for sequential computing.
#' @param args List. A list with arguments passed to method.
#' @param out.all If TRUE linear models will output results and p-values from \code{anova}/\code{drop1}, ds2/ds2x will run LRT and not Wald test, erq and erq2 will produce one p-value for the predictor, and limma will run F-tests. If FALSE will output results for 2. level of the \code{predictor}. If NULL (default) set as TRUE for multi-class predictors and FALSE otherwise
#' @param core.check If TRUE will make an interactive check that the amount of cores specified are desired. Only if \code{cores>20}. This is to ensure that the function doesn't automatically overloads a server with workers.  
#' @param verbose If TRUE will print informative messages
#' @details Currently implemented methods: see \code{testDA}
#' @return An object of class \code{DAPower}, which contains a list with 1: A data.frame with results, 2: alpha.p value, 3: alpha.q value
#' @import snow doSNOW foreach utils
#' @importFrom parallel detectCores
#' @importFrom pROC roc
#' @examples 
#' # Creating random count_table and predictor
#' set.seed(5)
#' mat <- matrix(rnbinom(1000, size = 0.5, mu = 500), nrow = 50, ncol = 20)
#' rownames(mat) <- 1:50
#' pred <- c(rep("Control", 10), rep("Treatment", 10))
#' 
#' # Running powerDA on Wilcoxon test to test it with different effect sizes
#' # This example uses 1 core (cores = 1). 
#' # Remove the cores argument to get it as high (and thereby fast) as possible.
#' res <- powerDA(data = mat, predictor = pred, test = "wil", cores = 1)
#' summary(res)
#' 
#' \donttest{
#' # Include a paired variable for dependent/blocked samples
#' subject <- rep(1:10, 2)
#' res <- powerDA(data = mat, predictor = pred, paired = subject, test = "ttt")
#' 
#' # Include covariates
#' covar1 <- rnorm(20)
#' covar2 <- rep(c("A","B"), 10)
#' res <- powerDA(data = mat, predictor = pred, 
#'                covars = list(FirstCovar = covar1, CallItWhatYouWant = covar2), test = "lrm")
#' 
#' # Data is absolute abundance
#' res <- powerDA(data = mat, predictor = pred, relative = FALSE, test = "ttt")
#' }
#' @export

powerDA <- function(data, predictor, paired = NULL, covars = NULL, test = NULL, effectSizes = c(2,4,8,16,32), alpha.p = 0.05, alpha.q = 0.1, p.adj = "fdr", R = 5, relative = TRUE, k = NULL, cores = (detectCores()-1), args = list(), out.all = NULL, core.check = TRUE, verbose = TRUE){

  stopifnot(exists("data"),exists("predictor"))
  # Check for servers
  if(core.check){
    if(cores > 20){
      ANSWER <- readline(paste("You are about to run testDA using",cores,"cores. Enter y to proceed "))
      if(ANSWER != "y") stop("Process aborted")
    }
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
  
  # Coerce data
  if(!is.null(paired)) paired <- as.factor(paired)
  count_table <- as.matrix(count_table)
  
  # Checks
  if(is.null(test)) stop("'test' has to be specified")
  if(relative) if(!isTRUE(all(unlist(count_table) == floor(unlist(count_table))))) stop("count_table must only contain integer values")
  if(min(count_table) < 0) stop("count_table contains negative values")
  if(sum(colSums(count_table) == 0) > 0) stop("Some samples are empty")
  if(ncol(count_table) != length(predictor)) stop("Number of samples in count_table does not match length of predictor")
  if(length(unique(predictor)) < 2) stop("predictor should have at least two levels")
  if(length(test) != 1) stop("'test' has to have length 1")
  
  if(verbose) message(paste("Running on",cores,"cores"))

  # Remove Features not present in any samples
  if(sum(rowSums(count_table) == 0) != 0) message(paste(sum(rowSums(count_table) == 0),"empty features removed"))
  count_table <- count_table[rowSums(count_table) > 0,]
  if(nrow(count_table) <= 15) warning("Dataset contains very few features") 
  
  # Spike vs no features
  if(is.null(k)){
    k <- rep(round(nrow(count_table)*0.02),3)
    if(sum(k) < 15){
      k <- c(5,5,5)
    } 
  } 
  if(sum(k) == nrow(count_table)) stop("Set to spike all features. Change k argument")
  if(sum(k) > nrow(count_table)) stop("Set to spike more features than are present in the data. Change k argument")
  if(sum(k) < 15 & sum(k) >= 10 & R <= 10) message("Few features spiked. Increase 'k' or set 'R' to more than 10 to ensure proper estimations")
  if(sum(k) < 10 & sum(k) >= 5 & R <= 20) message("Few features spiked. Increase 'k' or set 'R' to more than 20 to ensure proper estimations")                                  
  if(sum(k) < 5 & R <= 50) message("Very few features spiked. Increase 'k' or set 'R' to more than 50 to ensure proper estimations")
  if(sum(k) > nrow(count_table)/2) message("Set to spike more than half of the dataset, which might give unreliable estimates, Change k argument")
                                   
  # predictor
  if(verbose) if(any(is.na(predictor))) warning("Predictor contains NAs!")
  if(is.numeric(predictor[1])){
    num.pred <- TRUE
    if(verbose) message(paste("predictor is assumed to be a quantitative variable, ranging from",min(predictor, na.rm = TRUE),"to",max(predictor, na.rm = TRUE)))
    if(length(levels(as.factor(predictor))) == 2){
      ANSWER <- readline("The predictor is quantitative, but only contains 2 unique values. Are you sure this is correct? Enter y to proceed ")
      if(ANSWER != "y") stop("Wrap the predictor with as.factor(predictor) to treat it is a categorical variable")
    }
  } else {
    num.pred <- FALSE
    if(length(levels(as.factor(predictor))) > length(unique(predictor))) stop("predictor has more levels than unique values!")
    if(verbose) message(paste("predictor is assumed to be a categorical variable with",length(unique(predictor)),"levels:",paste(levels(as.factor(predictor)),collapse = ", ")))
  }
  if(!is.null(paired)){
    if(verbose) message(paste("The paired variable has",length(unique(paired)),"levels"))
  }
  
  # out.all
  if(is.null(out.all)){
    if(length(unique(predictor)) == 2) out.all <- FALSE
    if(length(unique(predictor)) > 2) out.all <- TRUE
    if(num.pred) out.all <- FALSE
  }
  
  # Covars
  if(!is.null(covars)){
    for(i in seq_along(covars)){
      if(verbose) if(any(is.na(covars[[i]]))) warning(names(covars)[i],"contains NAs!")
      if(is.numeric(covars[[i]][1])){
        if(verbose) message(paste(names(covars)[i],"is assumed to be a quantitative variable, ranging from",min(covars[[i]], na.rm = TRUE),"to",max(covars[[i]], na.rm = TRUE)))
      } else {
        if(verbose) message(paste(names(covars)[i],"is assumed to be a categorical variable with",length(unique(covars[[i]])),"levels:",paste(levels(as.factor(covars[[i]])),collapse = ", ")))
      }
    }
  }
  
  if(verbose) cat("Spikeing...\n")
  # Shuffle predictor
  if(is.null(paired)){
    rands <- lapply(seq_len(R),function(x) sample(predictor))
  } else {
    rands <- lapply(seq_len(R),function(x) unsplit(lapply(split(predictor,paired), sample), paired))
  }
  
  # Spikeins
  spikeds.l <- list()
  for(eff in seq_along(effectSizes)){
    spikeds.l[[eff]] <- lapply(seq_len(R),function(x) spikein(count_table, rands[[x]], effectSizes[eff],  k, num.pred, relative))
  }
  spikeds <- do.call(c, spikeds.l)
  count_tables <- lapply(seq_len(R*length(effectSizes)),function(x) spikeds[[x]][[1]])
  
  # Test list
  tests.par <- paste0(unlist(lapply(effectSizes,function(x) rep(x,R))),"-",rep(paste0(seq_len(R),"_",rep(test,R)),R))

  ### Run tests
  # Progress bar
  pb <- txtProgressBar(max = length(tests.par), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Start parallel
  if(cores == 1) {
    registerDoSEQ() 
  } else {
    cl <- parallel::makeCluster(cores)
    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
  }

  # Run the tests in parallel
  results <- foreach(i = tests.par, .options.snow = opts) %dopar% {

    # Extract run info
    what.run <- which(i == tests.par)
    run.no <- as.numeric(gsub(".*-","",gsub("_.*","",i)))
    i <- gsub(".*_","",i)

    if(!is.na(pmatch("zzz",i))){
      i <- "zzz"
    } 

    # Run tests
    res.sub <- tryCatch(switch(i,
                               zzz = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars), args)),
                               mva = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,relative, p.adj), args)),
                               wil = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired, relative,p.adj), args)),
                               ttt = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired, relative,p.adj), args)),
                               ttr = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired, relative,p.adj), args)),
                               ltt = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,relative,p.adj), args)),
                               ttc = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,p.adj), args)),
                               tta = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,p.adj), args)),
                               ltt2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,p.adj), args)),
                               neb = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,relative,out.all,p.adj), args)),
                               erq = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj), args)),
                               ere = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],p.adj), args)),
                               erq2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj), args)),
                               ere2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],p.adj), args)),
                               msf = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],p.adj), args)),
                               zig = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,p.adj), args)),
                               ds2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj), args)),
                               ds2x = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj), args)),
                               per = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired, relative,p.adj), args)),
                               bay = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]]), args)),
                               adx = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]]), args)),
                               lim = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,relative,out.all,p.adj), args)),
                               lic = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj), args)),
                               lia = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj), args)),
                               lli = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,relative,out.all,p.adj), args)),
                               lli2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj), args)),
                               kru = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],relative,p.adj), args)),
                               aoa = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,p.adj), args)),
                               aoc = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,p.adj), args)),
                               aov = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,relative,p.adj), args)),
                               lao = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,relative,p.adj), args)),
                               lao2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,p.adj), args)),
                               lrm = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars, relative,out.all,p.adj), args)),
                               lmc = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj), args)),
                               lma = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj), args)),
                               llm = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,relative,out.all,p.adj), args)),
                               llm2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj), args)),
                               spe = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],relative,p.adj), args)),
                               pea = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],relative,p.adj), args)),
                               poi = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,relative,out.all,p.adj), args)),
                               qpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,relative,out.all,p.adj), args)),
                               vli = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,covars,out.all,p.adj), args)),
                               zpo = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,relative,out.all,p.adj), args)),
                               znb = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],covars,relative,out.all,p.adj), args)),
                               fri = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,relative,p.adj), args)),
                               qua = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,relative,p.adj), args)),
                               sam = do.call(get(noquote(paste0("DA.",i))),c(list(count_tables[[what.run]],rands[[run.no]],paired,fdr.output = alpha.q), args))),
                        error = function(e) NULL)
    
    if(!i %in% c("sam","adx")){
      res.sub[is.na(res.sub$pval),"pval"] <- 1
      res.sub[is.na(res.sub$pval.adj),"pval.adj"] <- 1
    }
    
    return(res.sub)
    
  }
  names(results) <- tests.par
  
  # Split aldex
  if(test == "adx"){
    results.t <- lapply(results, function(x) x[,c(1:9,12,13)])
    results.w <- lapply(results, function(x) x[,c(1:7,10:13)])
    for(k in seq_along(results.t)){
      colnames(results.t[[k]]) <- c(colnames(results[[1]])[1:7],"pval","pval.adj","ordering","Feature")
      colnames(results.w[[k]]) <- c(colnames(results[[1]])[1:7],"pval","pval.adj","ordering","Feature")
      results.t[[k]]$Method <- "ALDEx2 t-test (adx)"
      results.w[[k]]$Method <- "ALDEx2 wilcox (adx)"
    }
    names(results.t) <- paste0(names(results.t),".t")
    names(results.w) <- paste0(names(results.w),".w")
    results <- c(results.t,results.w)
    tests.par <- names(results)
    spikeds <- c(spikeds,spikeds)
  }
  
  r <- NULL
  final.results <- foreach(r = tests.par, .combine = rbind) %do% {

    res.sub <- results[names(results) == r][[1]]
    
    # Make pseudo-pval for SAMseq
    if(test == "sam"){
      res.sub$pval <- 1/rank(res.sub$Score)
      res.sub$pval.adj <- 1
      res.sub[res.sub$Sig == "Yes","pval.adj"] <- 0 
    }
    
    # Confusion matrix
    totalPos <- nrow(res.sub[res.sub$pval <= alpha.p,])
    totalNeg <- nrow(res.sub[res.sub$pval > alpha.p,])
    truePos <- sum(res.sub[res.sub$pval <= alpha.p,"Feature"] %in% spikeds[[which(r == tests.par)]][[2]])
    falseNeg <- sum(res.sub[res.sub$pval > alpha.p,"Feature"] %in% spikeds[[which(r == tests.par)]][[2]])
    falsePos <- totalPos - truePos
    trueNeg <- totalNeg - falseNeg
    
    # FPR
    if(test == "sam"){
      fpr <- NA
    } else {
      if((falsePos + trueNeg) != 0){
        fpr <- falsePos / (falsePos + trueNeg)
      } else {
        fpr <- 0
      }
    }

    # Confusion matrix adjusted
    totalPos.adj <- nrow(res.sub[res.sub$pval.adj <= alpha.q,])
    truePos.adj <- sum(res.sub[res.sub$pval.adj <= alpha.q,"Feature"] %in% spikeds[[which(r == tests.par)]][[2]])
    falsePos.adj <- totalPos.adj - truePos.adj
    
    # FDR 
    if(totalPos.adj != 0){
      fdr <- falsePos.adj / totalPos.adj
    } else {
      fdr <- 0
    }
    
    # Spike detection rate (empircal power aka sensitivity)
    sdr <- truePos.adj / sum(k)
    
    # AUC
    test_roc <- NULL
     tryCatch(
       test_roc <- pROC::roc(as.numeric(res.sub$Feature %in% spikeds[[which(r == tests.par)]][[2]]) ~ res.sub$pval, auc=TRUE, direction = ">", quiet=TRUE),
       error = function(e) NULL)
     if(!is.null(test_roc)){
       auc <- as.numeric(test_roc$auc) 
     } else {
       auc <- 0.5
     }

    # Combine and return
    df.combined <- data.frame(Method = res.sub$Method[1],
                              Run = as.numeric(gsub(".*-","",gsub("_.*","",r))),
                              EffectSize = as.numeric(gsub("-.*","",gsub("_.*","",r))),
                              AUC = auc,
                              FPR = fpr,
                              FDR = fdr,
                              Power = sdr)
    rownames(df.combined) <- NULL

    return(df.combined)
    
  }
  final <- list(final.results,alpha.p,alpha.q)
  class(final) <- "DAPower"
  return(final)
}
