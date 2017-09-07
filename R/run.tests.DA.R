#' @export

run.tests.DA <- function(count_table, outcome, paired, tests, relative, args, cores){

  # Progress bar
  pb <- txtProgressBar(max = length(tests), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # Start parallel
  if(cores == 1) {
    registerDoSEQ() 
  } else {
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
  }

  results <- foreach(i = tests , .options.snow = opts) %dopar% {
  
    # Extract test arguments
    if(!all(names(args) %in% tests)) stop("One or more names in list with additional arguments does not match names of tests")
    for(j in seq_along(args)){
      assign(paste0(names(args)[j],".args"),args[[j]],pos=1)
    }
    test.args <- paste0(tests,".args")
    test.boo <- lapply(test.args,exists)
    for(l in seq_along(test.args)){
      if(test.boo[l] == FALSE) assign(test.args[l], list(),pos=1)
    }
  
    res.sub <- switch(i,
                    wil = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, relative),wil.args)),
                    ttt = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, relative),ttt.args)),
                    ltt = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired,relative),ltt.args)),
                    ltt2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired),ltt2.args)),
                    neb = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired),neb.args)),
                    erq = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired),erq.args)),
                    ere = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome),ere.args)),
                    msf = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome),msf.args)),
                    zig = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome),zig.args)),
                    ds2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired),ds2.args)),
                    per = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, relative),per.args)),
                    bay = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired),bay.args)),
                    adx = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome),adx.args)),
                    lim = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired,relative),lim.args)),
                    lli = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired,relative),lli.args)),
                    lli2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired),lli2.args)),
                    kru = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome, relative),kru.args)),
                    aov = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome, relative),aov.args)),
                    lao = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,relative),lao.args)),
                    lao2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome),lao2.args)),
                    lrm = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, relative),lrm.args)),
                    llm = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired,relative),llm.args)),
                    llm2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired),llm2.args)),
                    rai = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome),rai.args)),
                    spe = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,relative),spe.args)))
  
    res.sub[is.na(res.sub$pval),"pval"] <- 1

    return(res.sub)
  
  }
  names(results) <- tests

  # Split ALDEx2 results in t.test and wilcoxon
  if("adx" %in% tests){
    adx.t <- as.data.frame(results["adx"])[,c(1:7,12)]
    adx.w <- as.data.frame(results["adx"])[,c(1:7,12)]
    colnames(adx.t) <- gsub("adx.","",colnames(adx.t))
    colnames(adx.w) <- colnames(adx.t)
    adx.t$pval <- as.numeric(as.data.frame(results["adx"])$adx.we.ep)
    adx.w$pval <- as.numeric(as.data.frame(results["adx"])$adx.wi.ep)
    adx.t$Method <- "ALDEx2 t-test"
    adx.w$Method <- "ALDEx2 wilcox"
    results["adx"] <- NULL
    res.names <- names(results)
    results <- c(results,list(adx.t),list(adx.w))
    names(results) <- c(res.names,"adx.t","adx.w")
  }

  return(results)

}