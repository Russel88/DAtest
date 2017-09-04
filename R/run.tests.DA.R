#' @export

run.tests.DA <- function(count_table, outcome, paired, tests, relative, p.adj, args, cores){

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
                    wil = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, relative, p.adj),wil.args)),
                    ttt = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, relative, p.adj),ttt.args)),
                    ltt = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired,relative, p.adj),ltt.args)),
                    ltt2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, p.adj),ltt2.args)),
                    neb = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, p.adj),neb.args)),
                    erq = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, p.adj),erq.args)),
                    ere = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome, p.adj),ere.args)),
                    msf = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome, p.adj),msf.args)),
                    zig = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome, p.adj),zig.args)),
                    ds2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, p.adj),ds2.args)),
                    per = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, relative, p.adj),per.args)),
                    bay = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, p.adj),bay.args)),
                    adx = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome),adx.args)),
                    lim = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired,relative,p.adj),lim.args)),
                    lli = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired,relative,p.adj),lli.args)),
                    lli2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired,p.adj),lli2.args)),
                    kru = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome, relative, p.adj),kru.args)),
                    aov = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome, relative, p.adj),aov.args)),
                    lao = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,relative, p.adj),lao.args)),
                    lao2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome, p.adj),lao2.args)),
                    lrm = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, relative, p.adj),lrm.args)),
                    llm = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired,relative, p.adj),llm.args)),
                    llm2 = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,paired, p.adj),llm2.args)),
                    rai = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,p.adj),rai.args)),
                    spe = do.call(get(noquote(paste0("DA.",i))),c(list(count_table,outcome,relative,p.adj),spe.args)))
  
    res.sub[is.na(res.sub$pval),"pval"] <- 1
    res.sub[is.na(res.sub$pval.adj),"pval.adj"] <- 1
  
    return(res.sub)
  
  }
  if(cores != 1) stopCluster(cl)
  names(results) <- tests

  # Split ALDEx2 results in t.test and wilcoxon
  if("adx" %in% tests){
    adx.t <- as.data.frame(results["adx"])[,c(1:7,12)]
    adx.w <- as.data.frame(results["adx"])[,c(1:7,12)]
    colnames(adx.t) <- gsub("adx.","",colnames(adx.t))
    colnames(adx.w) <- colnames(adx.t)
    adx.t$pval <- as.numeric(as.data.frame(results["adx"])$adx.we.ep)
    adx.w$pval <- as.numeric(as.data.frame(results["adx"])$adx.wi.ep)
    adx.t$pval.adj <- p.adjust(adx.t$pval, method = p.adj)
    adx.w$pval.adj <- p.adjust(adx.w$pval, method = p.adj)
    adx.t$Method <- "ALDEx2 t-test"
    adx.w$Method <- "ALDEx2 wilcox"
    results["adx"] <- NULL
    res.names <- names(results)
    results <- c(results,list(adx.t),list(adx.w))
    names(results) <- c(res.names,"adx.t","adx.w")
  }

  return(results)

}