#' Prune tests argument for \code{testDA}
#'
#' @param tests A character vector with names of tests
#' @param predictor Vector with predictor variable
#' @param paired Factor with paired variable
#' @param covars Named list with covariables
#' @param relative Include tests that work with relative abundances (TRUE) or only absolute abundances (FALSE)
#' @param decimal Exclude tests that do not work with decimals (TRUE)
#' @export

prune.tests.DA <- function(tests, predictor, paired, covars, relative, decimal){

  # Prune test argument if packages are not installed
  if(!"baySeq" %in% rownames(installed.packages())) tests <- tests[tests != "bay"]
  if(!"ALDEx2" %in% rownames(installed.packages())) tests <- tests[tests != "adx"] 
  if(!"edgeR" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("ere","erq","ere2","erq2")]
  if(!"metagenomeSeq" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("msf","zig")]
  if(!"DESeq2" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("ds2","ds2x")]
  if(!"limma" %in% rownames(installed.packages())) tests <- tests[tests != c("lim","lli","lli2","vli")]
  if(!"statmod" %in% rownames(installed.packages())) tests <- tests[tests != c("lim","lli","lli2","vli")]
  if(!"RAIDA" %in% rownames(installed.packages())) tests <- tests[tests != "rai"]
  if(!"pscl" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("zpo","znb")]
  if(!"ancom.R" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("anc")]
  if(!"samr" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("sam")]
  if(!"mvabund" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("mvabund")]
  
  # Exclude tests that do not work with a paired argument
  if(!is.null(paired)){
    tests <- tests[!tests %in% c("qpo","zpo","znb","bay","adx","ere","ere2","msf","aov","lao","lao2","kru","rai","spe","pea")]
    # Exclude tests that only work with one value for each combination of predictor and paired arguments
    if(!all(table(paired,predictor) == 1)){
      tests <- tests[!tests %in% c("ttt","ltt","ltt2","wil","per","fri","qua","sam")]
    }
  } else {
    # Exclude if there is no paired
    tests <- tests[!tests %in% c("fri","qua")]
  }
  
  # Only include some tests if there are more than two levels in predictor
  if(length(levels(as.factor(predictor))) > 2){
    tests <- tests[tests %in% c("mva","sam","anc","qua","fri","znb","zpo","vli","qpo","poi","neb","erq","erq2","ds2","ds2x","lim","lli","lli2","aov","lao","lao2","kru","lrm","llm","llm2","spe","pea","zig")]
    # Exclude if only works for two-class paired
    if(!is.null(paired)){
      tests <- tests[!tests %in% c("sam")]
    } 
  } else {
    # Excluded tests if levels in predictor is exactly 2
    tests <- tests[!tests %in% c("aov","lao","lao2","kru","spe","pea","fri","qua","lrm","llm","llm2")]
  }
  
  # Only include specific tests if predictor is numeric
  if(is.numeric(predictor)){
    tests <- tests[tests %in% c("mva","sam","znb","zpo","vli","qpo","poi","neb","erq","erq2","lim","lli","lli2","lrm","llm","llm2","spe","pea")]
  } else {
    # Exclude if not numeric
    tests <- tests[!tests %in% c("spe","pea")]
  }
  
  # Exclude if relative is false
  if(relative == FALSE){
    tests <- tests[!tests %in% c("sam","anc","vli","ltt2","erq","ere","ere2","erq2","msf","zig","bay","ds2","ds2x","adx","lli2","lao2","llm2","rai")]
  } else {
    # Exclude if relative is TRUE
    tests <- tests[!tests %in% c("lrm","lim")]
  }
  
  # Exclude if decimal is TRUE
  if(decimal){
    tests <- tests[!tests %in% c("znb","zpo","qpo","poi","neb","mva")]
  }
  
  # Only include if covars are present
  if(!is.null(covars)){
    tests <- tests[tests %in% c("mva","znb","zpo","vli","qpo","poi","ds2","ds2x","neb","erq","erq2","zig","lrm","llm","llm2","lim","lli","lli2","aov","lao","lao2")]
  }
  
  return(tests)
  
}




