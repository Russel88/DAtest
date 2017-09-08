#' @export

prune.tests.DA <- function(tests, predictor, paired, relative){

  # Prune test argument if packages are not installed
  if(!"baySeq" %in% rownames(installed.packages())) tests <- tests[tests != "bay"]
  if(!"ALDEx2" %in% rownames(installed.packages())) tests <- tests[tests != "adx"] 
  if(!"MASS" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("neb")]
  if(!"lme4" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("neb")]
  if(!"edgeR" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("ere","erq")]
  if(!"metagenomeSeq" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("msf","zig")]
  if(!"DESeq2" %in% rownames(installed.packages())) tests <- tests[tests != "ds2"]
  if(!"limma" %in% rownames(installed.packages())) tests <- tests[tests != c("lim","lli","lli2")]
  if(!"RAIDA" %in% rownames(installed.packages())) tests <- tests[tests != "rai"]
  
  # Excluded tests that do not work with a paired argument
  if(!is.null(paired)){
    tests <- tests[!tests %in% c("bay","adx","ere","msf","zig","aov","lao","lao2","kru","rai","spe")]
  } 
  
  # Only include some tests if there are more than two levels in predictor
  if(length(levels(as.factor(predictor))) > 2){
    tests <- tests[tests %in% c("neb","erq","ds2","lim","lli","lli2","aov","lao","lao2","kru","lrm","llm","llm2","spe")]
  } else {
    # Excluded tests if levels in predictor is exactly 2
    tests <- tests[!tests %in% c("aov","lao","lao2","kru","lrm","llm","llm2","spe")]
  }
  
  # Only include specific tests if predictor is numeric
  if(is.numeric(predictor)){
    tests <- tests[tests %in% c("neb","erq","ds2","lim","lli","lli2","lrm","llm","llm2","spe")]
  } else {
    # Exclude if not numeric
    tests <- tests[!tests %in% c("spe")]
  }
  
  # Exclude if relative is false
  if(relative == FALSE){
    tests <- tests[!tests %in% c("ltt2","neb","erq","ere","msf","zig","bay","ds2","adx","lli2","lao2","llm2","rai")]
  }
  
  return(tests)
  
}




