pruneTests <- function(tests, predictor, paired, covars, relative, decimal, zeroes){

  # Prune test argument if packages are not installed
  if(!"baySeq" %in% rownames(installed.packages())) tests <- tests[tests != "bay"]
  if(!"ALDEx2" %in% rownames(installed.packages())) tests <- tests[tests != "adx"] 
  if(!"edgeR" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("ere","erq","ere2","erq2","vli")]
  if(!"metagenomeSeq" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("msf","zig")]
  if(!"DESeq2" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("ds2","ds2x")]
  if(!"limma" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("lim","lli","lli2","vli","lia","lic")]
  if(!"statmod" %in% rownames(installed.packages()) && is.null(paired)) tests <- tests[!tests %in% c("lim","lli","lli2","vli","lia","lic")]
  if(!"pscl" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("zpo","znb")]
  if(!"samr" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("sam")]
  if(!"mvabund" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("mva")]
  if(!"ANCOMBC" %in% rownames(installed.packages())) tests <- tests[!tests %in% c("abc")]
  
  # Exclude tests that do not work with a paired argument
  if(!is.null(paired)){
    tests <- tests[!tests %in% c("abc","qpo","zpo","znb","bay","adx","ere","ere2","msf","aov","lao","lao2","aoa","aoc","kru","rai","spe","pea")]
    # Exclude tests that only work with one value for each combination of predictor and paired arguments
    if(!all(table(paired,predictor) == 1)){
      tests <- tests[!tests %in% c("ttt","ttr","ltt","ltt2","wil","per","fri","qua","sam","tta","ttc")]
    }
    # Exclude if too few levels
    if(length(unique(paired)) < 5){
      tests <- tests[!tests %in% c("lrm","llm","llm2","lim","lli","lli2","vli","neb","poi","zig","lma","lmc","lia","lic")]
    }
    
  } else {
    # Exclude if there is no paired
    tests <- tests[!tests %in% c("fri","qua")]
  }
  
  # Only include some tests if there are more than two levels in predictor
  if(length(unique(predictor)) > 2){
    tests <- tests[tests %in% c("abc","bay","sam","qua","fri","znb","zpo","vli","qpo","poi","neb","erq","erq2","ds2","ds2x","lim","lli","lli2","aov","lao","lao2","aoa","aoc","kru","lrm","llm","llm2","spe","pea","zig","lma","lmc","lia","lic")]
    # Exclude if only works for two-class paired
    if(!is.null(paired)){
      tests <- tests[!tests %in% c("sam")]
    } 
  } else {
    # Excluded tests if levels in predictor is exactly 2
    tests <- tests[!tests %in% c("aov","lao","lao2","aoa","aoc","kru","spe","pea","fri","qua","lrm","llm","llm2","lma","lmc")]
  }
  
  # Only include specific tests if predictor is numeric
  if(is.numeric(predictor)){
    tests <- tests[tests %in% c("abc","mva","sam","znb","zpo","vli","qpo","poi","neb","erq","erq2","lim","lli","lli2","lrm","llm","llm2","spe","pea","lma","lmc","lia","lic")]
  } else {
    # Exclude if not numeric
    tests <- tests[!tests %in% c("spe","pea")]
  }
  
  # Exclude if relative is false
  if(relative == FALSE){
    tests <- tests[!tests %in% c("abc","sam","vli","ltt2","erq","ere","ere2","erq2","msf","zig","bay","ds2","ds2x","adx","lli2","lao2","aoa","aoc","llm2","rai","tta","ttc","lma","lmc","lia","lic")]
  } else {
    # Exclude if relative is TRUE
    tests <- tests[!tests %in% c("lrm","lim")]
  }
  
  # Exclude if decimal is TRUE
  if(decimal){
    tests <- tests[!tests %in% c("abc","znb","zpo","qpo","poi","neb","mva")]
  }
  
  # Only include if covars are present
  if(!is.null(covars)){
    tests <- tests[tests %in% c("mva","znb","zpo","vli","qpo","poi","ds2","ds2x","neb","erq","erq2","zig","lrm","llm","llm2","lim","lli","lli2","aov","lao","lao2","aoa","aoc","lma","lmc","lia","lic")]
  }
  
  # Exclude if no zeroes
  if(!zeroes){
    tests <- tests[!tests %in% c("znb","zpo")]
  }
  
  return(tests)
  
}




