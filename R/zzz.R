
.onLoad <- function(libname, pkgname){
  message("DAtest version 2.6.2")
  
  # Fix samr problem
  if(.Platform$OS.type == "windows"){
    if("samr" %in% rownames(installed.packages())){
      options(error = NULL)
    }
  }
}

