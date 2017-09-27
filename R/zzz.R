# Fix samr options problem

.onLoad <- function(libname, pkgname){
  options(error=NULL)
}

