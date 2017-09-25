# Fix samr options problem

.onLoad <- function(libname, pkgname){
  options(error=utils::dump.frames)
}

