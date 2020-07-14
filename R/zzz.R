.onAttach <- function(libname, pkgname){
  packageStartupMessage("DAtest version 2.7.17")
  if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() == "4.0.0") {
    snow::setDefaultClusterOptions(setup_strategy = "sequential")
  }
}