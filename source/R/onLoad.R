.onLoad <- function(libname, pkgname) {
  options("concordance.nCPU" = parallel::detectCores(),"concordance.nDraws"=1e4,"concordance.seed"=1,"concordance.verbose"=TRUE);
}
