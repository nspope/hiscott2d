
.onLoad <- function(libname, pkgname){
  Rcpp::loadModule("quadrature2d", TRUE)
  Rcpp::loadModule("quadrature1d", TRUE)
  Rcpp::loadModule("quadrature2drot", TRUE)
}
