// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


RcppExport SEXP _rcpp_module_boot_quadrature1d();
RcppExport SEXP _rcpp_module_boot_quadrature2d();
RcppExport SEXP _rcpp_module_boot_quadrature2drot();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_quadrature1d", (DL_FUNC) &_rcpp_module_boot_quadrature1d, 0},
    {"_rcpp_module_boot_quadrature2d", (DL_FUNC) &_rcpp_module_boot_quadrature2d, 0},
    {"_rcpp_module_boot_quadrature2drot", (DL_FUNC) &_rcpp_module_boot_quadrature2drot, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_hiscott2d(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
