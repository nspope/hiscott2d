#ifndef MVNINTERNAL_H
#define MVNINTERNAL_H

#include <RcppArmadillo.h> 

// [[Rcpp::depends(RcppArmadillo,mvtnorm)]]
// [[Rcpp::plugins(cpp11)]]

using ivec  = arma::ivec;
using uword = arma::uword;
using vec   = arma::vec;
using mat   = arma::mat;

vec normcdf (vec);

struct mvn_internal
{
  // hardcoded settings
  int    N      = 2,
         NU     = 0,
         INFORM = 0,
         RND    = 1,
         MAXPTS = 25000;
  double ABSEPS = 0.001,
         RELEPS = 0.;
  ivec   INFINF = {2, 2};
  vec    DELTA  = {0., 0.};

  mvn_internal (void)  
  {}

  double operator () (vec&, vec&, double);
  double operator () (vec&, vec&, double, ivec&);
};

#endif
