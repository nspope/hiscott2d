#include "mvn_internal.hpp"
#include <mvtnormAPI.h>

// [[Rcpp::depends(RcppArmadillo,mvtnorm)]]
// [[Rcpp::plugins(cpp11)]]

vec normcdf (vec input) // until RcppArmadillo 8.4 comes out
{
  for (auto& val : input)
    val = R::pnorm(val, 0., 1., 1., 0);
  return input;
}

double mvn_internal::operator () (vec&   LOWER,
                                  vec&   UPPER,
                                  double CORREL)
{
  return (*this)(LOWER, UPPER, CORREL, INFINF);
}

double mvn_internal::operator () (vec&   LOWER,
                                  vec&   UPPER,
                                  double CORREL,
                                  ivec&  INFIN)
{
  double value = 0.,
         error = 0.;

  mvtnorm_C_mvtdst (
      &N,            // dimension
      &NU,           // df parameter
      LOWER.memptr(),// normalized lower bound with -inf set to 0
      UPPER.memptr(),// normalized upper bound with inf set to 0
      INFIN.memptr(),// 2=both finite; -1=both infinite; 0=lower infinite; 1=upper infinite
      &CORREL,       // lower triangle of correlation matrix
      DELTA.memptr(),// location parameter
      &MAXPTS,       // maxpts (GenzBretz)
      &ABSEPS,       // abseps (GenzBretz)
      &RELEPS,       // releps (GenzBretz)
      &error,        // output
      &value,        // output
      &INFORM,       // inform (GenzBretz)
      &RND);         // rnd    (GenzBretz)

  return value;
}
