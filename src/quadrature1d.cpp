#include <RcppArmadillo.h> 

// [[Rcpp::depends(RcppArmadillo,mvtnorm)]]
// [[Rcpp::plugins(cpp11)]]

using uword = arma::uword;
using vec   = arma::vec;
using mat   = arma::mat;

struct quadrature1d
{
  const mat::fixed<3,3> P = {{1, -3,  2},
                             {0,  4, -4},
                             {0,  -1, 2}};
  vec::fixed<3> M = {0., 0., 0.};

  const uword     lev = 3,
                  np  = 7;

  quadrature1d (uword lev = 3) :
    lev (lev),
    np  (2*lev + 1)
  {}

  vec get_weights (void)
  {
    return M;

  }

  void compute_weights (double mu, double l)
  {
    double ss = sqrt(l);
    double aa = -mu/ss;
    double bb = (1.-mu)/ss;
    M[0] = R::pnorm(bb,0,1,1,0) - R::pnorm(aa,0,1,1,0);
    M[1] = mu*M[0] + ss*(R::dnorm(aa,0,1,0) - R::dnorm(bb,0,1,0));
    M[2] = mu*M[1] + ss*ss*M[0] - ss*R::dnorm(bb,0,1,0);
    M = P * M;
  }

  vec node (vec& fjs, 
            vec& fi, 
            vec& fis, 
            double lbp, 
            double lb, 
            double scp, 
            double sc, 
            double lm, 
            double le, 
            bool root = false)
  {
    /* quadrature at node */

    uword i, j, ii;
    double mu, lcb, v;

    vec fj = arma::zeros<vec>(np);

    double mx = arma::max(fi);
    fi -= mx;
    fi  = arma::exp(fi) % fis;

    v      = 1./(2.*sc);
    lm     = le * v * v * lm;
    lbp   *= v;
    scp   *= v;
    lb    *= v;

    mu = lbp;
    for (i=0; i<np; ++i)
    {
      lcb = lb;
      for (j=0; j<lev; ++j)
      {
        compute_weights(mu - lcb, lm);
        for (ii=0; ii<3; ++ii)
          fj[i] += fi[j*2 + ii] * M[ii];
        lcb += 1.;
      }
      if (root) goto finish;
      mu += scp;
    }
    finish:

    fjs %= arma::sign(fj);
    fj   = arma::log(arma::abs(fj));
    fj  += mx;

    return fj;
  }

  vec tip (double lbp,
           double scp,
           uword  data,
           double lm,
           double si,
           double le,
           bool root = false)
  {
    //TODO not sure this is quite right
    uword i;
    double v, mu;

    vec fj = arma::zeros<vec>(np);

    v = 1./sqrt(le*lm + si);

    lbp *= v;
    scp *= v;

    mu = lbp;
    for (i=0; i<np; ++i)
    {
      fj[i] = R::pnorm(0., mu, 1., !data, true);//TODO
      if (root) goto finish;
      mu += scp;
    }
    finish:

    return fj;
  }
};

RCPP_MODULE(quadrature1d) {
  Rcpp::class_<quadrature1d>( "quadrature1d" )
    .constructor<uword>()
    .method( "calculate", &quadrature1d::compute_weights )
    .method( "weights",   &quadrature1d::get_weights )
    .method( "node",      &quadrature1d::node )
    .method( "tip",       &quadrature1d::tip )
    ;
}
