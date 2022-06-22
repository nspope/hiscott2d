#include <RcppArmadillo.h> 
#include "mvn_internal.hpp"

// [[Rcpp::depends(RcppArmadillo,mvtnorm)]]
// [[Rcpp::plugins(cpp11)]]

using ivec  = arma::ivec;
using uword = arma::uword;
using vec   = arma::vec;
using mat   = arma::mat;

struct quadrature2drot
{
  /* 2d quadrature but with grid rotated such that it aligns with principle axes */

  const mat::fixed<9,9> P = {{1, -3,  2, -3,   9,  -6,  2,  -6,  4},
                             {0,  4, -4,  0, -12,  12,  0,   8, -8},
                             {0, -1,  2,  0,   3,  -6,  0,  -2,  4},
                             {0,  0,  0,  4, -12,   8, -4,  12, -8},
                             {0,  0,  0,  0,  16, -16,  0, -16, 16},
                             {0,  0,  0,  0,  -4,   8,  0,   4, -8},
                             {0,  0,  0, -1,   3,  -2,  2,  -6,  4},
                             {0,  0,  0,  0,  -4,   4,  0,   8, -8},
                             {0,  0,  0,  0,   1,  -2,  0,  -2,  4}};

  mat::fixed<2,3> Q;
  vec::fixed<9>   M = arma::zeros<vec>(9);
  vec::fixed<2>   m1, m2, a, b;

  const uword     lev = 3,
                  np  = 7;

  const double    err,
                  glb,
                  gsc;

  mvn_internal pmvnorm;

  quadrature2drot (uword lev, uword tips) :
    lev (lev),
    np  (2*lev + 1),
    err (pow(double(np-1),-4)*double(tips)/(2.*double(tips) - 1)),
    glb ( R::qnorm(err,0,1,1,0)),
    gsc ( (R::qnorm(err,0,1,0,0) - glb)/double(np-1))
  {}

  vec get_weights ()
  {
    return M;
  }

  void compute_weights (const vec& m, const double s)
  {
           m1 = arma::shift(m,1);
    double ss = sqrt(s);

    a = -m1 / ss;
    b = (1.-m1) / ss;

    Q.col(0) = normcdf(b) - normcdf(a);

    a = arma::normpdf(a);
    b = arma::normpdf(b);
    m2 = a - b;

    Q.col(1) = m1 % Q.col(0) + ss * m2;
    Q.col(2) = m1 % Q.col(1) + ss * (ss * Q.col(0) - b);

    // compute moments
    M[0] = arma::prod(Q.col(0));
    M[1] = M[0]*m[0] + ss*Q[0]*m2[1];
    M[2] = M[1]*m[0] + s*M[0] - ss*b[1]*Q[0];
    M[3] = M[0]*m[1] + ss*Q[1]*m2[0];
    M[4] = M[3]*m[0] + ss*Q[2]*m2[1];
    M[5] = M[4]*m[0] + s*M[3] - ss*b[1]*Q[2];
    M[6] = M[3]*m[1] + s*M[0] - ss*b[0]*Q[1];
    M[7] = M[6]*m[0] + ss*Q[4]*m2[1];
    M[8] = M[7]*m[0] + s*M[6] - ss*b[1]*Q[4];

    M = P * M;
  }

  vec node (vec&   fjs, // stores result (sign)
            vec&   fi,  // function evaluations
            vec&   fis, // sign of function evals
            double hp,  // height of parent node
            double h,   // height of child node
            vec&   gmu, // global mean 
            bool   root = false)   
  {      
    /* calculate weights and perform quadrature */

    if (fjs.n_elem != np*np ||
        fi.n_elem  != np*np ||
        fis.n_elem != np*np ||
        gmu.n_elem != 2     )
      Rcpp::stop("Wrong input dimensions");

    vec::fixed<2> lcb, // lower centroid bounds
                  mu;  // grid point of parent

    vec fj = arma::zeros<vec>(np*np);

    // underflow control
    double mx = arma::max(fi);
    fi -= mx;
    fi  = arma::exp(fi) % fis;

    // scaling
    vec lbp = 0.5 * glb * arma::ones<vec>(2) * sqrt(hp)/sqrt(h) / gsc; 
    vec  lb = 0.5 * glb * arma::ones<vec>(2) / gsc;
    double le  = h - hp;
    double scp = 0.5 * sqrt(hp)/sqrt(h);
    double  va = 0.25 * le/h * 1./(gsc * gsc);

    mu [1] = lbp [1];
    for (uword k=0; k<np; ++k)   //loop over grid points (parent)
    {
      mu [0] = lbp [0];
      for (uword l=0; l<np; ++l) //loop over grid points (parent)
      {
        lcb [1] = lb [1];
        for (uword i=0; i<lev; ++i)   //loop over grid cells
        {
          lcb [0] = lb [0];
          for (uword j=0; j<lev; ++j) //loop over grid cells
          {
            compute_weights(mu - lcb, va);
            for (uword ii=0; ii<3; ++ii) //loop over points in a cell
              for (uword jj=0; jj<3; ++jj)
                fj[np*k + l] += fi[ii*np + jj + 2*(np*i + j)] * M[ii*3 + jj];
            lcb [0] += 1.;
          }
          lcb [1] += 1.; 
        }
        if (root) goto finish;
        mu [0] += scp;
      }
      mu [1] += scp; 
    }

    finish:

    fjs %= arma::sign(fj);
    fj   = arma::log(arma::abs(fj));
    fj  += mx;

    return fj;
  }

  vec tip (double hp,   // height of parent
           double h,    // height of child
           vec&   gmu,  // global lower bound
           mat&   lsqrt,// sqrt of phylogenetic covariance matrix
           mat&   ssqrt,// srqt of error covariance matrix
           ivec&  data, // data for tip
           bool   root = false)   
  {      
    /* calculate tip probabilities */

    if (lsqrt.n_rows != 2 ||
        lsqrt.n_cols != 2 ||
        ssqrt.n_rows != 2 ||
        ssqrt.n_cols != 2 ||
        gmu.n_elem   != 2 ||
        data.n_elem != 2  ) 
      Rcpp::stop("Wrong input dimensions");

    vec fj = arma::zeros<vec>(np*np);

    double le = h - hp;                 // length of edge

    mat    cv = le * lsqrt * lsqrt.t() + ssqrt * ssqrt.t(); // covariance

    vec a  = {0., 0.}, // lower bounds (-inf set to 0)
        b  = {0., 0.}, // upper bounds (inf set to 0)
        mu = {0., 0.}; // mean in parent grid

    vec     v  = 1./arma::sqrt(cv.diag()); // sqrt precision
    double  co = arma::prod(v)*cv(1,0);    // correlation

    // scaling
    vec sci = {0., gsc},
        scj = {gsc, 0.},
        lbp = {glb, glb};
        sci = (sqrt(hp) * lsqrt * sci) % v;
        scj = (sqrt(hp) * lsqrt * scj) % v;
        lbp = (gmu + sqrt(hp) * lsqrt * lbp) % v;

    for (uword i=0; i<np; ++i)   //loop over grid points (parent)
    {
      mu = lbp + sci*double(i);
      for (uword j=0; j<np; ++j) //loop over grid points (parent)
      {
        for (uword k=0; k<2; ++k)
        {
          if (data[k]) a[k] = -mu[k];
          else         b[k] = -mu[k];
        }

        fj[np*i + j] = std::max(0., pmvnorm(a, b, co, data));

        if (root) goto finish;
        mu += scj;
      }
    }
    finish:

    return arma::log(fj);
  }

  vec tip2 (double hp,   // height of parent
            double h,    // height of child
            vec&   gmu,  // global lower bound
            mat&   lsqrt,// sqrt of phylogenetic covariance matrix
            mat&   ssqrt,// srqt of error covariance matrix
            int    cat,  // data for tip (binary)
            double cont, // data for tip (continuous)
            bool   root = false)   
  {      
    /* calculate tip probabilities given binary and continuous trait */

    if (lsqrt.n_rows != 2 ||
        lsqrt.n_cols != 2 ||
        ssqrt.n_rows != 2 ||
        ssqrt.n_cols != 2 ||
        gmu.n_elem   != 2 )
    {
      Rcpp::stop("Wrong input dimensions");
    }

    vec fj = arma::zeros<vec>(np*np);

    double le = h - hp;                 // length of edge

    mat cv = le * lsqrt * lsqrt.t() + ssqrt * ssqrt.t(); // covariance
    vec mu = {0., 0.}; // mean in parent grid

    vec     v  = 1./arma::sqrt(cv.diag()); // sqrt precision
    double  co = arma::prod(v)*cv(1,0),    // correlation
            vcat = sqrt(cv.at(0,0)),
            vcon = sqrt(cv.at(1,1));
    double  muc, vac;

    // scaling
    vec sci = {0., gsc},
        scj = {gsc, 0.},
        lbp = {glb, glb};
        sci = (sqrt(hp) * lsqrt * sci) % v;
        scj = (sqrt(hp) * lsqrt * scj) % v;
        lbp = (gmu + sqrt(hp) * lsqrt * lbp) % v;

    for (uword i=0; i<np; ++i)   //loop over grid points (parent)
    {
      mu = lbp + sci*double(i);
      for (uword j=0; j<np; ++j) //loop over grid points (parent)
      {
        muc = mu[0] + co*vcat/vcon*(cont - mu[1]);
        vac = sqrt((1. - co*co)*vcat*vcat);

        fj[np*i + j] = std::max(0.0, R::pnorm(0, muc, vac, 1-cat, 0) * R::dnorm(cont, mu[1], vcon, 0));

        if (root) goto finish;
        mu += scj;
      }
    }
    finish:

    return arma::log(fj);
  }
};

RCPP_MODULE(quadrature2drot) {
  Rcpp::class_<quadrature2drot>( "quadrature2drot" )
    .constructor<uword, uword>()
    .method( "calculate", &quadrature2drot::compute_weights )
    .method( "weights",   &quadrature2drot::get_weights )
    .method( "node",      &quadrature2drot::node )
    .method( "tip",       &quadrature2drot::tip )
    .method( "tip2",      &quadrature2drot::tip2 )
    ;
}

