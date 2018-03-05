#include <RcppArmadillo.h> 
#include <mvtnormAPI.h>

// [[Rcpp::depends(RcppArmadillo,mvtnorm)]]
// [[Rcpp::plugins(cpp11)]]

using ivec  = arma::ivec;
using uword = arma::uword;
using vec   = arma::vec;
using mat   = arma::mat;

vec normcdf (vec input) // until RcppArmadillo 8.4 comes out
{
  for (auto& val : input)
    val = R::pnorm(val, 0., 1., 1., 0);
  return input;
}

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

  // compute probability
  double operator () (vec&   LOWER,
                      vec&   UPPER,
                      double CORREL)
  {
    return (*this)(LOWER, UPPER, CORREL, INFINF);
  }

  double operator () (vec&   LOWER,
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
};

struct quadrature2d
{
  const mat::fixed<9,9> P = {{1, -3,  2, -3,   9,  -6,  2,  -6,  4},
                             {0,  4, -4,  0, -12,  12,  0,   8, -8},
                             {0, -1,  2,  0,   3,  -6,  0,  -2,  4},
                             {0,  0,  0,  4, -12,   8, -4,  12, -8},
                             {0,  0,  0,  0,  16, -16,  0, -16, 16},
                             {0,  0,  0,  0,  -4,   8,  0,   4, -8},
                             {0,  0,  0, -1,   3,  -2,  2,  -6,  4},
                             {0,  0,  0,  0,  -4,   4,  0,   8, -8},
                             {0,  0,  0,  0,   1,  -2,  0,  -2,  4}};

  mat::fixed<2,3> G, H;
  vec::fixed<9>   M = arma::zeros<vec>(9);
  vec::fixed<2>   m1, m2, ss, ga, gb, ha, hb;

  const uword     lev = 3,
                  np  = 7;

  const vec::fixed<2> zeros = {0., 0.},
                      ones  = {1., 1.};

  mvn_internal pmvnorm;

  const double tol = 0;

  quadrature2d (uword lev = 3, double tol = 0) :
    tol (tol),
    lev (lev),
    np  (2*lev + 1)
  {}

  vec get_weights ()
  {
    return M;
  }

  void compute_weights (const vec& m, const vec& s, double c)
  {
    // lower, upper bounds; correlation coefficient
    m1   = -m/arma::sqrt(s); 
    m2   = (1.-m)/arma::sqrt(s);
    M[0] = pmvnorm(m1, m2, c/arma::prod(arma::sqrt(s)));

//    if (M[0] < tol)
//    {
//      M.zeros();
//      return;
//    }

    G.zeros();
    H.zeros();

    m1 = arma::shift(m,1) - c * m / s;
    m2 = m1 + c / s;
    ss = arma::sqrt(arma::shift(s,1) - c * c / s);

    ga = -m1 / ss;
    gb = (1.-m1) / ss;
    ha = -m2 / ss;
    hb = (1.-m2) / ss;

    G.col(0) = normcdf(gb) - normcdf(ga);
    G.col(1) = m1 % G.col(0) + ss % (arma::normpdf(ga) - arma::normpdf(gb));
    G.col(2) = m1 % G.col(1) + ss % (ss % G.col(0) - arma::normpdf(gb));
    G.each_col() %= arma::normpdf(zeros, m, arma::sqrt(s));

    H.col(0) = normcdf(hb) - normcdf(ha);
    H.col(1) = m2 % H.col(0) + ss % (arma::normpdf(ha) - arma::normpdf(hb));
    H.col(2) = m2 % H.col(1) + ss % (ss % H.col(0) - arma::normpdf(hb));
    H.each_col() %= arma::normpdf(ones, m, arma::sqrt(s));

    H *= -1.;
    G += H;

    // compute moments
    M[1] = M[0]*m[0] + s[0]*G[0] + c*G[1];
    M[2] = M[1]*m[0] + s[0]*(M[0] + H[0]) + c*G[3];
    M[3] = M[0]*m[1] + c*G[0] + s[1]*G[1];
    M[4] = M[3]*m[0] + s[0]*G[2] + c*(M[0] + H[1]);
    M[5] = M[4]*m[0] + s[0]*(M[3] + H[2]) + c*(M[1] + H[3]);
    M[6] = M[3]*m[1] + c*G[2] + s[1]*(M[0] + H[1]);
    M[7] = M[6]*m[0] + s[0]*G[4] + c*(2.*M[3] + H[1]);
    M[8] = M[7]*m[0] + s[0]*(M[6] + H[4]) + c*(2.*M[4] + H[3]);

    M = P * M;

    if (tol > 0) // concievably a small accuracy sacrifice for a lot less storage
      M.elem(arma::find(arma::abs(M) < tol)).zeros();
  }

  mat node1 (vec    lbp, // lower bound of parent
             vec    lb,  // lower bound of child
             vec    scp, // step size for parent
             vec    sc,  // step size for child
             vec    lm,  // variances
             double co,  // covariance
             double le,  // length of edge
             bool   root = false)   
  {      
    /* calculate quadrature weights and return as matrix */

    if (lbp.n_elem != 2 ||
        lb.n_elem  != 2 ||
        scp.n_elem != 2 ||
        sc.n_elem  != 2 ||
        lm.n_elem  != 2 )
      Rcpp::stop("Wrong input dimensions");

    uword k, l, i, j, ii, jj, row, col;

    vec::fixed<2> lcb, // lower centroid bounds
                  mu;  // grid point of parent

    mat w = arma::zeros<mat>(np*np, np*np);

    // set up covariance
    vec v  = 1./(2.*sc);
    lm     = le * v % v % lm;
    co     = le * arma::prod(v) * co;

    lbp   %= v;
    scp   %= v;
    lb    %= v;

    mu [1] = lbp [1];
    for (k=0; k<np; ++k)   //loop over grid points (parent)
    {
      mu [0] = lbp [0];
      for (l=0; l<np; ++l) //loop over grid points (parent)
      {
        row     = np*k + l;
        lcb [1] = lb [1];
        for (i=0; i<lev; ++i)   //loop over grid cells
        {
          lcb [0] = lb [0];
          for (j=0; j<lev; ++j) //loop over grid cells
          {
            col = 2*(np*i + j);
            compute_weights(mu - lcb, lm, co);
            for (ii=0; ii<3; ++ii) //loop over points in a cell
              for (jj=0; jj<3; ++jj)
                w.at(row, ii*np + jj + col) += M[ii*3 + jj];
            lcb [0] += 1.;
          }
          lcb [1] += 1.; 
        }
        if (root) goto finish;
        mu [0] += scp [1];
      }
      mu [1] += scp [0]; 
    }
    finish:

    return w;
  }

  vec node2 (vec&   fjs, // stores result (sign)
             vec&   fi,  // function evaluations
             vec&   fis, // sign of function evals
             vec    lbp, // lower bound of parent
             vec    lb,  // lower bound of child
             vec    scp, // step size for parent
             vec    sc,  // step size for child
             vec    lm,  // variances
             double co,  // covariance
             double le,  // length of edge
             bool   root = false)   
  {      
    /* calculate weights and perform quadrature */

    if (fjs.n_elem != np*np ||
        fi.n_elem  != np*np ||
        fis.n_elem != np*np ||
        lbp.n_elem != 2     ||
        lb.n_elem  != 2     ||
        scp.n_elem != 2     ||
        sc.n_elem  != 2     ||
        lm.n_elem  != 2     )
      Rcpp::stop("Wrong input dimensions");

    vec::fixed<2> lcb, // lower centroid bounds
                  mu;  // grid point of parent

    vec fj = arma::zeros<vec>(np*np);

    double mx = arma::max(fi);
    fi -= mx;
    fi  = arma::exp(fi) % fis;

    // set up covariance
    vec v  = 1./(2.*sc);
    lm     = le * v % v % lm;
    co     = le * arma::prod(v) * co;

    lbp   %= v;
    scp   %= v;
    lb    %= v;

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
            compute_weights(mu - lcb, lm, co);
            for (uword ii=0; ii<3; ++ii) //loop over points in a cell
              for (uword jj=0; jj<3; ++jj)
                fj[np*k + l] += fi[ii*np + jj + 2*(np*i + j)] * M[ii*3 + jj];
            lcb [0] += 1.;
          }
          lcb [1] += 1.; 
        }
        if (root) goto finish;
        mu [0] += scp [1];
      }
      mu [1] += scp [0]; 
    }
    finish:

    fjs %= arma::sign(fj);
    fj   = arma::log(arma::abs(fj));
    fj  += mx;

    return fj;
  }

  vec node3 (vec&   fjs, // stores result (sign)
             vec&   fi,  // function evaluations
             vec&   fis, // sign of function evals
             mat&   w,   // quadrature weights
             bool   root = false)   
  {      
    /* perform quadrature with precomputed weights */

    if (fjs.n_elem != np*np ||
        fi.n_elem  != np*np ||
        fis.n_elem != np*np ||
        w.n_rows   != np*np ||
        w.n_cols   != np*np )
      Rcpp::stop("Wrong input dimensions");

    double mx = arma::max(fi);
    fi -= mx;
    fi  = arma::exp(fi) % fis;

    vec fj = w * fi;

    fjs %= arma::sign(fj);
    fj   = arma::log(arma::abs(fj));
    fj  += mx;

    return fj;
  }

  vec tip (vec    lbp,  // lower bound of parent
           vec    scp,  // step size for parent
           ivec&  data, // data for tip
           vec    lm,   // variances (phylogenetic)
           vec    si,   // variances (error)
           double co,   // covariance (phylogentic)
           double ci,   // covariance (error)
           double le,   // length of edge
           bool   root = false)   
  {      
    /* calculate tip probabilities */

    if (lbp.n_elem  != 2     ||
        scp.n_elem  != 2     ||
        data.n_elem != 2     ||
        si.n_elem   != 2     ||
        lm.n_elem   != 2     )
      Rcpp::stop("Wrong input dimensions");

    vec fj = arma::zeros<vec>(np*np);

    vec a  = {0., 0.},
        b  = {0., 0.},
        mu = {0., 0.};

    vec v  = 1./arma::sqrt(le*lm + si);
        co = arma::prod(v)*(le*co + ci);

    lbp %= v;
    scp %= v;

    mu [1] = lbp [1];
    for (uword k=0; k<np; ++k)   //loop over grid points (parent)
    {
      mu [0] = lbp [0];
      for (uword l=0; l<np; ++l) //loop over grid points (parent)
      {
        for (uword i=0; i<2; ++i)
        {
          if (data[i]) a[i] = -mu[i];
          else         b[i] = -mu[i];
        }
        fj[np*k + l] = std::max(0., pmvnorm(a, b, co, data));
        if (root) goto finish;
        mu [0] += scp [1];
      }
      mu [1] += scp [0]; 
    }
    finish:

    return arma::log(fj);
  }
};

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
};

RCPP_MODULE(quadrature2d) {
  Rcpp::class_<quadrature2d>( "quadrature2d" )
    .constructor<uword, double>()
//  .method( "calculate", &weights::compute_weights )
//  .method( "weights",   &weights::weight )
    .method( "node1",     &quadrature2d::node1 )
    .method( "node2",     &quadrature2d::node2 )
    .method( "node3",     &quadrature2d::node3 )
    .method( "tip",       &quadrature2d::tip )
    ;
}

RCPP_MODULE(quadrature2drot) {
  Rcpp::class_<quadrature2drot>( "quadrature2drot" )
    .constructor<uword, uword>()
    .method( "calculate", &quadrature2drot::compute_weights )
    .method( "weights",   &quadrature2drot::get_weights )
    .method( "node",      &quadrature2drot::node )
    .method( "tip",       &quadrature2drot::tip )
    ;
}

