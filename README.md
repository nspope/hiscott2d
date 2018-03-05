# hiscott2d

Calculate the log-likelihood for a two trait threshold model, with phylogenetic and sampling variance.

This package uses a modified form of the method of Hiscott et al. 2016. The quadrature rule used is a sort of Newton-Cotes with a Gaussian weighting function. The `level` argument determines the number of subintervals in one dimension: with `level=10` there are 100 2d intervals. Curse of dimensionality: it's slow -- implemented here mostly as a way to groundtruth the variational approximation used by the `epee` package.

Take the model given by,

```
set.seed(1)
Sigma  = diag(2)        # sampling covariance
Lambda = diag(2)+1      # phylogenetic covariance
mu     = c(-1, 1)       # ancestral trait "values"
tree   = ape::rcoal(25) # tree
```

Simulate some (coin-flip) data,

```
traits = matrix(sample(0:1, 50, replace=T), 2, 25)
```

The likelihood calculated using 15 intervals per dimension,

```
hiscott2d(data = traits, tree = tree, Lambda = Lambda, Sigma = Sigma, mu = mu, level = 15)
```

The approach used by `hiscott2d()` is to construct a quadrature grid aligned with the canonical basis. In practice, though, I don't think this makes much sense: I prefer to choose a basis such that the quadrature grid is aligned with the principle axes of the phylogenetic covariance. This should give the exact same results when the phylogenetic covariance is diagonal, but should be less wasteful when there is a non-zero correlation between the two traits. I've implemented this second approach as,

```
hiscott2drot(data = traits, tree = tree, Lambda = Lambda, Sigma = Sigma, mu = mu, level = 15)
```

Compare against Genz-Bretz,

```
genz2d(data = traits, tree = tree, Lambda = Lambda, Sigma = Sigma, mu = mu)
```


The univariate method is also included,

```
hiscott1d(data = traits[1,], tree = tree, Lambda = Lambda[1,1], Sigma = Sigma[1,1], mu = mu[1], level = 15)
genz1d(data = traits[1,], tree = tree, Lambda = Lambda[1,1], Sigma = Sigma[1,1], mu = mu[1])
```

Hiscott G, Fox C, Parry M, Bryant D. 2016. Efficient recycled algorithms for quantitative trait models on phylogenies. Genome Biology and Evolution 8: 1338-1350

