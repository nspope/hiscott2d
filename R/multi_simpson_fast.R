hiscott2d <- function(data, tree, Lambda=diag(2), Sigma=diag(2), mu=rep(0,2), level=3, weights=NULL)
{
  # Calculate loglikelihood with multivariate Hiscott algorithm.

  if (nrow(Lambda) != 2 || ncol(Lambda) != 2 ||
      nrow(Sigma)  != 2 || ncol(Sigma)  != 2 ||
      length(mu)   != 2 || nrow(data)   != 2 ||
      ncol(data)   != length(tree$tip.label) )
    stop ("Dimension mismatch")

  if (!is.null(weights) && class(weights) == "array")
    if (!all(dim(weights) == c((level*2+1)^2, (level*2+1)^2, tree$Nnode)))
      stop ("Weights of wrong dimension")

  ww     = new(quadrature2d, level)
  tree   = ape::reorder.phylo(tree, "postorder")
  nnodes = tree$Nnode
  ntips  = length(tree$tip.label)
  npoint = level*2 + 1
  len    = tree$edge.length
  root   = tree$edge[nrow(tree$edge),1]

  if (!is.null(weights) && class(weights) != "array")
    Fj = array(0, dim = c(npoint^2, npoint^2, nnodes)) # weights
  else
  {
    Fj  = matrix(0, npoint^2, nnodes) # log'd function evaluations
    Fjs = matrix(1, npoint^2, nnodes) # sign of function evaluations
  }

  # bounds
  hei = ape::node.depth.edgelength(tree)
  Lb  = matrix(0, 2, nnodes)
  Sc  = matrix(0, 2, nnodes)
  err = (npoint-1)^-4 * ntips/(2*ntips-1)
  sdv = sqrt(diag(Lambda))
  for (node in 1:nnodes)
  {
    Lb[,node] = qnorm(err, mu, sqrt(hei[node+ntips])*sdv)
    Sc[,node] = (qnorm(1-err, mu, sqrt(hei[node+ntips])*sdv) - Lb[,node])/(npoint-1)
  }

  # quadrature
  for (edge in 1:nrow(tree$edge))
  {
    i = tree$edge[edge,2]
    j = tree$edge[edge,1]-ntips

    if (i>ntips) 
    # internal node
    {
      if (!is.null(weights) && class(weights) == "array") 
        # weights supplied
        Fj[,j] = Fj[,j] + 
          ww$node3(Fjs[,j], 
                   Fj[,i-ntips], 
                   Fjs[,i-ntips],
                   weights[,,i-ntips],
                   j==root-ntips)
      else if (!is.null(weights) && class(weights) != "array") 
        # weights calculated with no quadrature
        Fj[,,i-ntips] = 
          ww$node1(Lb[,j],
                   Lb[,i-ntips],
                   Sc[,j],
                   Sc[,i-ntips],
                   diag(Lambda),
                   Lambda[1,2],
                   len[edge],
                   j==root-ntips)
      else 
        # weights calculated with quadrature
        Fj[,j] = Fj[,j] + 
          ww$node2(Fjs[,j], 
                   Fj[,i-ntips], 
                   Fjs[,i-ntips],
                   Lb[,j],
                   Lb[,i-ntips],
                   Sc[,j],
                   Sc[,i-ntips],
                   diag(Lambda),
                   Lambda[1,2],
                   len[edge],
                   j==root-ntips)
    } 
    else 
    # tip node
    {
      if (!is.null(weights) && class(weights) != "array") {}
      else 
        Fj[,j] = Fj[,j] + 
          ww$tip(Lb[,j],
                 Sc[,j],
                 data[,i], 
                 diag(Lambda), 
                 diag(Sigma), 
                 Lambda[1,2], 
                 Sigma[1,2], 
                 len[edge],
                 j==root-ntips)
    }
  }

  if (!is.null(weights) && class(weights) != "array")
    return (Fj) # return weights
  else
    return (Fj[1,j]) # return loglikelihood
}

genz2d <- function(data, tree, Lambda=diag(2), Sigma=diag(2), mu=rep(0,2))
{
  # Calculate loglikelihood with Genz-Bretz algorithm.

  if (nrow(Lambda) != 2 || ncol(Lambda) != 2 ||
      nrow(Sigma)  != 2 || ncol(Sigma)  != 2 ||
      length(mu)   != 2 || nrow(data)   != 2 ||
      ncol(data)   != length(tree$tip.label) )
    stop ("Dimension mismatch")

  tree   = ape::reorder.phylo(tree, "postorder")
  lower  = log(c(data))
  upper  = -log(1-c(data))
  W      = kronecker(ape::vcv(tree), Lambda) + kronecker(diag(ncol(data)), Sigma)
  m      = c(kronecker(rep(1,ncol(data)), mu))

  out    = mvtnorm::pmvnorm(lower, upper, mean=m, sigma=W)
  attr(out, "error") = log(out + c(-attr(out, "error"),attr(out, "error")))
  log(out)
}

hiscott1d <- function(data, tree, Lambda=1, Sigma=1, mu=0, level=3)
{
  # Calculate loglikelihood with univariate Hiscott algorithm.

  if (length(Lambda) != 1 || 
      length(Sigma)  != 1 || 
      length(mu)     != 1 || 
      length(data)   != length(tree$tip.label) )
    stop ("Dimension mismatch")

  ww     = new(quadrature1d, level)
  tree   = ape::reorder.phylo(tree, "postorder")
  nnodes = tree$Nnode
  ntips  = length(tree$tip.label)
  npoint = level*2 + 1
  len    = tree$edge.length
  root   = tree$edge[nrow(tree$edge),1]

  Fj  = matrix(0, npoint, nnodes) # log'd function evaluations
  Fjs = matrix(1, npoint, nnodes) # sign of function evaluations

  # bounds
  hei = ape::node.depth.edgelength(tree)
  Lb  = matrix(0, 1, nnodes)
  Sc  = matrix(0, 1, nnodes)
  err = (npoint-1)^-4 * ntips/(2*ntips-1)
  sdv = sqrt(diag(Lambda))
  for (node in 1:nnodes)
  {
    Lb[node] = qnorm(err, mu, sqrt(hei[node+ntips])*sdv)
    Sc[node] = (qnorm(1-err, mu, sqrt(hei[node+ntips])*sdv) - Lb[node])/(npoint-1)
  }

  # quadrature
  for (edge in 1:nrow(tree$edge))
  {
    i = tree$edge[edge,2]
    j = tree$edge[edge,1]-ntips

    if (i>ntips) 
    # internal node
    {
      Fj[,j] = Fj[,j] + 
        ww$node(Fjs[,j], 
                Fj[,i-ntips], 
                Fjs[,i-ntips],
                Lb[j],
                Lb[i-ntips],
                Sc[j],
                Sc[i-ntips],
                Lambda,
                len[edge],
                j==root-ntips)
    } 
    else 
    # tip node
    {
      Fj[,j] = Fj[,j] + 
        ww$tip(Lb[j],
               Sc[j],
               data[i], 
               Lambda, 
               Sigma, 
               len[edge],
               j==root-ntips)
    }
  }

  return (Fj[1,j]) # return loglikelihood
}

genz1d <- function(data, tree, Lambda=1, Sigma=1, mu=0)
{
  # Calculate loglikelihood with Genz-Bretz algorithm.

  if (length(Lambda) != 1 || 
      length(Sigma)  != 1 || 
      length(mu)     != 1 || 
      length(data)   != length(tree$tip.label) )
    stop ("Dimension mismatch")

  tree   = ape::reorder.phylo(tree, "postorder")
  lower  = log(c(data))
  upper  = -log(1-c(data))
  W      = ape::vcv(tree)*Lambda + diag(length(data))*Sigma
  m      = rep(mu,length(data))

  out    = mvtnorm::pmvnorm(lower, upper, mean=m, sigma=W)
  attr(out, "error") = log(out + c(-attr(out, "error"),attr(out, "error")))
  log(out)
}
