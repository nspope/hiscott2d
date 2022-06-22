# -------- algorithms for 2 traits -------#

hiscott2d_vanilla <- function(data, tree, Lambda=diag(2), Sigma=diag(2), mu=rep(0,2), level=3, weights=NULL, tol=0)
{
  # Calculate loglikelihood with multivariate Hiscott algorithm.

  if (nrow(Lambda) != 2 || ncol(Lambda) != 2 ||
      nrow(Sigma)  != 2 || ncol(Sigma)  != 2 ||
      length(mu)   != 2 || nrow(data)   != 2 ||
      ncol(data)   != length(tree$tip.label) )
  {
    stop("Dimension mismatch")
  }

  ww <- new(quadrature2d, level, tol)
  tree <- ape::reorder.phylo(tree, "postorder")
  nnodes <- tree$Nnode
  ntips <- length(tree$tip.label)
  npoint <- level*2 + 1
  len <- tree$edge.length
  root <- tree$edge[nrow(tree$edge),1]

  if (!is.null(weights) && class(weights) == "list")
  {
    if (length(weights) != tree$Nnode-1 || !all(sapply(weights, function(x) all(dim(x)==npoint^2))))
    {
      stop("Weights of wrong dimension")
    }
  }

  if (!is.null(weights) && class(weights) != "list")
  {
    Fj <- rep(list(), nnodes) # weights
  }
  else
  {
    Fj <- matrix(0, npoint^2, nnodes) # log'd function evaluations
    Fjs <- matrix(1, npoint^2, nnodes) # sign of function evaluations
  }

  # bounds
  hei <- ape::node.depth.edgelength(tree)
  Lb <- matrix(0, 2, nnodes)
  Sc <- matrix(0, 2, nnodes)
  err <- (npoint-1)^-4 * ntips/(2*ntips-1)
  sdv <- sqrt(diag(Lambda))
  for (node in 1:nnodes)
  {
    Lb[,node] <- qnorm(err, mu, sqrt(hei[node+ntips])*sdv)
    Sc[,node] <- (qnorm(1-err, mu, sqrt(hei[node+ntips])*sdv) - Lb[,node])/(npoint-1)
  }

  # quadrature
  for (edge in 1:nrow(tree$edge))
  {
    i <- tree$edge[edge,2]
    j <- tree$edge[edge,1]-ntips

    if (i>ntips) 
    # internal node
    {
      if (!is.null(weights) && class(weights) == "list") 
      {
        # weights supplied
        mm <- max(Fj[,i-ntips])
        Fjt <- as.vector(weights[[i-ntips-1]] %*% (exp(Fj[,i-ntips] - mm) * Fjs[,i-ntips]))
        Fjs[,j] <- Fjs[,j] * sign(Fjt)
        Fj[,j] <- Fj[,j] + log(abs(Fjt)) + mm
      # b/c of sparse matrix format, do this in R
      #    ww$node3(Fjs[,j], 
      #             Fj[,i-ntips], 
      #             Fjs[,i-ntips],
      #             weights[[i-ntips]],
      #             j==root-ntips)
      }
      else if (!is.null(weights) && class(weights) != "list") 
      {
        # weights calculated with no quadrature
        Fj[[i-ntips-1]] <-
          Matrix::Matrix(
          ww$node1(Lb[,j],
                   Lb[,i-ntips],
                   Sc[,j],
                   Sc[,i-ntips],
                   diag(Lambda),
                   Lambda[1,2],
                   len[edge],
                   j==root-ntips))
      }
      else 
      {
        # weights calculated with quadrature
        Fj[,j] <- Fj[,j] + 
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
    } 
    else 
    # tip node
    {
      if (!is.null(weights) && class(weights) != "list") {}
      else 
      {
        Fj[,j] <- Fj[,j] + 
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
  }

  if (!is.null(weights) && class(weights) != "list")
  {
    return (Fj) # return weights
  }
  else
  {
    return (Fj[1,j]) # return loglikelihood
  }
}

hiscott2d <- function(data, tree, Lambda=diag(2), Sigma=diag(2), mu=rep(0,2), level=3)
{
  # Calculate loglikelihood with multivariate Hiscott algorithm.
  # Uses rotated quadrature grid. Allows for mixed-type traits.

  if (nrow(Lambda) != 2 || ncol(Lambda) != 2 ||
      nrow(Sigma)  != 2 || ncol(Sigma)  != 2 ||
      length(mu)   != 2 || nrow(data)   != 2 ||
      ncol(data)   != length(tree$tip.label) )
  {
    stop("Dimension mismatch")
  }

  # Determine type of trait
  type <- apply(data, 1, function(x) all(x >= 0) & all(abs(x) < .Machine$double.eps | abs(1-x) < .Machine$double.eps))
  cat("Trait types:", c("continuous", "binary")[as.numeric(type)+1], "\n", sep=" ")
  if (all(!type))
  {
    stop("Both traits continuous")
  }
  if (type[2] & !type[1])
  {
    stop("Not working with mixed continuous/binary")
    Lambda <- Lambda[c(2,1),c(2,1)]
    Sigma <- Sigma[c(2,1),c(2,1)]
    mu <- mu[c(2,1)]
    data <- data[c(2,1),]
  }

  tree <- ape::reorder.phylo(tree, "postorder")
  nnodes <- tree$Nnode
  ntips <- length(tree$tip.label)
  npoint <- level*2 + 1
  len <- tree$edge.length
  root <- tree$edge[nrow(tree$edge),1]
  ww <- new(quadrature2drot, level, ntips)

  Fj <- matrix(0, npoint^2, nnodes) # log'd function evaluations
  Fjs <- matrix(1, npoint^2, nnodes) # sign of function evaluations

  # scaling etc.
  hei <- ape::node.depth.edgelength(tree)
  Lambda <- eigen(Lambda)
  Lambda <- Lambda$vectors %*% diag(sqrt(Lambda$values))
  Sigma <- eigen(Sigma)
  Sigma <- Sigma$vectors %*% diag(sqrt(Sigma$values))

  # quadrature
  for (edge in 1:nrow(tree$edge))
  {
    i <- tree$edge[edge,2]
    j <- tree$edge[edge,1]

    if (i>ntips) 
    {
      # internal node
      Fj[,j-ntips] <- Fj[,j-ntips] + 
        ww$node(Fjs[,j-ntips], 
                Fj[,i-ntips], 
                Fjs[,i-ntips],
                hei[j],
                hei[i],
                mu,
                j==root)
    }
    else 
    {
      # tip node
      if (all(type))
      {
        Fj[,j-ntips] <- Fj[,j-ntips] + 
          ww$tip(hei[j],
                 hei[i],
                 mu,
                 Lambda,
                 Sigma,
                 data[,i],
                 j==root)
      }
      else
      {
        Fj[,j-ntips] <- Fj[,j-ntips] + 
          ww$tip2(hei[j],
                  hei[i],
                  mu,
                  Lambda,
                  Sigma,
                  as.integer(data[1,i]),
                  data[2,i],
                  j==root)
      }
    }
  }

  return (Fj[1,j-ntips]) # return loglikelihood
}

genz2d <- function(data, tree, Lambda=diag(2), Sigma=diag(2), mu=rep(0,2))
{
  # Calculate loglikelihood with Genz-Bretz algorithm.

  if (nrow(Lambda) != 2 || ncol(Lambda) != 2 ||
      nrow(Sigma)  != 2 || ncol(Sigma)  != 2 ||
      length(mu)   != 2 || nrow(data)   != 2 ||
      ncol(data)   != length(tree$tip.label) )
  {
    stop("Dimension mismatch")
  }

  type <- apply(data, 1, function(x) all(x >= 0) & all(abs(x) < .Machine$double.eps | abs(1-x) < .Machine$double.eps))
  cat("Trait types:", c("continuous", "binary")[as.numeric(type)+1], "\n", sep=" ")
  if (all(!type))
  {
    stop("Both traits continuous")
  }
  if (type[2] & !type[1])
  {
    Lambda <- Lambda[c(2,1),c(2,1)]
    Sigma <- Sigma[c(2,1),c(2,1)]
    mu <- mu[c(2,1)]
    data <- data[c(2,1),]
  }

  tree <- ape::reorder.phylo(tree, "postorder")
  
  if (all(type))
  {
    lower <- log(c(data))
    upper <- -log(1-c(data))
    W <- kronecker(ape::vcv(tree), Lambda) + kronecker(diag(ncol(data)), Sigma)
    m <- c(kronecker(rep(1,ncol(data)), mu))
    out <- mvtnorm::pmvnorm(lower, upper, mean=m, sigma=W)
  }
  else
  {
    k <- seq(1, ncol(data)*2, 2)
    lower <- log(c(data[1,]))
    upper <- -log(1-c(data[1,]))
    W <- kronecker(ape::vcv(tree), Lambda) + kronecker(diag(ncol(data)), Sigma)
    m <- c(kronecker(rep(1,ncol(data)), mu))
    Wc <- W[k,k] - W[k,-k] %*% solve(W[-k,-k]) %*% W[-k,k]
    mc <- m[k] + c(W[k,-k] %*% solve(W[-k,-k]) %*% (c(data[2,]) - m[-k]))
    out <- mvtnorm::pmvnorm(lower, upper, mean=mc, sigma=Wc) * mvtnorm::dmvnorm(c(data[2,]), m[-k], sigma=W[-k,-k])
  }

  attr(out, "error") <- log(out + c(-attr(out, "error"),attr(out, "error")))
  log(out)
}

# ------------ algorithms for 1 trait --------------#

hiscott1d <- function(data, tree, Lambda=1, Sigma=1, mu=0, level=3)
{
  # Calculate loglikelihood with univariate Hiscott algorithm.

  if (length(Lambda) != 1 || 
      length(Sigma)  != 1 || 
      length(mu)     != 1 || 
      length(data)   != length(tree$tip.label) )
  {
    stop ("Dimension mismatch")
  }

  ww <- new(quadrature1d, level)
  tree <- ape::reorder.phylo(tree, "postorder")
  nnodes <- tree$Nnode
  ntips <- length(tree$tip.label)
  npoint <- level*2 + 1
  len <- tree$edge.length
  root <- tree$edge[nrow(tree$edge),1]

  Fj <- matrix(0, npoint, nnodes) # log'd function evaluations
  Fjs <- matrix(1, npoint, nnodes) # sign of function evaluations

  # bounds
  hei <- ape::node.depth.edgelength(tree)
  Lb <- matrix(0, 1, nnodes)
  Sc <- matrix(0, 1, nnodes)
  err <- (npoint-1)^-4 * ntips/(2*ntips-1)
  sdv <- sqrt(Lambda)
  for (node in 1:nnodes)
  {
    Lb[node] <- qnorm(err, mu, sqrt(hei[node+ntips])*sdv)
    Sc[node] <- (qnorm(1-err, mu, sqrt(hei[node+ntips])*sdv) - Lb[node])/(npoint-1)
  }

  # quadrature
  for (edge in 1:nrow(tree$edge))
  {
    i <- tree$edge[edge,2]
    j <- tree$edge[edge,1]-ntips

    if (i>ntips) 
    # internal node
    {
      Fj[,j] <- Fj[,j] + 
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
      Fj[,j] <- Fj[,j] + 
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

genz1d <- function(data, tree, Lambda=1, Sigma=1, mu=0, algorithm=mvtnorm::GenzBretz())
{
  # Calculate loglikelihood with Genz-Bretz algorithm.

  if (length(Lambda) != 1 || 
      length(Sigma) != 1 || 
      length(mu) != 1 || 
      length(data) != length(tree$tip.label) )
  {
    stop("Dimension mismatch")
  }

  tree <- ape::reorder.phylo(tree, "postorder")
  lower <- log(c(data))
  upper <- -log(1-c(data))
  W <- ape::vcv(tree)*Lambda + diag(length(data))*Sigma
  m <- rep(mu,length(data))

  out <- mvtnorm::pmvnorm(lower, upper, mean=m, sigma=W, algorithm=algorithm)
  attr(out, "error") <- log(out + c(-attr(out, "error"),attr(out, "error")))
  log(out)
}
