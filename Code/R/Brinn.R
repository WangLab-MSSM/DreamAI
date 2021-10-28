#' title: IRNN algorithm for matrix completion
#' Arguments:
#' fun_sg: {scad_sg, mcp_sg}, the supergradient of the selected nonconvex surrogate function of L0 norm
#' y: a vector with all entry values
#' IDX: a vector with all entry positions
#' m: the number of rows present in the incomplete dataset (proteins)
#' n: the number of columns present in the incomplete dataset (samples)
#' gamma: 'fun_sg''s parameter
#' x_initial: the initial matrix
#'
#' return: an imputed matrix
#'
mc.IRNN <- function(fun_sg, y, IDX, m, n, gamma, x_initial)
{
  err <- 1e-6
  normfac <- 1
  insweep <- 200
  tol <- 1e-5
  decfac <- 0.9

  mu <- 1.1 * normfac
  x <- as.vector(as.matrix(x_initial))
  x.matrix <- x_initial
  lambdaInit <- decfac * max(abs(Mfun(m*n, IDX, y, 2)))
  lambda <- lambdaInit
  f_current <- sqrt(sum((y-Mfun(m*n, IDX, x, 1))^2)) + lambda * sum(abs(x))

  while (lambda > lambdaInit * tol)
  {
    for (ins in 1:insweep)
    {
      f_previous = f_current
      x = x + Mfun(m*n, IDX, y - Mfun(m*n, IDX, x, 1), 2) / mu
      x.matrix = matrix(x, nrow=m)
      ans.svd <- svd(x.matrix)
      sigma = ans.svd$d
      U = ans.svd$u
      V = ans.svd$v
      w = fun_sg(sigma, gamma, lambda)
      sigma = sigma - w/mu
      svp = length(which(sigma > 0))
      if (svp > 1) {
        sigma = sigma[1:svp]
        x.matrix = U[,1:svp] %*% diag(sigma) %*% t(V[,1:svp])
      } else if (svp == 1) {
        x.matrix = tcrossprod(U[,1], V[,1]) * sigma[1]
      } else { break }
      x = as.vector(x.matrix)
      f_current = sqrt(sum((y-Mfun(m*n, IDX, x, 1))^2)) + lambda * sum(abs(x))
      if (abs(f_current - f_previous)/abs(f_current + f_previous) < tol) {break}
    }
    if (sqrt(sum((y-Mfun(m*n, IDX, x, 1))^2)) < err) {break}
    lambda = decfac * lambda
  }
  return(x.matrix)
}

# the supergradients of popular nonconvex surrogate functions of L0-norm
scad_sg <- function(x, gamma, lambda)
{
  x <- abs(x)
  y <- rep(0, length(x))
  ind <- which(x <= lambda)
  y[ind] <- lambda
  ind <- which(x>lambda & x<=gamma*lambda)
  y[ind] <- (gamma*lambda-x[ind]) / (gamma-1)
  return (y)
}

mcp_sg <- function(x, gamma, lambda)
{
  x = abs(x)
  y = rep(0, length(x))
  ind = which(x < gamma * lambda)
  y[ind] = lambda - x[ind]/gamma
  return(y)
}

Mfun <- function(N, IDX, x, mode)
{
  if (mode == 1) {y <- x[IDX]}
  else if (mode == 2)
  {
    y <- rep(0, N)
    y[IDX] <- x
  }
  return (y)
}

#'
#' title: the selection process of parameter gamma using Cross Validation
#'
#' Arguments:
#' x: a matrix or data frame with either NAs or 0s as missings.
#'
#' return: the best performed (with lowest MSE) gamma among
#'         (1e-3, 1e-2, 1e-1, 1, 10, 20, ..., 90, 100, 120, 140, ..., 480, 500)
#'
gamma.CV <- function (x)
{
  m <- dim(x)[1]
  n <- dim(x)[2]
  x[is.na(x)] <- 0
  allmean <- mean(x[x != 0])
  x <- as.matrix(x)
  x.mis <- which(x == 0)
  x.vec <- 1:(m*n)
  x.nomis <- x.vec[is.na(match(x.vec, x.mis))]
  l <- min(length(x.mis), length(x.nomis)/2)
  x.moremis <- sample(x.nomis, l)
  x.obs.vec <- matrix(x, ncol = 1)
  x.obs.vec[x.moremis] <- 0
  x.new <- matrix(x.obs.vec, nrow = m, ncol = n)
  x.mean <- t(apply(x.new, 1, function(z) {
    if (sum(z[z!=0]) > 0) {
      z[z==0] <- mean(z[z!=0])
    }
    else {
      z[z==0] <- allmean
    }
    z}))
  IDX <- which(x.new != 0)
  y <- x.new[x.new != 0]

  glist <- c(1e-3, 1e-2, 0.1, 1, seq(10, 100, 10), seq(120, 500, 20))
  gmse <- rep(NA, length(glist))

  for (i in 1:length(glist))
  {
    XRec <- as.data.frame(mc.IRNN(scad_sg, y, IDX, m, n, glist[i], x.mean))
    gmse[i] <- CV.mse(XRec, x, x.moremis)
  }
  id <- order(gmse, decreasing = FALSE)[1]

  return (glist[id])
}

CV.mse <- function(x.imp, x.true, v.posi)
{
  ans <- mean((matrix(as.matrix(x.imp),ncol = 1)[v.posi]-matrix(as.matrix(x.true),ncol = 1)[v.posi])^2)
  return (ans)
}

#'
#' @title Imputation of Missing Protein Abundances using IRNN-SCAD Algorithm
#' @description The function impute.Birnn imputes a dataset with missing values or NA's using IRNN-SCAD Algorithm
#' @param x dataset in the form of a matrix or data frame with either NAs or 0s as missings
#' @param gamma the parameter of  the nonconvex surrogate function of L0 norm.
#'              The default \code{gamma = 50} works well for a 7000 * 80 protein-sample matrix.
#'              Can be automatically fit by setting \code{CV = TRUE}.
#' @param CV a logical value indicating whether to fit the best gamma with Cross Validation.
#'           The process can be time consuming.
#'
#' @return the imputed version of the dataset
#'
#' @examples
#' \dontrun{
#' data(datapnnl)
#' data<-datapnnl.rm.ref[1:100,1:21]
#' impute.Birnn(data=as.matrix(data), gamma = 50, CV = FALSE)
#' }
#' @export
#'
impute.Birnn <- function (x, gamma = gamma, CV = CV)
{
  x[is.na(x)] <- 0
  if (CV) {
    g = gamma.CV(x)
  }
  else {
    g = gamma
  }

  allmean <- mean(x[x != 0])
  x.mean <- t(apply(x,1,function(y){
    if (sum(y[y!=0])>0) {
      y[y==0] <- mean(y[y!=0])
    }
    else {
      y[y==0] <- allmean
    }
    y}))
  m <- dim(x)[1]
  n <- dim(x)[2]
  IDX <- which(x != 0)
  y <- x[x != 0]
  XRec <- as.data.frame(mc.IRNN(scad_sg, y, IDX, m, n, g, x.mean))
  colnames(XRec) <- colnames(x)
  rownames(XRec) <- rownames(x)
  return(XRec)
}




