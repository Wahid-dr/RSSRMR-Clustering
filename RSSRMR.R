
#' Robust Self-Representation Sparse Reconstruction and Manifold Regularization
#' @param x A featureS matrix of dimension m and p
#' @param Wt a weight matrix
#' @param  L A laplacian matrix
#' @param alpha regularization parameter
#' @param beta menofold parameter
#' @param epsilon convergence threshold
#' @param maxites tolerance limit
#' @return a FS vector
#' @export
#' @examples
#' set.seed(412)
#'
#' cluster1 <- matrix(
#'   rnorm(12 * 5, mean = 2, sd = 0.5),
#'   nrow = 12,
#'   ncol = 5
#' )
#'
#' cluster2 <- matrix(
#'   rnorm(13 * 5, mean = 7, sd = 0.5),
#'   nrow = 13,
#'   ncol = 5
#' )
#'
#' X <- rbind(cluster1, cluster2)
#'
#' ws <- runif(25)
#' wd <- diag(ws)
#'
#' ls <- runif(25)
#' lp <- diag(ls)
#'
#' RSSRMR(
#'   x = X,
#'   Wt = wd,
#'   L = lp,
#'   alpha = 1.0,
#'   beta = 1.0,
#'   epsilon = 0.001,
#'   maxites = 5
#' )
##### proposed function;
RSSRMR <- function(x,Wt,L,alpha=1.0,beta=1.0,epsilon=0.001,maxites=50){
  n <- dim(x)[1];
  d <- dim(x)[2];
  Wt <- diag(n);
  #### Laplacian matrix;
  L <- matrix(1,nrow=n,ncol=n);
  ######### Initialization
  opt.A <- matrix(1,nrow=d,ncol=d);
  opt.G <- diag(d);
  delta<- 1;
  step <- 0;
  while(delta >= epsilon & step <= maxites){
    step<- step + 1;
    old.A <- opt.A;
    old.G <- opt.G;
    Q <- (t(x)%*%Wt%*%x + alpha*t(x)%*%L%*%x + beta*old.G);
    inv.Q <- solve(Q);
    Q1 <- (t(x)%*%Wt%*%x);
    opt.A <- inv.Q%*%Q1;
    aii <- c(diag(opt.A));
    opt.G <- diag(1/(2*aii));
    delta <- max(abs(old.A - opt.A));
    opt.A = opt.A;
  }
  print("Algorithm Converged")
  values <- list("Optimum.A"= opt.A,"Optimum.G"=opt.G)
  return(values)
}
