#' @export png.ZeroReplace.simple
png.ZeroReplace.simple <- function(x, delta=0.01){
  x[abs(x)<1e-12] <- delta
  x/sum(x)
}

#' @export png.ZeroReplace.additive
png.ZeroReplace.additive <- function(x, delta=1e-6){
  x=ifelse( abs(x)<1e-12, 0, x)

  D=length(x);  Z=sum(abs(x==0))

  x[x!=0] = x[x!=0] - delta*(Z+1)*Z / D^2
  x[x==0] = delta*(Z+1)*(D-Z) / D^2

  return(x)
}

#' @export png.ZeroReplace.multiplicative
png.ZeroReplace.multiplicative <- function(x, delta){
  x=ifelse( abs(x)<1e-12, 0, x)
  c <- sum(x)
  delta_k <- rep(delta, sum(x==0))

  x[x!=0] <- (1-sum(delta_k)/c)*x[x!=0]
  x[x==0] <- (delta_k)

  return(x)
}


#' @export png.tpca
png.tpca <- function(X, nrank=2){
  d0 <- sum(svd(X)$d>1e-10)
  if( nrank > (d0-1) ) stop("nrank should be smaller than rank(X)-1")

  n=nrow(X); p=ncol(X); mu=colMeans(X)
  uhat <- prcomp(X)$x[,1:nrank,drop=F]
  vhat <- prcomp(X)$rotation[,1:nrank,drop=F]
  xhat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat, vhat)

  return( list(mu=mu, uhat=uhat, vhat=vhat, xhat=xhat, X=X) )
}
