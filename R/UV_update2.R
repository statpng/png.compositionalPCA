#' @export V_update2
V_update2 <- function(X, Uhat, Vhat, Uk, kappa=1e-8){
  if(FALSE){
    X=X; Uhat=Unew[,1:(r-1)]; Vhat=Vhat; Uk=Unew[,r]; kappa=kappa
  }
  
  n=nrow(X); p=ncol(X); r=NCOL(Uhat)
  mu=colMeans(X)
  
  C <- tcrossprod(rep(1,n),mu) + tcrossprod(Uhat, Vhat)
  
  Vnew <- rep(0,p)
  Dmat <- diag(rep(1,p))
  dvec <- t(X-C) %*% Uk/sum(Uk^2)
  
  solve_V <- function(Vhat, lbmat=NULL, ubmat=NULL, kappa=0) {
    p=NROW(Vhat);  r=NCOL(Vhat)
    lb <- if(!is.null(lbmat)) apply(lbmat,2,max) else NULL
    ub <- if(!is.null(ubmat)) apply(ubmat,2,min) else NULL
    
    bvec <- c(rep(0,r+1), lb-kappa, -(ub+kappa))
    Amat <- cbind(rep(1,p), Vhat, if (!is.null(lbmat)) diag(rep(1,p)) else NULL, if (!is.null(ubmat)) -diag(rep(1,p)) else NULL)
    Vnew <- try( solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=r+1, factorized=FALSE)$solution, silent=TRUE )
    if(inherits(Vnew, "try-error")){
      return( "No solution" )
    } else {
      return( Vnew/norm(Vnew,"2") )
    }
    # bvec <- c(rep(0,1), rep(-1e-12,2*r), lb-kappa, -(ub+kappa))
    # Amat <- cbind(rep(1,p), Vhat, -Vhat, diag(rep(1,p)), -diag(rep(1,p)))
    # 
    # solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=1, factorized=FALSE)$solution
  }
  
  positiveindex = which(Uk>0); negativeindex = which(Uk<0)
  n.pos=length(positiveindex); n.neg=length(negativeindex)
  
  if(n.pos > 0 && n.neg > 0) {
    Vnew <- solve_V(Vhat=Vhat,
                    lbmat=-diag(1/Uk[positiveindex],n.pos,n.pos) %*% C[positiveindex,],
                    ubmat=-diag(1/Uk[negativeindex],n.neg,n.neg) %*% C[negativeindex,],
                    kappa=kappa)
  }
  
  if (n.pos>0 && n.neg==0) {
    Vnew <- solve_V(Vhat=Vhat,
                    lbmat=-diag(1/Uk[positiveindex],n.pos,n.pos) %*% C[positiveindex,], 
                    kappa=kappa)
  }
  
  if (n.pos==0 && n.neg>0) {
    Vnew <- solve_V(Vhat=Vhat,
                    ubmat=-diag(1/Uk[negativeindex],n.neg,n.neg) %*% C[negativeindex,], 
                    kappa=kappa)
  }
  
  return(Vnew)
}


#' @export Solve_U_SP
Solve_U_SP <- function(x, mu, v, gamma=0){
  if(FALSE){
    x=X[i,]; mu=C[i,]; v=Vold; gamma=gamma/(it^(1/2))
  }
  
  t=sum((x-mu)*v)
  
  indexminus=which(v<0);  indexplus=which(v>0)
  lminus=length(indexminus);  lplus=length(indexplus)
  
  if (lminus>0 & lplus>0){
    m <- max( -(mu/v)[indexplus]);  M=min( -(mu/v)[indexminus])
    U <- min(max(m,t),M)
  }
  
  if (lminus>0 & lplus==0){
    m <- max( -(mu/v)[indexminus])
    U <- max(m,t)
  }
  
  if (lminus==0 & lplus>0){
    M <- min( -(mu/v)[indexplus])
    U <- min(t,M)
  }
  
  if (lminus==0 & lplus==0){
    U <- t
  }
  
  return(U * (1-gamma))
}

#' @export Solve_U_GP
Solve_U_GP <- function(x, mu, V, gamma=0){
  if(FALSE){
    x=X[i,]; mu=C[i,]; V=cbind(Vhat, Vk_old); gamma=gamma/(it^(1/2))
  }
  
  x=as.numeric(x)
  mu=as.numeric(mu)
  
  U <- solve.QP(Dmat=crossprod(V), dvec=t(V) %*% (x - mu), Amat=t(V), bvec=-mu - 1e-8)$solution
  # U <- solve.QP(Dmat=crossprod(V), dvec=t(V) %*% (x - mu), Amat=t(V), bvec=-mu)$solution
  return(U * (1-gamma))
}

#' @export UV_update2
UV_update2 <- function(X, Uhat, Vhat, maxit=500, eps=1e-6, kappa=1e-4, gamma=0, phi=0.01, V.init=c("PC","random"), verbose=TRUE){
  if(FALSE){
    Uhat=U_total; Vhat=V_total
  }
  
  
  n=nrow(X); p=ncol(X); r=NCOL(Vhat) + 1
  mu=colMeans(X)
  
  C <- tcrossprod(rep(1,n), mu)
  # C <- tcrossprod(rep(1,n), mu) + X %*% tcrossprod(Vhat)
  
  
  # Initial V0
  V0.PC <- prcomp(X-(tcrossprod(rep(1,n), mu) + X %*% tcrossprod(Vhat)))$rotation[,1]
  if(V.init == "PC"){
    V0 <- V0.PC
  } else {
    V0 <- (V0.PC + qr.Q(qr(cbind(1,rnorm(ncol(X)))))[,2]*phi) %>% {./norm(.,"2")}
  }
  
  # png.angle(V0.PC, V0)[[1]] %>% print
  
  crit.path <- est.path <- NULL
  for( it in 1:maxit ){
    if(verbose) print(paste0("it: ", it))
    
    if( it == 1 ){
      Vk_old <- V0
    } else {
      Vk_old <- Vk_new
    }
    
    
    
    # U-update
    Unew <- t(sapply(1:n, function(i) Solve_U_GP(x=X[i,], mu=C[i,], V=cbind(Vhat, Vk_old), gamma=gamma/(it^(1/2)))))
    # V-update
    Vk_new <- V_update2(X, Uhat=Unew[,1:(r-1)], Vhat=Vhat, Uk=Unew[,r], kappa=kappa)
    if( Vk_new[1] == "No solution" ){
      Vk_new <- Vk_old
      NoSolution <- TRUE
    } else {
      NoSolution <- FALSE
    }
    
    crit <- min( sqrt(mean((Vk_old-Vk_new)^2)), sqrt(mean((Vk_old+Vk_new)^2)) )
    
    crit.path[it] <- crit
    # est.path[[it]] <- list(uhat=Unew, vhat=cbind(Vhat,Vk_new))
    
    if( crit < eps ) break
  }
  
  if(NoSolution) crit.path[it] <- 0.1
  
  Vnew <- cbind(Vhat, Vk_new)
  
  # U-update for Vnew
  Unew <- t(sapply(1:n, function(i) Solve_U_GP(x=X[i,], mu=C[i,], V=Vnew, gamma=0)))
  # Unew <- t(apply(X, 1, function(x) Solve_U_GP(x=x, mu=mu, V=Vnew, gamma=0)))
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(Unew, Vnew)
  
  loss <- sqrt(mean((X - xhat)^2))
  
  # return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, est.path=est.path) )
  return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, loss=loss) )
}
