#' @export UV_update
UV_update <- function(X, Vhat, maxit=500, eps=1e-6, kappa=1e-4, gamma=0){
  if(FALSE){
    # UV_update(X, V_total, eps=eps)
    Vhat=V_total;  eps=1e-6;  maxit=500
  }
  
  if(FALSE){
    set.seed(1)
    n=50; p=4; r=2; mu=rep(1,p)/p
    U=matrix(0,n,r)
    for( i in 1:n ){
      U[i,]=runif(r,-0.1,0.1)
    }
    v1=runif(p,0,1); v1=v1-sum(v1*1/p)
    v2=runif(p,0,1); v2=v2-sum(v2*1/p); v2=v2-sum(v1*v2)/sum(v1*v1)*v1
    v1=v1/norm(v1,"2"); v2=v2/norm(v2,"2")
    V=cbind(v1,v2)
    E=matrix(0,n,p)
    for( i in 1:n ){
      E[i,]=runif(p,-0.1/p,0.1/p); E[i,]=E[i,]-sum(E[i,]/p)
    }
    
    X=tcrossprod(rep(1,n),mu)+tcrossprod(U%*%diag(c(4,1),r,r),V)+E
    X=t(apply(X,1,png.proj2simplex))
    Vhat=V[,1]
    
    # png.quaternary(X)
    # png.quaternary3d(X)
    
    UV_update(X, Vhat)
  }
  
  
  
  
  
  n=nrow(X); p=ncol(X); r=NCOL(Vhat)+1
  mu=colMeans(X);
  C <- tcrossprod(rep(1,n),mu) + X %*% tcrossprod(Vhat)
  
  # Vhat + Initial V0
  V0=prcomp( X - C )$rotation[,1]
  
  
  # U-update
  Unew <- matrix(0,n,r)
  for( i in 1:n ){
    x=as.vector(X[i,]); V=cbind(Vhat,V0)
    
    Unew[i,] <- solve.QP(Dmat=crossprod(V),
                         dvec=t(V) %*% (x-C[i,]),
                         Amat=t(V),
                         bvec=-C[i,]-1e-8 )$solution
  }
  
  Unew <- Unew * (1-gamma)
  
  # V-update
  Vk_new <- V_update(X, Uhat=Unew[,1:(r-1)], Vhat=Vhat, Uk=Unew[,r], kappa=kappa)
  
  
  est.path <- crit.path <- NULL
  for( it in 1:maxit ){
    if(verbose) print(paste0("it: ", it))
    
    Vk_old <- Vk_new
    U_old <- Unew
    
    # U-update
    Unew <- matrix(0,n,r)
    for( i in 1:n ){
      x=as.vector(X[i,]);  V=cbind(Vhat,Vk_new)
      
      Unew[i,] <- solve.QP(Dmat=diag(rep(1,r)), 
                           dvec=t(V) %*% (x-C[i,]), 
                           Amat=t(V), 
                           bvec=-C[i,]-1e-8 )$solution
    }
    Unew <- Unew * (1-gamma/it)
    
    # V-update
    Vk_new <- V_update(X, Uhat=Unew[,1:(r-1)], Vhat=Vhat, Uk=Unew[,r], kappa=kappa)
    
    est.path[[it]] <- list(uhat=Unew, vhat=cbind(Vhat,Vk_new))
    crit.path[it] <- min( sqrt(mean((Vk_old-Vk_new)^2)), sqrt(mean((Vk_old+Vk_new)^2)) )
    
    if( crit.path[it] < eps ){
      break;
    }
  }
  
  Vnew <- cbind(Vhat, Vk_new)
  
  Unew <- matrix(0,n,r);
  for( i in 1:n ){
    Unew[i,] <- solve.QP(Dmat=diag(rep(1,r)),
                         dvec=t(Vnew) %*% (as.vector(X[i,])-mu),
                         Amat=t(Vnew),
                         bvec=-mu-kappa )$solution
  }
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(Unew,Vnew)
  
  return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, est.path=est.path, crit.path=crit.path) )
  
}
















#' @export png.rank1
png.rank1 <- function(X, maxit=500, eps=1e-6, kappa=1e-4, gamma=0){
  
  require(quadprog)
  
  n=nrow(X); p=ncol(X);
  mu=colMeans(X); 
  C=tcrossprod(rep(1,n), mu)
  
  # Initial V0
  V0=prcomp(X-tcrossprod(rep(1,n),mu))$rotation[,1]
  Vnew=V0
  
  est.path <- crit.path <- NULL
  for( it in 1:maxit ){
    Vold<-Vnew
    
    Unew<-rep(0,n)
    for(i in 1:n){
      # Unew[i] <- Solve_U_SP(X[i,], C[i,],Vold)
      Unew[i]<-onedimconvexprojection(C[i,],X[i,],Vold)
    }
    Unew <- Unew * (1-gamma/sqrt(it))
    
    Vnew <- V_update2(X, Uhat=matrix(0,n,1), Vhat=matrix(0,p,1), Uk=as.matrix(Unew), kappa=kappa)
    if(norm(Vnew,"2")!=0) Vnew <- Vnew/sqrt(sum(Vnew^2))
    
    est.path[[it]] <- list(uhat=Unew, vhat=Vnew)
    crit.path[it] <- min( sqrt(mean((Vnew-Vold)^2)), sqrt(mean((Vnew+Vold)^2)) )
    
    if( crit.path[it] < eps ){
      break;
    }
    
  }
  
  Unew<-rep(0,n)
  for(i in 1:n){
    Unew[i]<-onedimconvexprojection(C[i,],X[i,],Vnew)
  }
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(Unew,Vnew)
  
  return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, est.path=est.path, crit.path=crit.path) )
  
}


















#' @export png.rank12
png.rank12 <- function(X, maxit=500, eps=1e-6, kappa=1e-4, gamma=0, phi=0.01, V.init=c("PC","random"), verbose=TRUE){
  
  require(quadprog)
  
  n=nrow(X); p=ncol(X);
  mu=colMeans(X); 
  C=tcrossprod(rep(1,n), mu)
  
  # Initial V0
  V0.PC <- prcomp(X-C)$rotation[,1]
  if(V.init == "PC"){
    V0 <- V0.PC
  } else {
    V0 <- (V0.PC + qr.Q(qr(cbind(1,rnorm(ncol(X)))))[,2]*phi) %>% {./norm(.,"2")}
  }
  
  # png.angle(V0.PC, V0)[[1]] %>% print
  
  
  est.path <- crit.path <- NULL
  for( it in 1:maxit ){
    if(verbose) print(paste0("it: ", it))
    
    if( it == 1 ){
      Vold <- V0
    } else {
      Vold <- Vnew
    }
    
    
    Unew <- sapply(1:n, function(i) Solve_U_SP(x=X[i,], mu=C[i,], v=Vold, gamma=gamma/(it^(1/2))))
    
    Vnew <- V_update2(X, Uhat=matrix(0,n,1), Vhat=matrix(0,p,1), Uk=as.matrix(Unew), kappa=kappa)
    
    crit <- min( sqrt(mean((Vnew-Vold)^2)), sqrt(mean((Vnew+Vold)^2)) )
    
    crit.path[it] <- crit
    # est.path[[it]] <- list(uhat=Unew, vhat=Vnew)
    
    if( crit < eps ) break
  }
  
  
  Unew <- sapply(1:n, function(i) Solve_U_SP(x=X[i,], mu=C[i,], v=Vnew, gamma=0))
  
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(Unew,Vnew)
  
  loss <- sqrt(mean((X - xhat)^2))
  
  # return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, est.path=est.path) )
  return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, loss=loss) )
  
}
