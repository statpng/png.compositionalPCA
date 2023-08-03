#' @export V_update
V_update <- function(X, Uhat, Vhat, Uk, kappa=1e-8){
  if(FALSE){
    Uhat=Unew[,1:(r-1)]
    Vhat=Vhat
    Uk=Unew[,r]
    kappa=kappa
  }
  if(FALSE){
    set.seed(12)
    nrank=3
    X <- sim.simplex(50,500,5,eta=0.1)$X2
    png.gppca_qp(X, nrank, kappa=1e-4, gamma=0.2)
    
    biplot(prcomp(X))
    (colMeans(X)+prcomp(X)$rot[,1])[1:10,1:2] %>% lines(col="red", lwd=2)
    prcomp(X)$rot %>% head
    
      
    Uhat = prcomp(X)$x[,1]
    Uk = prcomp(X)$x[,2]
    Vhat = prcomp(X)$rotation[,1]
    # Vhat = prcomp( t(apply(X, 1, png.ZeroReplace.additive)) )$rotation[,2]
    # Vhat = qr.Q(qr(cbind(1,rnorm(ncol(X)))))[,2]
    
    V_update(X, Uhat, Vhat, Uk, kappa=1e-8)
    V_update2(X, Uhat, Vhat, Uk, kappa=1e-8)
    #
      
    
    Vhat = qr.Q(qr(cbind(1,prcomp(X)$rotation[,1]+rnorm(p)*0.01 )))[,2]
    Vhat = qr.Q(qr(cbind(1,rnorm(ncol(X)))))[,2]
      
    Uk = prcomp(X)$x[,2]
    kappa=1e-8
    
    fit <- UV_update(X, Vhat)
    fit$crit.path %>% plot(type="l")
    UV_update(X, Vhat)$vhat %>% head
    UV_update2(X, Vhat)$vhat %>% head
  }
  if(FALSE){
    set.seed(1)
    n=50; p=4; r=2; mu=rep(1,p)/p
    U=matrix(0,n,r)
    for( i in 1:n ){
      U[i,]=rnorm(r,0,1)
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
    Uhat=U[,1]; Vhat=V[,1]
    Uk=U[,2]
  }
  
  
  
  n=nrow(X); p=ncol(X); r=NCOL(Uhat)
  mu=colMeans(X);
  
  C <- tcrossprod(rep(1,n),mu) + Uhat %*% t(Vhat)
  
  Vnew <- rep(0,p)
  Dmat <- diag(rep(1,p))
  dvec <- t(X-C) %*% Uk/sum(Uk^2)
  
  positiveindex = which(Uk>0); negativeindex = which(Uk<0)
  n.pos <- length(positiveindex)
  n.neg <- length(negativeindex)
  if( n.pos>0 && n.neg>0 ){
    lbmat=-diag(1/Uk[positiveindex],n.pos,n.pos) %*% C[positiveindex,]
    lb=apply(lbmat,2,max)
    
    ubmat <- -diag(1/Uk[negativeindex],n.neg,n.neg) %*% C[negativeindex,]
    ub<-apply(ubmat,2,min)
    
    bvec <- c(rep(0,r+1),(lb-kappa),-(ub+kappa))
    
    Amat <- cbind(rep(1,p), Vhat, diag(rep(1,p)), -diag(rep(1,p)))
    
    Vnew <- solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=r+1, factorized=FALSE )$solution
  }
  
  if (n.pos>0 && n.neg==0){
    lbmat<- -diag(1/Uk[positiveindex],n.pos,n.pos) %*% C[positiveindex,] 
    
    lb<-apply(lbmat,2,max)
    
    bvec<-c(rep(0,r+1),(lb-kappa))
    
    Amat<-cbind(rep(1,p), Vhat, diag(rep(1,p)))
    
    Vnew<-solve.QP(Dmat, dvec, Amat, bvec, meq=r+1, factorized=FALSE )$solution
  }
  
  if (n.pos==0 && n.neg>0){
    ubmat<- -diag(1/Uk[negativeindex],n.neg,n.neg) %*% C[negativeindex,] 
    
    ub<-apply(ubmat,2,min)
    
    bvec<-c(rep(0,r+1),-(ub+kappa))
    
    Amat<-cbind(rep(1,p), Vhat, -diag(rep(1,p)))
    
    Vnew<-solve.QP(Dmat, dvec, Amat, bvec, meq=r+1,factorized=FALSE )$solution
  }
  
  # if (n.pos==0 && n.neg==0){
  #   # dvec <- t(X-Xhat) %*% Uk
  #   # bvec<-c(rep(0,r+1))
  #   # Amat<-cbind(rep(1,p), Vhat)
  #   # 
  #   # Vnew<-solve.QP(Dmat, dvec, Amat, bvec, meq=r+1,factorized=FALSE )$solution
  #   Vnew <- prcomp(X)$rotation[,1]
  # }
  
  
  if(norm(Vnew,"2")>1e-10){
    Vnew <- Vnew/norm(Vnew,"2")
  }
  
  
  return(Vnew)
}
