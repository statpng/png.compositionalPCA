if(FALSE){
  rm(list=ls())
  devtools::load_all("./png.compositionalPCA")
}


#' @export png.projection
png.projection <- function(X, fit, method=c("ppca_qp", "gppca_qp")){
  if(FALSE){
    X=Xtrain; method=fit$method
  }
  
  n=nrow(X); p=ncol(X); r=ncol(fit$vhat)
  mu=fit$mu
  vhat=fit$vhat
  
  if( method == "ppca_qp" ){
    
    uhat <- matrix(0,n,r)
    for( k in 1:r ){
      chat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat[,1:(k-1)], vhat[,1:(k-1)])
      for( i in 1:n ){
        uhat[i,k] <- onedimconvexprojection(chat[i,], X[i,], vhat[,k])
      }
      it=fit$fit.path[[k]]$it;  gamma=fit$params$gamma
      uhat[,k] <- uhat[,k] * (1-gamma/it)
    }
    
  } else if( method == "gppca_qp" ){
    
    uhat <- matrix(0,n,r)
    for( i in 1:n ){
      uhat[i,] <- multidimconvexprojection(mu, X[i,], vhat)
    }
    it=fit$fit.path[[r]]$it;  gamma=fit$params$gamma
    uhat <- uhat * (1-gamma/it)
  }
  
  xhat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat, vhat)
  
  
  return(xhat)
}




#' @export png.gppca_qp
png.gppca_qp <- function(X, nrank=2, maxit=500, eps=1e-4, kappa=1e-6, gamma=0.5, save.est.path=FALSE){
  if(FALSE){
    nrank=2; epsilon=1e-4; maxit=100; kappa=1e-8
  }
  
  if(FALSE){
    n=150; p=4; r=2;
    seed <- 5
    data <- sim.simplex2(n=n,p=p,r=r,seed=seed) #1, 5, 10, 17
    X <- data$X2
    
    fit1 <- png.ppca_qp(X, nrank=2, eps=1e-12)
    fit2 <- png.gppca_qp(X, nrank=3, maxit=500, gamma=1, kappa=1e-4, eps=1e-12)
    
    fit1$fit.path[[2]]$crit.path %>% plot(type="l")
    fit2$fit.path[[2]]$crit.path %>% plot(type="l")
    fit2$fit.path[[3]]$crit.path %>% plot(type="l")
    
    
    
    # train
    n=150; p=100; r=5;
    seed <- 5
    train <- sim.simplex2(n=n,p=p,r=r,snr=1,sigma=1/p,
                         seed.U=seed, seed.V=seed) #1, 5, 10, 17
    test <- sim.simplex2(n=10*n,p=p,r=r,snr=1,sigma=1/p,
                         seed.U=seed+1, seed.V=seed)
    
    fit <- purrr::map(1:10, ~png.gppca_qp(train$X2, nrank=.x, epsilon=1e-4, kappa=1e-4))
    purrr::map_dbl(fit, ~{norm(.x$X - .x$xhat, "F")})
    purrr::map_dbl(fit, ~{norm(test$X2 - proj(test$X2, .x$vhat), "F")})
    
  }
  
  
  
  n=nrow(X); p=ncol(X); r=nrank
  mu=colMeans(X);
  
  
  if( p < r ) stop("r should be less than or equal to p !")
  
  
  start <- proc.time()
  
  U_total <- V_total <- NULL
  fit.path <- NULL
  for( iter in 1:nrank ){
    if( iter == 1 ){
      fit.rank1 <- png.rank1(X, maxit=maxit, eps=eps, gamma=gamma)
      
      U_total <- cbind(U_total, fit.rank1$uhat)
      V_total <- cbind(V_total, fit.rank1$vhat)
      
      fit.path[[iter]] <- fit.rank1
    } else {
      fit.UV <- UV_update(X, V_total, maxit=maxit, eps=eps, kappa=kappa, gamma=gamma)
      
      U_total <- fit.UV$uhat
      V_total <- fit.UV$vhat
      
      if(!save.est.path){
        fit.UV$est.path <- NULL
      }
      
      fit.path[[iter]] <- fit.UV
    }
  }
  
  end <- proc.time()
  
  colnames(U_total) <- colnames(V_total) <- NULL
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(U_total, V_total)
  xhat <- ifelse(abs(xhat)<1e-10,0,xhat)
  
  
  params=list(nrank=nrank, 
              maxit=maxit, 
              eps=eps, 
              kappa=kappa, 
              gamma=gamma)
  
  return(list(mu=mu, uhat=U_total, vhat=V_total, xhat=xhat, X=X, fit.path=fit.path, fit.rank1=fit.rank1, maxit=maxit, time=end-start, params=params, method="gppca_qp"))
  
}










#' @export png.rank1
png.rank1 <- function(X, maxit=100, eps=1e-4, gamma=0.2){
  
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
      Unew[i]<-onedimconvexprojection(C[i,],X[i,],Vold)
    }
    Unew <- Unew * (1-gamma/it)
    
    Vnew <- V_update(X, Uhat=matrix(0,n,1), Vhat=matrix(0,p,1), Uk=as.matrix(Unew))
    Vnew <- Vnew/sqrt(sum(Vnew^2))
    
    est.path[[it]] <- list(uhat=Unew, vhat=Vnew)
    crit.path[it] <- min( sqrt(mean((Vnew-Vold)^2)), sqrt(mean((Vnew+Vold)^2)) )
    
    if( crit.path[it] < eps ){
      break;
    }
    
  }
  
  # Unew<-rep(0,n)
  # for(i in 1:n){
  #   Unew[i]<-onedimconvexprojection(C[i,],X[i,],Vnew)
  # }
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(Unew,Vnew)
  
  return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, est.path=est.path, crit.path=crit.path) )
  
}




















#' @export UV_update
UV_update <- function(X, Vhat, maxit=100, eps=1e-4, kappa=1e-8, gamma=0.2){
  if(FALSE){
    # UV_update(X, V_total, eps=eps)
    Vhat=V_total;  eps=1e-4;  maxit=100
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
  
  # Vhat + Initial V0
  V0=prcomp( X - tcrossprod(rep(1,n),mu) - X %*% tcrossprod(Vhat) )$rotation[,1]
  
  # U-update
  Unew <- matrix(0,n,r)
  for( i in 1:n ){
    x=as.vector(X[i,]); V=cbind(Vhat,V0)
    Unew[i,] <- solve.QP(Dmat=diag(rep(1,r)), 
                         dvec=t(V) %*% (x-mu), 
                         Amat=t(V), 
                         bvec=-mu )$solution
  }
  
  # V-update
  Vk_new <- V_update(X, Uhat=Unew[,1:(r-1)], Vhat=Vhat, Uk=Unew[,r], kappa=kappa)
  
  est.path <- crit.path <- NULL
  for( it in 1:maxit ){
    Vk_old <- Vk_new
    U_old <- Unew
    
    # U-update
    Unew <- matrix(0,n,r)
    for( i in 1:n ){
      x=as.vector(X[i,]);  V=cbind(Vhat,Vk_new)
      
      Unew[i,] <- solve.QP(Dmat=diag(rep(1,r)), 
                           dvec=t(V) %*% (x-mu), 
                           Amat=t(V), 
                           bvec=-mu )$solution
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
  
  # Unew <- matrix(0,n,r);
  # for( i in 1:n ){
  #   Unew[i,] <- solve.QP(Dmat=diag(rep(1,r)), 
  #                        dvec=t(Vnew) %*% (as.vector(X[i,])-mu), 
  #                        Amat=t(Vnew), 
  #                        bvec=-mu )$solution
  # }
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(Unew,Vnew)
  
  return( list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, est.path=est.path, crit.path=crit.path) )
  
}







#' @export V_update
V_update <- function(X, Uhat, Vhat, Uk, kappa=1e-8){
  if(FALSE){
    Uhat=matrix(0,n,1);  Vhat=matrix(0,p,1);  Uk=as.matrix(Unew)
    Uhat=Unew[,1:(r-1)];  Uk=Unew[,r]
    kappa=1e-10
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
  
  Xhat <- tcrossprod(rep(1,n),mu) + Uhat %*% t(Vhat)
  # png.plot.pca3d(mu,Uhat,Vhat)
  
  Vnew <- rep(0,p)
  Dmat <- diag(rep(1,p))
  dvec <- t(X-Xhat) %*% Uk/sum(Uk^2)
  # png.plot.pca3d(mu,Uhat,dvec)
  # png.plot.pca3d(mu,cbind(Uhat,Uk),cbind(Vhat,dvec))
  
  positiveindex = which(Uk>0); negativeindex = which(Uk<0)
  n.pos <- length(positiveindex)
  n.neg <- length(negativeindex)
  if( n.pos>0 && n.neg>0 ){
    lbmat<- -diag(1/Uk[positiveindex],n.pos,n.pos) %*% Xhat[positiveindex,]
    lb<-apply(lbmat,2,max)
    
    ubmat <- -diag(1/Uk[negativeindex],n.neg,n.neg) %*% Xhat[negativeindex,]
    ub<-apply(ubmat,2,min)
    
    bvec <- c(rep(0,r+1),(lb-kappa),-(ub+kappa))
    
    Amat <- cbind(rep(1,p), Vhat, diag(rep(1,p)), -diag(rep(1,p)))
    
    Vnew <- solve.QP(Dmat, dvec, Amat, bvec, meq=r+1, factorized=FALSE )$solution
    # idx_tmp <- c(1:2,(1:10)[-c(7,10)]+2,(1:10)[-c(7,10)]+12)
    # idx_tmp <- c(1:2,(1:10)[-c(10)]+2,(1:10)[-c(10)]+12)
    # solve.QP(Dmat, dvec, Amat[,idx_tmp], bvec[idx_tmp], meq=r+1, factorized=FALSE )$solution
    # 
    # idx_tmp <- c(1:2,(1:p)[-c(1)]+2,(1:p)[-c(1)]+p+2)
    # solve.QP(Dmat, dvec, Amat[,idx_tmp], bvec[idx_tmp], meq=r+1, factorized=FALSE )$solution
  }
  
  if (n.pos>0 && n.neg==0){
    lbmat<- -diag(1/Uk[positiveindex],n.pos,n.pos) %*% Xhat[positiveindex,] 
    
    lb<-apply(lbmat,2,max)
    
    bvec<-c(rep(0,r+1),(lb-kappa))
    
    Amat<-cbind(rep(1,p), Vhat, diag(rep(1,p)))
    
    Vnew<-solve.QP(Dmat, dvec, Amat, bvec, meq=r+1, factorized=FALSE )$solution
  }
  
  if (n.pos==0 && n.neg>0){
    ubmat<- -diag(1/Uk[negativeindex],n.neg,n.neg) %*% Xhat[negativeindex,] 
    
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
