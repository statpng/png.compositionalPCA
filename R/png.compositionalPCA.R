#' @export onedimconvexprojection
onedimconvexprojection<-function(c,x,v){
  # one-dimensional convex projection function
  
  # c,x is an element of S^{p} and v is a unit vector perpendicular with 1=(1,1,...,1)\in \mathbb{R}^{p}
  
  # direct projection score
  t=sum((x-c)*v)
  
  # define index set
  indexminus<-which(v<0)
  indexplus<-which(v>0)
  
  lminus<-length(indexminus)
  lplus<-length(indexplus)
  
  # calculate convec projection score case by case
  if (lminus>0 & lplus>0){
    m=max( -(c/v)[indexplus])
    M=min( -(c/v)[indexminus])
    return(min(max(m,t),M))
  }
  
  if (lminus>0 & lplus==0){
    m=max( -(c/v)[indexplus])
    return(max(m,t))
  }
  
  if (lminus==0 & lplus>0){
    M=min( -(c/v)[indexminus])
    return(min(t,M))
  }
  
  if (lminus==0 & lplus==0){
    return(t)
  }
  
}

# # example of onedim convex projection
# {
#   p<-100
#   c<-rep(1,p)/p
#   x<-runif(p,0,1)
#   x<-x/sum(x)
#   v<-rnorm(p,0,1)
#   v<-v-sum(v*1)/p
#   v<-v/sqrt(sum(v^2))
#   
#   t<-onedimconvexprojection(c,x,v)
# }




#' @export multidimconvexprojection
multidimconvexprojection<-function(c,x,V){
  # c,x is an element of S^{p} and V=(V1, V2, ..., Vr) is r number of orthonormal vectors perpendicular with 1=(1,1,...,1)\in \mathbb{R}^{p}, given by V\in \mathbb{R}^{p \times r}
  
  if(FALSE){
    set.seed(1)
    n=1000; p=3; r=2
    X <- sim.simplex2(n=n, p=p, r=r)$X2
    
    
    mu <- colMeans(X)
    vhat <- prcomp(X)$rotation[,1:(ncol(X)-1)]
    x <- X[1,]
    
    mu + tcrossprod( multidimconvexprojection(mu,x,vhat[,1:2]), vhat[,1:2] )
    
    c=C[i,]; x=X[i,]; V=Vold
  }
  
  
  require(quadprog)
  
  r<-ncol(V)
  U<-solve.QP(Dmat=diag(rep(1,r)), dvec=t(V) %*% (x-c), Amat=t(V), bvec=-c )$solution
  return(U)
  
}











#' @export png.pca
png.pca <- function(X, nrank=2){
  d0 <- sum(svd(X)$d>1e-10)
  # if( nrank > (d0-1) ) stop("nrank should be smaller than or equal to rank(X)-1")
  
  n=nrow(X); p=ncol(X); mu=colMeans(X)
  uhat <- prcomp(X)$x[,1:nrank,drop=F]
  vhat <- prcomp(X)$rotation[,1:nrank,drop=F]
  xhat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat, vhat)
  
  return( list(mu=mu, uhat=uhat, vhat=vhat, xhat=xhat, X=X, method="pca") )
}


#' @export png.lrpca
png.lrpca <- function(X, nrank=2, zero.replace=NULL, delta=1e-6){
  
  n=nrow(X); p=ncol(X); 
  
  if(!is.null(zero.replace)){
    f <- switch(zero.replace, 
                "simple"=png.ZeroReplace.simple,
                "additive"=png.ZeroReplace.additive,
                "multiplicative"=png.ZeroReplace.multiplicative)
    Xnew <- t(apply(X, 1, f, delta=delta))
  } else {
    Xnew <- X
  }
  
  mu = png.iclr( colMeans( log(Xnew) ) )
  
  Xclr <- t(apply(Xnew,1,png.clr))
  fit <- png.pca(Xclr, nrank=nrank)
  uhat <- fit$uhat
  vhat <- lapply( seq(-20,20,0.5), function(z){
    z*fit$vhat
  })
  xhat <- fit$xhat %>% {t(apply(.,1,png.iclr))}
  
  return( list(mu=mu, uhat=uhat, vhat=vhat, logvhat=fit$vhat, xhat=xhat, X=X, Xnew=Xnew, method="lrpca", zero.replace=zero.replace, delta=delta) )
}


#' @export png.ppca
png.ppca <- function(X, nrank=2){
  
  if(FALSE){
    set.seed(2)
    n=500; p=4; r=2
    X <- sim.simplex2(n=n, p=p, r=r, snr=1, d=1, d0=10)$X2
    png.quaternary3d(X,vhat=png.ppca(X,2)$vhat, xhat=png.ppca(X,2)$xhat)
    
    png.ppca(X,2)$xhat %>% {sum(.<0)}
    
    png.ppca(X,2)$xhat[png.ppca(X,2)$xhat<0]
  }
  
  
  d0 <- sum(svd(X)$d>1e-10)
  
  if( nrank > (d0-1) ) stop("nrank should be smaller than rank(X)-1")
  
  
  n=nrow(X); p=ncol(X)
  mu <- colMeans(X)
  vhat <- prcomp(X)$rotation[,1:nrank,drop=F]
  
  
  uhat <- NULL
  for( j in 1:nrank ){
    if(j==1) Xnew=X
    v <- vhat[,j,drop=F]
    uhat_j <- apply(Xnew,1,function(x) onedimconvexprojection(mu, x, v))
    uhat <- cbind(uhat, uhat_j)
    Xnew <- Xnew #- tcrossprod(uhat,vhat[,1:j])
  }
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(uhat,vhat)
  
  return( list(mu=mu, uhat=uhat, vhat=vhat, xhat=xhat, X=X, method="ppca") )
}




#' @export png.gppca
png.gppca <- function(X, nrank=2, V=prcomp(X)$rotation[,1:nrank,drop=F]){
  
  if(FALSE){
    set.seed(3)
    n=500; p=4; r=2
    X <- sim.simplex2(n=n, p=p, r=r, delta=runif(1,0.8,1.2))$X2
    png.ppca(X,2)$xhat %>% {sum(.<(-1e-10))}
    png.gppca(X,2)$xhat %>% {sum(.<(-1e-10))}
    nrank=1
  }
  
  
  d0 <- sum(svd(X)$d>1e-10)
  
  if( nrank > (d0-1) ) stop("nrank should be smaller than rank(X)-1")
  
  
  mu <- colMeans(X)
  n <- nrow(X)
  
  uhat <- apply(X,1,function(x){
    multidimconvexprojection(mu,x,V)
  }) %>% {matrix(., n, nrank, byrow=T)}
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(uhat, V)
  
  return( list(mu=mu, uhat=uhat, vhat=V, xhat=xhat, X=X, method="gppca") )
}








# png.SQP <- function(X,V,U,V0,eta=10^-4){
#   
#   # X is n\times p data matrix and u\in S^{p}, V=(\hat{V}_{1}, \cdots, \hat{V}_{k-1}) are given direction vectors, U=(\hat{U}_{1}, \cdots, \hat{U}_{k-1}) are given projection scores, V0 is an initial value of V_{k}
#   
#   require(quadprog)
#   
#   count<-0
#   n<-nrow(X)
#   p<-ncol(X)
#   d<-ncol(U)
#   mu <- colMeans(X)
#   
#   
#   C<-tcrossprod(rep(1,n),mu)+U%*%t(V)
#   
#   Vold<-V0
#   
#   while(TRUE){
#     count<-count+1
#     Unew<-rep(0,n)
#     for(i in 1:n){
#       Unew[i]<-onedimconvexprojection(C[i,],X[i,],Vold)
#     }
#     
#     V_update(X, Uhat=U, Vhat=as.matrix(V), Uk=as.matrix(Unew))
#     
#     Vnew<-rep(0,p)
#     
#     Dmat<- diag(rep(1,p))
#     dvec<-t(X-C)%*%Unew/sum(Unew^2)
#     
#     positiveindex<-which(Unew>0)
#     negativeindex<-which(Unew<0)
#     
#     if (length(positiveindex)>0 && length(negativeindex)>0){
#       lbmat<- -diag(1/Unew[positiveindex]) %*% C[positiveindex,] 
#       
#       lb<-apply(lbmat,2,max)
#       
#       ubmat<- -diag(1/Unew[negativeindex]) %*% C[negativeindex,] 
#       
#       ub<-apply(ubmat,2,min)
#       
#       
#       bvec<-c(rep(0,d+1),lb,-ub)
#       
#       Amat<-cbind(one, V, diag(rep(1,p)) , -diag(rep(1,p)))
#       
#       Vnew<-solve.QP(Dmat, dvec, Amat, bvec, meq=d+1,factorized=FALSE )$solution
#     }
#     
#     if (length(positiveindex)>0 && length(negativeindex)==0){
#       lbmat<- -diag(1/Unew[positiveindex]) %*% C[positiveindex,] 
#       
#       lb<-apply(lbmat,2,max)
#       
#       bvec<-c(rep(0,d+1),lb)
#       
#       Amat<-cbind(one, V, diag(rep(1,p)))
#       
#       Vnew<-solve.QP(Dmat, dvec, Amat, bvec, meq=d+1,factorized=FALSE )$solution
#     }
#     
#     if (length(positiveindex)==0 && length(negativeindex)>0){
#       ubmat<- -diag(1/Unew[negativeindex]) %*% C[negativeindex,] 
#       
#       ub<-apply(ubmat,2,min)
#       
#       bvec<-c(rep(0,d+1),-ub)
#       
#       Amat<-cbind(one, V, -diag(rep(1,p)))
#       
#       Vnew<-solve.QP(Dmat, dvec, Amat, bvec, meq=d+1,factorized=FALSE )$solution
#     }
#     
#     Vnew<-Vnew/sqrt(sum(Vnew^2))
#     
#     if (min( sum((Vnew-Vold)^2) ,sum((Vnew+Vold)^2) )<eta){
#       break
#     }
#     Vold<-Vnew
#     
#     
#   }
#   
#   Vold
#   Unew
#   count
#   result<-list()
#   result[[1]]<-Unew
#   result[[2]]<-Vold
#   result[[3]]<-count
#   names(result)<-c("U","V","count")
#   return(result)
#   
# }







#' @export png.ppca_qp
png.ppca_qp <- function(X, nrank=2, maxit=500, eps=1e-6, kappa=1e-4, gamma=0, phi=0.01, V.init=c("PC","random"), verbose=TRUE){
  if(FALSE){
    set.seed(12)
    nrank=5
    X <- sim.simplex(500,50,r=nrank,eta=0.2)$X2
    fit1 <- png.ppca_qp(X, nrank=nrank, V.init="PC")
    fit2 <- png.gppca_qp(X, nrank=nrank, gamma=0.5, V.init="PC")
    
    fit1 %>% png.crit.path()
    fit2 %>% png.crit.path()
    
    fit1$fit.path[[3]]$xhat[1:5,1:5]
    
    fit1$fit.path[[3]]$xhat[1:5,1:5]
    fit2$fit.path[[3]]$xhat[1:5,1:5]
    
    png.quarternary3d(X)
    png.plot.pca3d(fit$mu, fit$uhat, fit$vhat)
    fit$vhat
    
    nrank=r; kappa=1e-6; maxit=500; eps=1e-6; gamma=0.5; save.est.path = TRUE
  }
  if(FALSE){
    nrank=5; maxit=500; eps=1e-6; kappa=1e-4; gamma=0.5; V.init=c("PC","random")[1]
  }
  
  library(quadprog)
  
  n=nrow(X); p=ncol(X); r=nrank
  mu=colMeans(X);
  
  
  Uhat <- Vhat <- NULL
  fit.path <- NULL
  for( iter in 1:nrank ){
    if( iter == 1 ){
      fit.rank1 <- png.rank12(X, maxit=maxit, eps=eps, kappa=kappa, gamma=gamma, phi=phi, V.init=V.init, verbose=verbose)
      
      Uhat <- cbind(Uhat, fit.rank1$uhat)
      Vhat <- cbind(Vhat, fit.rank1$vhat)
      
      fit.path[[iter]] <- fit.rank1
    } else {
      
      fit.UkVk <- update_UkVk(X=X, Uhat=Uhat, Vhat=Vhat, maxit=maxit, eps=eps, kappa=kappa, gamma=gamma, phi=phi, V.init=V.init, verbose=verbose)
      
      Uhat <- cbind(Uhat, fit.UkVk$uhat)
      Vhat <- cbind(Vhat, fit.UkVk$vhat)
      
      # if(!save.est.path){
      #   fit.UkVk$est.path <- NULL
      # }
      
      fit.path[[iter]] <- fit.UkVk
    }
  }
  
  xhat <- tcrossprod(rep(1,n),mu) + tcrossprod(Uhat, Vhat)
  
  params=list(nrank=nrank, 
              maxit=maxit, 
              eps=eps, 
              kappa=kappa, 
              gamma=gamma)
  
  return(list(mu=mu, uhat=Uhat, vhat=Vhat, xhat=xhat, X=X, fit.path=fit.path, maxit=maxit, method="ppca", params=params))
  
}













#' @export update_UkVk
update_UkVk <- function(X, Uhat, Vhat, maxit=500, eps=1e-6, kappa=1e-4, gamma=0, phi=0.01, V.init=c("PC","random"), verbose=TRUE){
  
  # X is n\times p data matrix and u\in S^{p}, V=(\hat{V}_{1}, \cdots, \hat{V}_{k-1}) are given direction vectors, U=(\hat{U}_{1}, \cdots, \hat{U}_{k-1}) are given projection scores, V0 is an initial value of V_{k}
  
  if(FALSE){
    X=X
    V=V_total
    U=U_total
    V0=V0
    eta=eta
  }
  
  require(quadprog)
  
  n=nrow(X); p=ncol(X); d=NCOL(Uhat)
  mu=colMeans(X)
  
  C <- tcrossprod(rep(1,n), mu) + tcrossprod(Uhat, Vhat)
  # C <- tcrossprod(rep(1,n), mu) + X %*% tcrossprod(Vhat)
  
  # Initial V0
  V0.PC <- prcomp(X-(tcrossprod(rep(1,n), mu) + X %*% tcrossprod(Vhat)))$rotation[,1]
  if(V.init == "PC"){
    V0 <- V0.PC
  } else {
    V0 <- (V0.PC + qr.Q(qr(cbind(1,rnorm(ncol(X)))))[,2]*phi) %>% {./norm(.,"2")}
  }
  
  png.angle(V0.PC, V0)[[1]] %>% print
  
  
  # est.path <- 
  crit.path <- NULL
  for( it in 1:maxit ){
    if(verbose) print(paste0("it: ", it))
    
    if( it == 1 ){
      Vold <- V0
    } else {
      Vold <- Vnew
    }
    
    Unew <- sapply(1:n, function(i) Solve_U_SP(x=X[i,], mu=C[i,], v=Vold, gamma=gamma/(it^(1/2))))
    Vnew <- V_update2(X, Uhat=Uhat, Vhat=as.matrix(Vhat), Uk=as.matrix(Unew), kappa=kappa)
    if( Vnew[1] == "No solution" ){
      Vnew <- Vold
      NoSolution <- TRUE
    } else {
      NoSolution <- FALSE
    }
    
    # est.path[[it]] <- list(uhat=cbind(Uhat, Unew), vhat=cbind(Vhat, Vnew))
    crit.path[it] <- min( sqrt(mean((Vnew-Vold)^2)), sqrt(mean((Vnew+Vold)^2)) )
    
    if( crit.path[it] < eps ) break
    
  }
  
  if(NoSolution) crit.path[it] <- 0.1
  
  
  Unew <- sapply(1:n, function(i) Solve_U_SP(x=X[i,], mu=C[i,], v=Vnew, gamma=0))
  
  xhat <- C + tcrossprod(Unew,Vnew)
  
  loss <- sqrt(mean((X - xhat)^2))
  
  # result <- list(uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, est.path=est.path)
  result <- list(xhat=xhat, uhat=Unew, vhat=Vnew, it=it, crit.path=crit.path, loss=loss)
  
  return(result)
  
}
