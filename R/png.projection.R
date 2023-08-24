#' @export png.projection
png.projection <- function(X, fit, nrank=NULL, method=c("pca","ppca", "gppca", "lrpca")){
  if(FALSE){
    X=Xtrain; method=fit$method
  }
  
  n=nrow(X); p=ncol(X); 
  
  if(is.null(nrank)){
    r=NCOL(fit$uhat)
  } else {
    r=nrank
    if(method == "lrpca"){
      fit$logvhat <- fit$logvhat[,1:r,drop=F]
    } else {
      fit$vhat <- fit$vhat[,1:r,drop=F]
    }
  }
  
  
  mu=fit$mu
  vhat=fit$vhat
  
  if( method == "pca" ){
    
    xhat2 <- (X-tcrossprod(rep(1,n),fit$mu)) %*% tcrossprod(fit$vhat)
    xhat <- t(apply(xhat2,1,png.proj2simplex))
    
    return(xhat)
    
  } else if( method == "ppca" ){
    
    uhat <- matrix(0,n,r)
    for( k in 1:r ){
      if(k-1 == 0){
        chat <- tcrossprod(rep(1,n), mu)
      } else {
        chat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat[,1:(k-1)], vhat[,1:(k-1)])
      }
      
      for( i in 1:n ){
        # uhat[i,k] <- onedimconvexprojection(chat[i,], as.vector(X[i,]), vhat[,k])
        uhat[i,k] <- Solve_U_SP(as.vector(X[i,]), chat[i,], vhat[,k], gamma=0)
      }
      # it=fit$fit.path[[k]]$it;  gamma=fit$params$gamma
      # uhat[,k] <- uhat[,k] * (1-gamma/it)
    }
    
  } else if( method == "gppca" ){
    
    uhat <- matrix(0,n,r)
    for( i in 1:n ){
      # uhat[i,] <- multidimconvexprojection(mu, as.vector(X[i,]), vhat)
      uhat[i,] <- Solve_U_GP(as.vector(X[i,]), mu, vhat, gamma=0)
    }
    # it=fit$fit.path[[r]]$it;  gamma=fit$params$gamma
    # uhat <- uhat * (1-gamma/it)
    
  } else if( method == "lrpca" ){
    
    f <- switch(fit$zero.replace, 
                "simple"=png.ZeroReplace.simple,
                "additive"=png.ZeroReplace.additive,
                "multiplicative"=png.ZeroReplace.multiplicative)
    
    Xnew <- t(apply(X, 1, f, delta=fit$delta))
    Xclr <- t(apply(Xnew, 1, png.clr))
    
    mu <- t(apply(fit$Xnew, 1, png.clr)) %>% colMeans()
    vhat <- fit$logvhat
    uhat <- (Xclr - tcrossprod(rep(1,n),mu)) %*% vhat
    
  }
  
  xhat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat, vhat)
  
  if( method == "lrpca" ){
    xhat <- xhat %>% {t(apply(.,1,png.iclr))}
  }
  
  return(xhat)
}













#' @export png.projection2
png.projection2 <- function(X, fit, nrank=NULL, method=c("pca_proj","ppca", "gppca", "lrpca")){
  if(FALSE){
    X=Xtrain; method=fit$method
  }
  
  n=nrow(X); p=ncol(X); 
  
  if(is.null(nrank)){
    r=NCOL(fit$uhat)
  } else {
    r=nrank
    if(method == "lrpca"){
      fit$logvhat <- fit$logvhat[,1:r,drop=F]
    } else {
      fit$vhat <- fit$vhat[,1:r,drop=F]
    }
  }
  
  
  mu=fit$mu
  vhat=fit$vhat
  
  if( method == "pca_proj" ){
    
    xhat2 <- (X-tcrossprod(rep(1,n),fit$mu)) %*% tcrossprod(fit$vhat)
    xhat <- t(apply(xhat2,1,png.proj2simplex))
    
    return(list(mu=mu, uhat=uhat, vhat=vhat, xhat=xhat))
    
  } else if( method == "ppca" ){
    
    uhat <- matrix(0,n,r)
    for( k in 1:r ){
      if(k-1 == 0){
        chat <- tcrossprod(rep(1,n), mu)
      } else {
        chat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat[,1:(k-1)], vhat[,1:(k-1)])
      }
      
      for( i in 1:n ){
        # uhat[i,k] <- onedimconvexprojection(chat[i,], as.vector(X[i,]), vhat[,k])
        uhat[i,k] <- Solve_U_SP(as.vector(X[i,]), chat[i,], vhat[,k], gamma=0)
      }
      # it=fit$fit.path[[k]]$it;  gamma=fit$params$gamma
      # uhat[,k] <- uhat[,k] * (1-gamma/it)
    }
    
  } else if( method == "gppca" ){
    
    uhat <- matrix(0,n,r)
    for( i in 1:n ){
      # uhat[i,] <- multidimconvexprojection(mu, as.vector(X[i,]), vhat)
      uhat[i,] <- Solve_U_GP(as.vector(X[i,]), mu, vhat, gamma=0)
    }
    # it=fit$fit.path[[r]]$it;  gamma=fit$params$gamma
    # uhat <- uhat * (1-gamma/it)
    
  } else if( method == "lrpca" ){
    
    f <- switch(fit$zero.replace, 
                "simple"=png.ZeroReplace.simple,
                "additive"=png.ZeroReplace.additive,
                "multiplicative"=png.ZeroReplace.multiplicative)
    
    Xnew <- t(apply(X, 1, f, delta=fit$delta))
    Xclr <- t(apply(Xnew, 1, png.clr))
    
    mu <- t(apply(fit$Xnew, 1, png.clr)) %>% colMeans()
    vhat <- fit$logvhat
    uhat <- (Xclr - tcrossprod(rep(1,n),mu)) %*% vhat
    
  }
  
  xhat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat, vhat)
  
  if( method == "lrpca" ){
    xhat <- xhat %>% {t(apply(.,1,png.iclr))}
  }
  
  return(list(mu=mu, uhat=uhat, vhat=vhat, xhat=xhat))
}
