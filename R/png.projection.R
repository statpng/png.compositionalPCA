#' @export png.projection
png.projection <- function(X, fit, method=c("pca","ppca", "gppca", "lrpca")){
  if(FALSE){
    X=Xtrain; method=fit$method
  }
  
  n=nrow(X); p=ncol(X); r=NCOL(fit$uhat)
  mu=fit$mu
  vhat=fit$vhat
  
  if( method == "pca" ){
    
    xhat2 <- (X-tcrossprod(rep(1,n),fit$mu)) %*% tcrossprod(fit$vhat)
    xhat <- t(apply(xhat2,1,png.proj2simplex))
    
    return(xhat)
    
  } else if( method == "ppca" ){
    
    uhat <- matrix(0,n,r)
    for( k in 1:r ){
      chat <- tcrossprod(rep(1,n), mu) + tcrossprod(uhat[,1:(k-1)], vhat[,1:(k-1)])
      for( i in 1:n ){
        uhat[i,k] <- onedimconvexprojection(chat[i,], as.vector(X[i,]), vhat[,k])
      }
      # it=fit$fit.path[[k]]$it;  gamma=fit$params$gamma
      # uhat[,k] <- uhat[,k] * (1-gamma/it)
    }
    
  } else if( method == "gppca" ){
    
    uhat <- matrix(0,n,r)
    for( i in 1:n ){
      uhat[i,] <- multidimconvexprojection(mu, as.vector(X[i,]), vhat)
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
