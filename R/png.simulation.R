#' @export EstimationPerformance
EstimationPerformance <- function(){
  
  lapply(1:3, function(nrank){
    X = cpca.sim1(n=500,p=4,r=2)$X
    muhat = rep(1,nrow(X)) %*% t(colMeans(X))
    Xhat = with( prcomp(X), muhat+tcrossprod(x[,1:nrank], rotation[,1:nrank]) )
    prcomp(X)
    norm(X-Xhat,"F")
  })
  
  
  lapply(1:3, function(nrank){
    data = cpca.sim2(n=500,p=4,r=3,rho=0.5)
    X = data$X
    muhat = rep(1,nrow(X)) %*% t(colMeans(X))
    fit.pca = prcomp(X)
    Xhat = with( fit.pca, muhat+tcrossprod(x[,1:nrank], rotation[,1:nrank]) )
    norm(X-Xhat,"F")
    
    # png.utils::png.angle(fit.pca$rotation, data$V)
    # png.utils::png.CompareSubspace(fit.pca$rotation, data$V)
    # fit.pca$rotation
    # data$V
  })
  
}


