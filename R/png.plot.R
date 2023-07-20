#' @export png.pca.convergence
png.pca.convergence <- function(fit){
  if(FALSE){
    fit <- fit5_1
  }
  n=nrow(fit$X); p=ncol(fit$X); r=NCOL(fit$vhat)
  mu=fit$mu
  maxit=fit$maxit
  
  crit.array <- array(NA, 
                      dim=c(maxit-1, r, 3), 
                      dimnames=list(paste0("iter=",1:(maxit-1)), 
                                    paste0("rank=",1:r), 
                                    c("uhat", "vhat", "xhat") ))
  for( k in 1:r ){
    est.path <- fit$fit.path[[k]]$est.path
    
    if( length(est.path) == 1 ){ 
      crit.array[1,k,1] <- 1e-16
      crit.array[1,k,2] <- 1e-16
      crit.array[1,k,3] <- 1e-16
      
      next;
    }
    
    for( kk in 1:(length(est.path)-1) ){
      # kk
      uhat_old <- as.matrix(est.path[[kk]]$uhat)
      vhat_old <- as.matrix(est.path[[kk]]$vhat)
      xhat_old <- tcrossprod(rep(1,n),mu)+tcrossprod(uhat_old,vhat_old)
      
      # kk+1
      uhat_new <- as.matrix(est.path[[kk+1]]$uhat)
      vhat_new <- as.matrix(est.path[[kk+1]]$vhat)
      xhat_new <- tcrossprod(rep(1,n),mu)+tcrossprod(uhat_new,vhat_new)
      
      crit.array[kk,k,1] <- sqrt(sum((uhat_new-uhat_old)^2)/(n*r))
      crit.array[kk,k,2] <- sqrt(sum((vhat_new-vhat_old)^2)/(n*r))
      crit.array[kk,k,3] <- sqrt(sum((xhat_new-xhat_old)^2)/(n*r))
    }
  }
  
  out <- crit.array[,,] %>% plyr::adply(1:3)
  colnames(out) <- c("iter", "rank", "est", "value")
  
  out2 <- subset(out, !is.na(value))
  
  print( out2 %>% group_by(est, rank) %>% slice_tail(n=1) )
  
  return( out2 )
}


#' @export png.pca.plot_convergence
png.pca.plot_convergence <- function(fit, maxit=100){
  if(is.data.frame(fit)){
    df <- fit
  } else {
    df <- png.pca.convergence(fit)
  }
  
  df[,1] <- as.numeric(df[,1])
  
  df %>%
    # filter(est=="xhat") %>% 
    filter(iter>1, iter<=maxit) %>% 
    ggplot() +
    geom_line(aes(iter, value, color=rank)) + 
    ggsci::scale_color_lancet() +
    # scale_y_sqrt() +
    # scale_y_continuous(trans='log2') +
    xlab("Iteration") + ylab("Changes") +
    facet_grid(est~., scales="free") +
    theme_bw()
}










png.sim.simplex.test <- function(params){
  TestData <- params %>% 
    { sim.simplex(n=.["n"], 
                  p=.["p"], 
                  r=.["r"], 
                  snr=.["snr"], 
                  d=.["d"],
                  d0=.["d0"],
                  seed.U=.["seed.U"]*2,
                  seed.V=.["seed.V"],
                  eta=.["eta"]) }
  
  TestData
}



#' @export png.angle
png.angle <- function(true, est){
  # true: n x p; est: n x p
  
  # The largest principal angle
  qt <- qr.Q(qr(true))
  qe <- qr.Q(qr(est))
  fit.svd <- svd( crossprod(qe, qt) )
  theta <- acos(fit.svd$d |> round(12))
  
  # theta[1] * 180 / pi (in degree)
  list( max = theta[1] * 180 / pi, Grassmannian = norm( theta, "2" ) * 180 / pi )
}




#' @export png.pca.criteria
png.pca.criteria <- function(fit, data, n.test){
  if(FALSE){
    n=500; p=50; r=5; snr=2; eta=0.1/log(p); seed=1
    data <- sim.simplex(n=n,p=p,r=r,snr=snr,d=10,d0=0.01,seed=seed,eta=eta)
    
    fit <- png.ppca_qp(data$X2, nrank=r, kappa=1e-5, maxit=500, eps=1e-6, gamma=0.5, save.est.path = TRUE)
    
    png.pca.criteria(fit, data, n.test=n*10)
    
  }
  

  data$params["n"] <- n.test
  true <- data %>% { list(Xtrain=.$X2,
                          Xtest=png.sim.simplex.test(.$params)$X2,
                          V=.$V) }
  
  Xtrain=true$Xtrain
  Xtest=true$Xtest
  Vtrue=true$V
  
  vhat <- fit$vhat
  xhat <- fit$xhat
  xhat_train <- png.projection(Xtrain, fit, method=fit$method)
  xhat_test <- png.projection(Xtest, fit, method=fit$method)
  
  xhat[1:5,1:5]
  xhat_train[1:5,1:5]
  
  rmse.Xtrain <- sqrt(mean((Xtrain-xhat)^2))
  rmse.Xtest <- sqrt(mean((Xtest-xhat)^2))
  Pangle.V <- png.angle(Vtrue, vhat)$max
  Gangle.V <- png.angle(Vtrue, vhat)$Grassmannian
  OutOfSimplex <- mean(apply(xhat,1,function(x) any(x < -1e-10)))
  # Out-of-simplex Sample Percentage
  
  c(rmse.Xtrain=rmse.Xtrain,
    rmse.Xtest=rmse.Xtest,
    Pangle.V=Pangle.V,
    Gangle.V=Gangle.V,
    OutOfSimplex=OutOfSimplex)
  
}
