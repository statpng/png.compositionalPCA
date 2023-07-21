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
















#' @export sim01_convergence
sim01_convergence <- function(n, p, r, snr, eta, kappa, gamma, seed){
  if(FALSE){
    n=100; p=4; r=2; snr=2; eta=0.1/log(p); seed=1
    n=100; p=50; r=5; snr=5; eta=0.0; seed=1
  }
  
  data <- sim.simplex(n=n,p=p,r=r,snr=snr,d=10,d0=0.01,seed=seed,eta=eta/log(p))
  X <- data$X2
  
  fit4 <- png.ppca_qp(X, nrank=r, kappa=kappa, maxit=500, eps=1e-6, gamma=gamma, save.est.path=TRUE)
  fit5 <- png.gppca_qp(X, nrank=r, kappa=kappa, maxit=500, eps=1e-6, gamma=gamma, save.est.path=TRUE)
  
  # png.pca.plot_convergence(fit4)
  # png.pca.plot_convergence(fit5)
  
  return( list(fit4, fit5) %>% map(~png.pca.convergence(.x)) )
  
}





#' @export sim02_kappa_gamma
sim02_kappa_gamma <- function(n, p, r, snr, eta, kappa=1e-2, gamma=0.5, seed){
  if(FALSE){
    idx=1; r=5; snr=2
    n=grid$n.seq[idx]
    p=grid$p.seq[idx]
    eta=grid$eta.seq[idx]
    kappa=grid$kappa.seq[idx]
    gamma=grid$gamma.seq[idx]
    seed=grid$seed.seq[idx]
  }
  data <- sim.simplex(n=n,p=p,r=r,snr=snr,d=10,d0=0.01,seed=seed,eta=eta/log(p))
  X <- data$X2
  
  # kappa: 1e-0, 1e-2, 1e-4, 1e-6, 1e-8
  # gamma: 0, 0.1, 0.2, 0.3, 0.4, 0.5
  fit_ppca <- try(png.ppca_qp(X, nrank=r, kappa=kappa, maxit=500, eps=1e-6, gamma=gamma))
  fit_gppca <- try(png.gppca_qp(X, nrank=r, kappa=kappa, maxit=500, eps=1e-6, gamma=gamma))
  res1 <- try(png.pca.criteria(fit_ppca, data, n.test=n*10))
  res2 <- try(png.pca.criteria(fit_gppca, data, n.test=n*10))
  
  list(ppca=res1, gppca=res2)
}



#' @export png.list.replace
png.list.replace <- function(params, List){
  # list_to_be_replaced
  if(FALSE){
    List=list(kappa.seq=1e-2, gamma.seq=0.5)
  }
  
  for( i in 1:length(List) ){
    params[ names(List)[i] ] <- List[[i]]
  }
  params
}




#' @export run.sim
run.sim <- function(f.sim, params){
  library(parallel)
  library(purrr)
  
  if(FALSE){
    {
      n.seq <- c(50,100,500)
      p.seq <- c(100,200)
      eta.seq <- c(0, 0.1)
      seed.seq <- c(1:50)
      kappa.seq <- 10^(-seq(0,8,2))
      gamma.seq <- seq(0, 0.5, 0.1)
      
      params <- list(n.seq=n.seq,
                         p.seq=p.seq,
                         eta.seq=eta.seq,
                         kappa.seq=kappa.seq,
                         gamma.seq=gamma.seq,
                         seed.seq=seed.seq)
    }
    
    
    
    
    run.sim(sim01_convergence, 
            params=params %>% 
              png.list.replace(list(n.seq=50, 
                                    kappa.seq=1e-2, 
                                    gamma.seq=0.5)))
    
    run.sim(sim02_kappa_gamma, 
            params=params)
    
    
  }  
  
  title <- deparse(substitute(f.sim))
  
  
  grid <- expand.grid(params)
  
  
  result <- mcmapply(function(n, p, r, snr, eta, kappa, gamma, seed){
    print(seed)
    f.sim(n=n, p=p, r=r, snr=snr, eta=eta/log(p), kappa=kappa, gamma=gamma, seed=seed)
    },
    n=grid$n.seq,
    p=grid$p.seq,
    eta=grid$eta.seq,
    kappa=grid$kappa.seq,
    gamma=grid$gamma.seq,
    seed=grid$seed.seq,
    r=5, snr=2,
    mc.cores = parallel::detectCores()-1)
  
  save(result, file=paste0("./",title,".RData"))
}




