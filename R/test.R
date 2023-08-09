if(FALSE){
  rm(list=ls())
  devtools::load_all("./png.compositionalPCA")
}


#' @export png.test
png.test <- function(){
  
  n <- 100; p <- 50; r <- 5
  data <- sim.simplex(n=n, p=p, r=r, snr=0.1, d=100, eta=0.2, seed=1)
  gamma.seq <- c(10^(-seq(1,6,1)),0)
  
  fit1 <- purrr::map(gamma.seq, ~{
    print(.x);  png.ppca_qp(data$X2, nrank=r, gamma=.x, V.init="PC")
  })
  fit2 <- purrr::map(gamma.seq, ~{
    print(.x);  png.gppca_qp(data$X2, nrank=r, gamma=.x, V.init="PC")
  })
  
  fit1[[4]] %>% png.crit.path()
  fit2[[4]] %>% png.crit.path()
  fit2[[5]] %>% png.crit.path()
  fit2[[6]] %>% png.crit.path()
  fit2[[7]] %>% png.crit.path()
  
  
  fit1 %>% sapply(function(fit){
    sqrt(mean((data$X0 - fit$xhat)^2))
  })
  fit2 %>% sapply(function(fit){
    sqrt(mean((data$X0 - fit$xhat)^2))
  })
  
  fit0 <- png.gppca_qp(data$X2, nrank=r, kappa=0, maxit=1000, gamma=0, V.init="PC") 
  fit <- png.gppca_qp(data$X2, nrank=r, kappa=0, maxit=1000, gamma=0.5, V.init="random", phi=0.2)
  
  fit <- png.gppca_qp(data$X2, nrank=r, kappa=0, gamma=0.01, V.init="PC") 
  fit %>% png.crit.path()
  fit %>% {sqrt(mean((data$X0 - .$xhat)^2))}
  
  
  #
  #
  #
  
  png.angle(fit0$vhat[,1:2], fit$vhat[,1:2])
  
  
  fit2[[1]] %>% png.crit.path()
  #
  
  
  
  GRID <- expand.grid(gamma.seq=gamma.seq, kappa.seq=kappa.seq)
  
  fit1 <- lapply(1:nrow(GRID), function(i){
    print(i)
    kappa=GRID$kappa.seq[i];  gamma=GRID$gamma.seq[i[]]
    png.ppca_qp(data$X2, nrank=r, kappa=kappa, gamma=gamma, V.init="PC", verbose=FALSE)
  })
  
  fit_grid2 <- lapply(1:nrow(GRID), function(i){
    print(i)
    kappa=GRID$kappa.seq[i];  gamma=GRID$gamma.seq[i[]]
    png.gppca_qp(data$X2, nrank=r, kappa=kappa, gamma=gamma, V.init="PC", verbose=FALSE)
  })
  
  fit1 %>% sapply(function(fit){
    sqrt(mean((data$X0 - fit$xhat)^2))
  })
  fit_grid2 %>% sapply(function(fit){
    sqrt(mean((data$X0 - fit$xhat)^2))
  })
  
  fit1[[2]] %>% png.crit.path()
  fit2[[30]] %>% png.crit.path()
  
  fit <- png.gppca_qp(data$X2, nrank=r, kappa=0, gamma=0.8, V.init="PC", verbose=FALSE)
  list(fit) %>% sapply(function(fit){
    sqrt(mean((data$X0 - fit$xhat)^2))
  })
  fit %>% png.crit.path(remove = 1:10)
  
  fit1 %>% png.crit.path()
  fit2[[1]] %>% png.crit.path()
  
  fit2[[1]] %>% png.crit.path()
  
  
  
}
