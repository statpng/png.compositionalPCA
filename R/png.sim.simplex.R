# if(FALSE){
#   
#   n=1000; p=4; r=2; snr=100; seed=1
#   X <- sim.simplex(n,p,r,snr=snr,d=10,d0=0.1,seed=1028,eta=0.1)$X2
#   mean(X<1e-10)
#   
#   fit <- png.ppca(X,2)
#   png.quaternary3d(X, fit, size=1)
#   #
#   #
#   
#   sim.simplex(n,p,r,snr=snr,d=10,d0=0.01,seed=rpois(1,1000))$X2 %>% 
#     {png.quaternary3d(., png.ppca(.,1))}
#   
#   {
#     sim.simplex2(n,p,r,snr=snr,d=10,d0=0.01,seed=rpois(1,1000))$X2 %>% 
#       {png.quaternary3d(., png.ppca(.,1))}
#   }
#   
#   {
#     lim <- c(-0.5,0.5)/10
#     sim.simplex(n,p=100,r,snr=snr,d=10,d0=0.1,seed=rpois(1,1000))$X2 %>% 
#       {prcomp(.)$x[,1:2]} %>% plot(xlim=lim,ylim=lim)
#   }
#   
#   {
#     lim <- c(-0.5,0.5)/3
#     sim.simplex2(n,p=200,r=3,snr=snr,d=10,d0=0.01,seed=rpois(1,1000))$X2 %>% 
#       {prcomp(.)$x[,1:2]} %>% plot(xlim=lim,ylim=lim)
#   }
#   #
#   #
#   
#   n=1000; p=4; r=2; snr=2; seed=1
#   X <- sim.simplex2(n,p,r,snr=100,d=10,d0=1,seed=rpois(1,1000))$X2
#   mean(X<1e-10)
#   png.quaternary3d(X, png.ppca(X,1))
#   
#   fit1 <- png.ppca(X)
#   fit2 <- png.ppca_qp(X)
#   fit3 <- png.gppca_qp(X)
#   png.quaternary(X, png.ppca(X))
#   png.quaternary(X, png.ppca_qp(X))
#   png.quaternary(X, png.gppca_qp(X))
#   
#   fit1$vhat
#   fit2$vhat
#   fit3$vhat
#   
#   
#   #2
#   n=100; p=100; r=2; snr=1; seed=1
#   X <- sim.simplex2(n,p,r,snr,seed+2,d=20)$X2
#   fit1 <- png.ppca(X)
#   fit2 <- png.ppca_qp(X)
#   fit3 <- png.gppca_qp(X)
#   fit1$vhat %>% head
#   fit2$vhat %>% head
#   fit3$vhat %>% head
#   
# }



#' @export calc_snr
calc_snr <- function(delta, mu,U,V, seed1){
  n=nrow(U); p=nrow(V)
  set.seed(seed1)
  E <- matrix(runif(n*p,-delta,delta),n,p) %*% (diag(p) - 1/p*matrix(1,p,p))
  X <- tcrossprod(rep(1,n),mu) + tcrossprod(U,V) + E
  sum(tcrossprod(U,V)^2) / sum(E^2)
}

#' @export calc_snr2
calc_snr2 <- function(delta, mu, X0, seed1){
  n=nrow(X0); p=ncol(X0)
  UV <- X0 - tcrossprod(rep(1,n),mu) # scale(X0, center=TRUE, scale=FALSE)
  
  set.seed(seed1)
  E <- matrix(runif(n*p,-delta,delta),n,p) %*% (diag(p) - 1/p*matrix(1,p,p))
  X <- UV + E
  sum(UV^2) / sum(E^2)
}

#' @export png.NormalVector
png.NormalVector <- function(x,v){
  x-sum(x*v)/norm(v,"2") * v/norm(v,"2")
}



#' @export sim.simplex
sim.simplex <- function(n, p, r, snr=2, d=10, d0=0, seed=1, seed.U=seed, seed.V=seed, alpha=NULL, eta=0, verbose=FALSE){
  # Simulate compositional data
  
  if(FALSE){
    snr=2; d=10; d0=0.1; seed=1; seed.U=seed; seed.V=seed; alpha=NULL
  }
  
  if(FALSE){
    n=500; p=4; r=1; d=10
    
    line=-4.5
    pdf(file="Figure-sim.simplex.pdf", width=10, height=6)
    par(mfrow=c(2,3), mar = rep(0.2, 4), mai=c(0.4,0,0,0), omi=c(0,0,0,0))
    snr=2; eta=-0.1
    sim.simplex(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    snr=2; eta=0
    sim.simplex(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    snr=2; eta=0.1
    sim.simplex(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    
    snr=1; eta=0.1
    sim.simplex(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    snr=2; eta=0.1
    sim.simplex(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    snr=5; eta=0.1
    sim.simplex(n=n, p=p, r=r, snr=snr, d=d, eta=eta, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$V, mu=.$mu, use.par=F)}
    mtext(text=sprintf("[snr, eta] = [%d, %1.1f]", snr, eta), side=3, line=line)
    
    mtext(text=sprintf("[n, p, r, d] = [%d, %d, %d, %d]", n,p,r,d), side=3, outer=TRUE, line=-2) # adj=0.02
    
    dev.off()
  }
  
  
  
  
  
  
  if(is.null(alpha)){
    alpha <- rep(10,p)
  }
  
  if(sum(alpha)==1){
    mu <- alpha
  } else {
    set.seed(seed.U)
    mu <- as.vector(dirmult::rdirichlet(1, alpha=alpha))
  }
  
  
  set.seed(seed.V)
  V <- qr.Q(qr(cbind(1, do.call("cbind", lapply(1:r, function(x) rnorm(p))) )))[,-1,drop=F]
  
  
  # D <- sapply(1:r, function(k) d/k * (k<=ceiling(r/2)) + d0)
  D <- sapply(1:r, function(k) d/k + d0)
  
  set.seed(seed.U)
  U <- matrix(0,n,r)
  for( k in 1:r ){
    v = V[,k];
    start=onedimconvexprojection(mu,-100*v,v); end=onedimconvexprojection(mu,100*v,v)
    U[,k]=truncnorm::rtruncnorm(n, a=start-eta, b=end+eta, mean=0, sd=D[k])
  }
  
  
  delta_snr <- optimize(function(x, mu=mu,U=U,V=V, seed2) norm(calc_snr(delta=x, mu=mu,U=U,V=V, seed1=seed2)-snr, "2"), c(0,1), tol=1e-10, mu=mu,U=U,V=V, seed2=seed.U)$minimum
  
  
  set.seed(seed.U)
  E <- matrix(runif(n*p,-delta_snr,delta_snr),n,p) %*% (diag(p) - 1/p*matrix(1,p,p))
  X <- tcrossprod(rep(1,n),mu) + tcrossprod(U,V) + E
  X2 <- t(apply(X, 1, png.proj2simplex))
  X0 <- (tcrossprod(rep(1,n),mu) + tcrossprod(U,V)) %>% {
    t(apply(., 1, png.proj2simplex))
  }
  
  snr_out <- sum(tcrossprod(U,V)^2) / sum(E^2)
  if(verbose) print(paste0("delta=",round(delta_snr,4), "; snr=", round(snr_out,4)))
  
  idx_not_sum_to_unit <- which(apply(X2,1,sum)>0.9999 & apply(X2,1,sum)<1)
  X2[idx_not_sum_to_unit,] <- t(apply(X2[idx_not_sum_to_unit,,drop=F],1,function(x) x/sum(x)))
  
  if(verbose) print(paste0("seed=", seed))
  
  
  params <- list(n=n, p=p, r=r, snr=snr, d=d, d0=d0, seed=seed, seed.U=seed.U, seed.V=seed.V, alpha=alpha, eta=eta)
  
  result <- list(mu=mu, U=U, D=D, V=V, E=E, X0=X0, X=X, X2=X2, params=params)
  
    
  return( result )
  
}






#' @export sim.LogNormal
sim.LogNormal <- function(n, p, r, snr=2, d=10, d0=0, seed=1, seed.U=seed, seed.V=seed, verbose=FALSE){
  # Simulate compositional data
  
  if(FALSE){
    n=500; p=4; r=1; snr=2; d=10; d0=0; seed=1; seed.U=seed; seed.V=seed; verbose=FALSE
    
    data <- sim.LogNormal(n=n, p=p, r=r, snr=5, d=10, d0=0, seed=1, verbose=TRUE)
    data %>% {
      Xnew <- png.lrpca(.$X2, nrank=2, zero.replace="simple", delta=1e-10)$Xnew
      png.quaternary3d(Xnew, vhat=.$Vlist, mu=.$mu)
    }
    
    data$X2 %>% apply(1,function(x) any(x==0)) %>% mean
    
    if(FALSE){
      n=500; p=4; r=1;
      
      line=-4.5
      pdf(file="Figure-sim.LogNormal.pdf", width=10, height=6)
      par(mfrow=c(2,3), mar = rep(0.2, 4), mai=c(0.4,0,0,0), omi=c(0,0,0,0))
      snr=2; d=5
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[log(snr), d] = [%d, %d]", snr, d), side=3, line=line)
      snr=2; d=10
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[log(snr), d] = [%d, %d]", snr, d), side=3, line=line)
      snr=2; d=100
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[log(snr), d] = [%d, %d]", snr, d), side=3, line=line)
      
      snr=1; d=10
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[log(snr), d] = [%d, %d]", snr, d), side=3, line=line)
      snr=2; d=10
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[log(snr), d] = [%d, %d]", snr, d), side=3, line=line)
      snr=5; d=10
      sim.LogNormal(n=n, p=p, r=r, snr=snr, d=d, d0=0, seed=1, verbose=TRUE) %>% {png.quaternary(X=.$X2+1e-16, vhat=.$Vlist, mu=.$mu, use.par=F)}
      mtext(text=sprintf("[log(snr), d] = [%d, %d]", snr, d), side=3, line=line)
      
      mtext(text=sprintf("[n, p, r] = [%d, %d, %d]", n,p,r), side=3, outer=TRUE, line=-2) # adj=0.02
      
      dev.off()
    }
    
    
  }
  
  snr <- exp(snr)
  
  mu <- rep(0,p)
  
  set.seed(seed.V)
  V <- qr.Q(qr( do.call("cbind", lapply(1:r, function(x) rnorm(p))) ))
  

  D <- sapply(1:r, function(k) log(d/k + d0))
  
  set.seed(seed.U)
  U <- matrix(0,n,r)
  for( k in 1:r ){
    U[,k]=rnorm(n, mean=0, sd=D[k])
  }
  
  UZ <- matrix(0,n,r)
  for(i in 1:n){
    for(k in 1:r){
      UZ[i,k] <- rbinom(1, size=1, prob=exp(-0.01*U^2)[i,k])
    }
  }
  
  U[UZ==0] <- 20*sign(U[UZ==0])
  print(mean(UZ==0))
  X0 <- tcrossprod(rep(1,n),mu) + tcrossprod(U,V)
  X02 <- t(apply(X0,1,function(x) exp(x)/sum(exp(x))))
  #
  mu0 <- png.iclr(mu)
  delta_snr <- optimize(function(x, mu=mu,U=U,V=V, seed2) norm(calc_snr2(delta=x, mu=mu0, X0=X02, seed1=seed2)-snr, "2"), c(0,1), tol=1e-10, mu=mu,U=U,V=V, seed2=seed.U)$minimum
  #
  set.seed(seed.U)
  E <- matrix(runif(n*p,-delta_snr,delta_snr),n,p) %*% (diag(p) - 1/p*matrix(1,p,p))
  X <- X02 + E
  
  X2 <- t(apply(X, 1, png.proj2simplex))
  #
  snr_out <- sum((X02 - tcrossprod(rep(1,n),mu0) )^2) / sum(E^2)
  if(verbose) print(paste0("delta=",round(delta_snr,4), "; snr=", round(snr_out,4), "; log(snr)=", round(log(snr_out),4) ))
  #
  # X02 %>% png.quaternary3d()
  # X2 %>% png.quaternary3d()
  
  # X2 %>% {
  #   fit1 <- png.lrpca(., zero.replace="simple", nrank=2, delta=1e-10)
  #   fit2 <- png.lrpca(.[!apply(., 1,function(x)any(x>0.8)),], zero.replace="simple", nrank=2, delta=1e-10)
  #   png.quaternary3d(X=fit1$Xnew,
  #                    vhat=fit2$vhat,
  #                    mu=fit2$mu)
  # }
  
  
  idx_not_sum_to_unit <- which(apply(X2,1,sum)>0.9999 & apply(X2,1,sum)<1)
  X2[idx_not_sum_to_unit,] <- t(apply(X2[idx_not_sum_to_unit,,drop=F],1,function(x) x/sum(x)))
  
  if(verbose) print(paste0("seed=", seed))
  
  
  params <- list(n=n, p=p, r=r, snr=snr, d=d, d0=d0, seed=seed, seed.U=seed.U, seed.V=seed.V)
  
  # mu <- png.iclr(mu)
  
  result <- list(mu=mu, U=U, D=D, V=V, Vlist=lapply(seq(-20,20,1), function(z) z*V), E=E, X=X, X2=X2, params=params)
  
  
  return( result )
  
}






#' @export sim.simplex.test
sim.simplex.test <- function(params){
  TestData <- params %>% 
    { sim.simplex(n=.[["n"]], 
                  p=.[["p"]], 
                  r=.[["r"]], 
                  snr=.[["snr"]], 
                  d=.[["d"]],
                  d0=.[["d0"]],
                  seed.U=.[["seed.U"]]*123,
                  seed.V=.[["seed.V"]],
                  eta=.[["eta"]]) }
  
  TestData
}





#' @export sim.simplex2
sim.simplex2 <- function(n, p, r, snr=2, d=10, d0=0.1, seed=1, seed.U=seed, seed.V=seed, alpha=NULL){
  # Simulate compositional data

  if(FALSE){
    snr=100; d=10; d0=1; seed=1; seed.U=seed; seed.V=seed; alpha=NULL
  }
  
  
  if(is.null(alpha)){
    alpha <- rep(10,p)
  }
  
  mu <- as.vector(dirmult::rdirichlet(1, alpha=alpha))
  
  set.seed(seed.V)
  V <- qr.Q(qr(cbind(1, do.call("cbind", lapply(1:r, function(x) rnorm(p))) )))[,-1,drop=F]
  
  
  D <- sapply(1:r, function(k) d/k * (k<=ceiling(r/2)) + d0)
  
  set.seed(seed.U)
  U <- matrix(0,n,r)
  for( k in 1:r ){
    v = V[,k];
    start=onedimconvexprojection(mu,-100*v,v); end=onedimconvexprojection(mu,100*v,v)
    U[,k]=truncnorm::rtruncnorm(n, a=start, b=end, mean=0, sd=0.025)
  }
  
  U <- U %*% diag(D,r,r)
  
  delta_snr <- optimize(function(x, mu=mu,U=U,V=V, seed2) sum(calc_snr(delta=x, mu=mu,U=U,V=V, seed1=seed2)-snr)^2, c(0,1), tol=1e-10, mu=mu,U=U,V=V, seed2=seed.U)$minimum
  
  
  set.seed(seed.U)
  E <- matrix(runif(n*p,-delta_snr,delta_snr),n,p) %*% (diag(p) - 1/p*matrix(1,p,p))
  X <- tcrossprod(rep(1,n),mu) + tcrossprod(U,V) + E
  
  {
    X2 <- X
    for( k in 1:r ){
      X2 <- t(apply(X2,1,function(x) png.NormalVector(x-mu,V[,k])+mu+tcrossprod(onedimconvexprojection(mu,x,V[,k]),V[,k]) ))
    }
    
    # {
    #   X2 <- X
    #   X2 <- t(apply(X2,1,function(x) mu+tcrossprod(multidimconvexprojection(mu,x,V),V) ))
    #   
    #   png.quaternary3d(X2+E)
    #   png.quaternary3d(tcrossprod(rep(1,n),mu) + tcrossprod(U,V) )
    # }
    
    # png.quaternary(X2, fit=list(mu=mu,vhat=V[,1,drop=F]))
    X2 <- t(apply(X, 1, png.proj2simplex))
    X2 <- ifelse(abs(X2)<1e-15,0,X2)
    # png.quaternary(X2, fit=list(mu=mu,vhat=V[,1,drop=F]))
  }
  
  
  snr_out <- sum(tcrossprod(U,V)^2) / sum(E^2)
  if(verbose) print(paste0("delta=",round(delta_snr,4), "; snr=", round(snr_out,4)))
  
  idx_not_sum_to_unit <- which(apply(X2,1,sum)>0.99909 & apply(X2,1,sum)<1)
  X2[idx_not_sum_to_unit,] <- t(apply(X2[idx_not_sum_to_unit,,drop=F],1,function(x) x/sum(x)))
 
  if(verbose) print(paste0("seed=", seed))
  return( list(mu=mu, U=U, D=D, V=V, E=E, X=X, X2=X2) )
  
}



