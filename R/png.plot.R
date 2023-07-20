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

