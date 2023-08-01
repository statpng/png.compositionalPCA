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
    n=50; p=50; r=5; snr=2; eta=0.1/log(p); seed=1
    data <- sim.simplex(n=n,p=p,r=r,snr=snr,d=10,d0=0.01,seed=seed,eta=eta)
    
    fit1 <- png.ppca_qp(data$X2, nrank=r, kappa=1e-6, maxit=2000, eps=1e-6, gamma=0.5, save.est.path = TRUE)
    fit2 <- png.gppca_qp(data$X2, nrank=r, kappa=1e-6, maxit=2000, eps=1e-6, gamma=0.5, save.est.path = TRUE)
    
    fit1 %>% png.pca.criteria(data, n.test=n*5)
    fit2 %>% png.pca.criteria(data, n.test=n*5)
    
  }
  
  if(FALSE){
    fit=fit_ppca
    data
    n.test=200
  }
  

  data$params["n"] <- n.test
  true <- data %>% { list(Xtrain=.$X2,
                          Xtest=sim.simplex.test(.$params)$X2,
                          V=.$V) }
  
  Xtrain=true$Xtrain
  Xtest=true$Xtest
  Vtrue=true$V
  
  vhat <- fit$vhat
  xhat_train <- fit$xhat
  xhat_test <- png.projection(Xtest, fit, method=fit$method)
  
  # png.projection(Xtrain, fit, method=fit$method)[1:5,1:5]
  # xhat_train[1:5,1:5]
  
  obj.Xtrain <- sqrt(mean((Xtrain-xhat_train)^2))
  obj.Xtest <- sqrt(mean((Xtest-xhat_test)^2))
  
  rmse.Xtrain <- sqrt(mean((Xtrain-xhat_train)^2))
  rmse.Xtest <- sqrt(mean((Xtest-xhat_test)^2))
  Pangle.V <- png.angle(Vtrue, vhat)$max
  Gangle.V <- png.angle(Vtrue, vhat)$Grassmannian
  # Out-of-simplex Sample Percentage
  OutOfSimplex <- mean(apply(xhat_train,1,function(x) any(x < -1e-8)))
  Sparsity <- mean( abs(xhat_train) < 1e-12 )
  Convergence <- sapply(fit$fit.path, function(fit_rank) fit_rank$it < fit$maxit)
  # fit$fit.path[[4]]$crit.path[-(1:100)] %>% plot(type="l")
  # fit$fit.path[[4]]$crit.path %>% tail
  
  
  c(rmse.Xtrain=rmse.Xtrain,
    rmse.Xtest=rmse.Xtest,
    Pangle.V=Pangle.V,
    Gangle.V=Gangle.V,
    OutOfSimplex=OutOfSimplex,
    Sparsity=Sparsity,
    Convergence=Convergence)
  
}





#' @export png.CompositionalPlot
png.CompositionalPlot <- function(pseq, xhat=NULL, title="", legend.ncol=10){
  if(FALSE){
    pseq <- pseq_list_total_xhat$urine$Phylum
    xhat <- fit4$xhat
  }
  
  if(!is.null(xhat)){
    pseq@otu_table <- t(xhat) %>% otu_table(taxa_are_rows = TRUE )
  }
  rownames(pseq@otu_table) <- pseq@tax_table[,"unique"]
  
  n.taxa <- length(taxa(pseq))
  # otu.sort <- "abundance"
  otu.sort <- top_taxa(pseq, n=n.taxa)
  # otu.sort <- rev(c(rev(names(sort(rowSums(abu)))[seq(1, nrow(abu), 2)]), 
  # names(sort(rowSums(abu)))[seq(2, nrow(abu), 2)]))
  
  TopTaxa <- pseq@tax_table[ top_taxa(pseq, n=1), "unique" ]
  sample.sort <- rev(sample_names(pseq)[order(abundances(pseq)[TopTaxa, 
  ])])
  
  
  library(RColorBrewer)
  # set.seed(1)
  Taxa_cols = brewer.pal.info[brewer.pal.info$category == 'qual',] %>% { unlist(mapply(brewer.pal, .$maxcolors, rownames(.))) }# %>% sample()
  
  
  abu <- abundances(pseq)
  dfm <- psmelt(otu_table(abu, taxa_are_rows = TRUE))
  names(dfm) <- c("Tax", "Sample", "Abundance")
  
  dfm$Sample <- factor(dfm$Sample, levels=sample.sort)
  dfm$Tax <- factor(dfm$Tax, levels=otu.sort)
  
  {
    Tax <- Sample <- Abundance <- NULL
    dfm <- dfm %>% arrange(Tax)
    dfm$Tax <- factor(dfm$Tax, levels = unique(dfm$Tax))
    p <- ggplot(dfm, aes(x=Sample, y=Abundance, fill=Tax)) + 
      geom_bar(position="stack", stat="identity") + 
      scale_x_discrete(labels = dfm$xlabel, 
                       breaks = dfm$Sample)
    p <- p + labs(y = "Abundance")
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                              hjust = 0))
    p <- p + guides(fill = guide_legend(reverse = FALSE))
    # if (!is.null(group_by)) {
    #   p <- p + facet_grid(. ~ Group, drop = TRUE, space = "free", 
    #                       scales = "free")
    # }
    p
  }
  
  
  
  p +
    # scale_fill_manual("", breaks=taxa, values = Taxa_cols, na.value = "black") +
    scale_fill_manual("", values = Taxa_cols, na.value = "black") +
    scale_y_continuous(label = scales::percent) +
    # hrbrthemes::theme_ipsum(grid="Y") +
    theme_bw(base_size=14) +
    labs(x = "Samples", y = "Relative Abundance",
         subtitle = title) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "bottom") +
    guides(fill = guide_legend(ncol = legend.ncol))
  
}




#' @export png.MultiCompositionalPlot
png.MultiCompositionalPlot <- function(p.list, title="", legend.ncol=10){
  if(FALSE){
    p.list <- list(ppca=p1, gppca=p2, lrpca=p3)
  }
  
  
  p.names <- names(p.list)
  data.list <- lapply(p.list, function(x) x$data)
  
  
  TopTaxa.list <- lapply(data.list, function(x) x %>% group_by(Tax) %>% summarise(avg=mean(Abundance)) %>% arrange(desc(avg)) %>% slice_head(n=1) %>% .$Tax)
  
  
  OrderedSample.list <- lapply(1:length(data.list), function(i) data.list[[i]] %>% filter(Tax == TopTaxa.list[[i]]) %>% arrange(desc(Abundance)) %>% .$Sample)
  
  tmp <- gtools::mixedsort(as.character(unique(OrderedSample.list[[1]])))
  sample.order <- lapply(OrderedSample.list, function(y){
    purrr::map_int( tmp, ~which( y %in% .x ) )
  }) %>% {Reduce("+", .)/length(.)} %>% order
  #
  
  for( i in 1:length(data.list) ){
    data.list[[i]]$Sample <- factor(data.list[[i]]$Sample, levels=tmp[sample.order])
  }
  
  
  dfm <- NULL
  for(i in 1:length(data.list) ){
    dfm <- rbind.data.frame(dfm, cbind.data.frame(Type=p.names[i], data.list[[i]]))
  }
  
  
  n.taxa <- length(unique(p.list[[1]]$data$Tax))
  TopTaxa <- dfm %>% group_by(Tax) %>% summarise(avg=mean(Abundance)) %>% slice_head(n=n.taxa) %>% .$Tax
  otu.sort <- TopTaxa
  
  
  library(RColorBrewer)
  # set.seed(1)
  Taxa_cols = brewer.pal.info[brewer.pal.info$category == 'qual',] %>% { unlist(mapply(brewer.pal, .$maxcolors, rownames(.))) }# %>% sample()
  
  
  dfm$Type <- factor(dfm$Type, levels=p.names)
  dfm$Sample <- factor(dfm$Sample, levels=sample.sort)
  dfm$Tax <- factor(dfm$Tax, levels=otu.sort)
  
  {
    Tax <- Sample <- Abundance <- NULL
    dfm <- dfm %>% arrange(Tax)
    dfm$Tax <- factor(dfm$Tax, levels = unique(dfm$Tax))
    p <- ggplot(dfm, aes(x=Sample, y=Abundance, fill=Tax)) + 
      geom_bar(position="stack", stat="identity") + 
      scale_x_discrete(labels = dfm$xlabel, 
                       breaks = dfm$Sample)
    p <- p + labs(y = "Abundance")
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                              hjust = 0))
    p <- p + guides(fill = guide_legend(reverse = FALSE))
    p <- p + facet_grid( Type ~ ., drop = TRUE, space = "free",
                        scales = "free")
    
    p +
      # scale_fill_manual("", breaks=taxa, values = Taxa_cols, na.value = "black") +
      scale_fill_manual("", values = Taxa_cols, na.value = "black") +
      scale_y_continuous(label = scales::percent) +
      # hrbrthemes::theme_ipsum(grid="Y") +
      theme_bw(base_size=14) +
      labs(x = "Samples", y = "Relative Abundance",
           subtitle = title) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.position = "bottom") +
      guides(fill = guide_legend(ncol = legend.ncol))
  }
  
  
}

