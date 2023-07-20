#' @export png.ternary.init
png.ternary.init <- function(X, cex=0.8, color="grey50"){
  # par(mfrow = c(1,1), mar = rep(0.5, 4))
  par(mar = rep(0.5, 4))
  TernaryPlot("A","B","C", grid.lines = 5, grid.lty = "dotted",
              grid.minor.lines = 1, grid.minor.lty = "dotted")
  AddToTernary(text, X, "o", cex = cex, font = 2, col=color)
}

#' @export png.ternary.PC_axis
png.ternary.PC_axis <- function(mu, vhat, color="grey20", cex=1){
  
  # fit.StartEnd <- png.loading2StartEnd(mu, vhat)
  # for( k in 1:length(Start) ){
  #   TernaryLines(fit.StartEnd[[k]], col = color)
  #   TernaryText(fit.StartEnd[[k]][2,], paste0("PC",k), cex = cex, col = "darkblue", font = 2)
  # }
  
  Start <- apply(vhat,2,function(v)onedimconvexprojection(mu, v*-10, v))
  End <- apply(vhat,2,function(v)onedimconvexprojection(mu, v*10, v))
  for( k in 1:length(Start) ){
    TernaryLines(rbind((mu+Start[k]*vhat[,k]), (mu+End[k]*vhat[,k])), col = color)
    TernaryText((mu+End[k]*vhat[,k]), paste0("PC",k), cex = cex, col = "darkblue", font = 2)
  }
  
  return(list(Start=Start, End=End))
}

#' @export png.ternary.point
png.ternary.point <- function(tuple, shape="x", cex, color="red"){
  AddToTernary(text, tuple, shape, cex = cex, font = 1, col=color)
}

#' @export png.ternary.ternary2coord
png.ternary.ternary2coord <- function(x){
  Ternary::TernaryToXY(x)
  # a=x[1]; b=x[2]; c=x[3]
  # c(1/2*(2*b+c)/(a+b+c), sqrt(3)/2*(c)/(a+b+c))
}

#' @export png.ternary.coord2ternary
png.ternary.coord2ternary <- function(x){
  Ternary::XYToTernary(x[1], x[2])
  # a=x[1]; b=x[2]; c=x[3]
  # c(1/2*(2*b+c)/(a+b+c), sqrt(3)/2*(c)/(a+b+c))
}

#' @export png.rotmat
png.rotmat <- function(angle){
  theta = angle * 180/pi
  matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
}



#' @export png.ternary
png.ternary <- function(X, vhat=NULL, xhat=NULL, cex=1.0){
  if(FALSE){
    cex=0.5
  }
  
  if(FALSE){
    n=150; p=3; r=2;
    seed <- 2
    X <- sim.simplex2(n=n,p=p,r=r,snr=5, seed=1)$X2 #1, 5, 10, 17
    fit <- png.gppca_qp(X, nrank=2, epsilon=1e-8)
    fit <- png.lrpca(png.ZeroReplace.additive(X), nrank=2)
    
    png.ternary(X, vhat=fit$vhat)
    png.ternary(X, fit=fit, print.xhat=TRUE)
    png.plot.pca(fit,cex=0.8)
    
    png.quaternary(X, fit=fit)
    png.quaternary(fit$xhat, fit=fit)
    
    png.quaternary3d(X, fit=fit)
    png.quaternary3d(fit$xhat, fit=fit)
  }
  
  mu=colMeans(X)
  
  png.ternary.init(X, color="grey40", cex=cex*0.6)
  
  if(!is.null(vhat)){
    
    if(is.list(vhat)){
      
      for( jj in 1:NCOL(vhat[[1]]) ){
        for( ii in 1:(length(vhat)-1) ){
          vhat.i1 = vhat[[ii]][,jj]
          vhat.i2 = vhat[[ii+1]][,jj]
          # vhat.start <- png.iclr( vhat.i1 )
          # vhat.end <- png.iclr( vhat.i2 )
          # vhat.start <- png.iclr( log(mu)+vhat.i1 )
          # vhat.end <- png.iclr( log(mu)+vhat.i2 )
          vhat.start <- png.iclr( colMeans(log(X))+vhat.i1 )
          vhat.end <- png.iclr( colMeans(log(X))+vhat.i2 )
          
          
          TernaryLines(rbind(vhat.start,vhat.end), col = "darkblue")
        }
        TernaryText(vhat.end, paste0("PC",jj), cex = cex*0.8, col = "darkblue", font = 2)
      }
      
    } else {
      
      png.ternary.PC_axis(mu, vhat, color="darkblue", cex=cex*0.8)
      
    }
  }
  
  if(!is.null(xhat)){
    AddToTernary(text, xhat, "x", font = 1.5, col="darkred", cex=cex*0.6)
    text(x = 0, y = -0.15, "o: data points;  x: fitted values", col = "black", cex=cex*0.8)
  }
  
  
}
