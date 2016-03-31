
is.square<-function(M) nrow(M)==ncol(M)

pQF<-function(x, M, method=c("ssvd","lanczos","satterthwaite"), neig=100, tr2.sample.size=500, q=NULL,  convolution.method=c("saddlepoint","integration")){
   method<-match.arg(method)
   conv.method<-match.arg(convolution.method)

## sanity check
   if(method != "satterthwaite"){
     mindim<-min(dim(M))
     if (mindim < neig) stop("Can't have more eigenvalues than min(nrow,ncol)")
   }
  
## sparse matrix
   if (inherits(M, "matrixfree")){
     if (is.null(q)) q<- 3
     switch(method, ssvd=pchisqsum_smf(x,M$mult, M$tmult, M$ncol, M$trace, n=neig, q=q, tr2.sample.size=tr2.sample.size, method=conv.method),
     lanczos= pchisqsum_lmf(x,M$mult, M$tmult, M$ncol,M$nrow, M$trace, n=neig,  tr2.sample.size=tr2.sample.size, method=conv.method),
     satterthwaite= pchisqsum_spsatt(x,M$mult, M$tmult, M$ncol, M$trace,tr2.sample.size=tr2.sample.size))

   } else { ## square matrix H=G^TG
   if (is.square(M)){
      if (is.null(q)) q<- 1
      switch(method, ssvd=pchisqsum_ssvd(x,M, n=neig, q=q,  method=conv.method),
      lanczos= pchisqsum_partial(x,M, n=neig,  method=conv.method),
      satterthwaite= pchisqsum_partial(x,M, n=0,  method=conv.method))
   } else { 
   if (is.null(q)) q<- 3
      switch(method, ssvd=pchisqsum_rsvd(x,M, n=neig, q=q,  method=conv.method, tr2.sample.size=tr2.sample.size),
      lanczos= pchisqsum_rpartial(x,M, n=neig,  method=conv.method, tr2.sample.size=tr2.sample.size),
      satterthwaite= pchisqsum_Gsatt(x,M,  tr2.sample.size=tr2.sample.size))
   }
   }
}

print.matrixfree<-function(x,...){cat("implicit matrix object:",x$nrow,"x",x$ncol,"\n")}

dim.matrixfree<-function(x){c(x$nrow,x$ncol)}


sparse.matrixfree<-function(M){

  rval<-list(
    mult=function(X) M%*%X,
	
    tmult=function(X) crossprod(M,X),	
    trace=sum(M^2),
    ncol=ncol(M),
    nrow=nrow(M)
  )
 class(rval)<-"matrixfree"
 rval
}



SKAT.matrixfree.default<-function(G,weights=function(maf) dbeta(maf,1,25),model=NULL){
  center<-colMeans(G)
  ww<-weights(center/2)
  spG<-Matrix(G,sparse=TRUE)%*%Diagonal(x=ww)
  cntr<-colMeans(spG)

  rval<-list(
    mult=function(X){
	 t(t(spG%*%X)-colSums(cntr*as.matrix(X)))/sqrt(2)
	},
	
    tmult=function(X){
	(crossprod(spG,X)- outer(cntr,colSums(as.matrix(X))))/sqrt(2)
	}	,	
    trace=(sum(spG^2)-sum(cntr^2)*nrow(G))/2,
    ncol=ncol(G),
    nrow=nrow(G)
  )
 class(rval)<-"matrixfree"
 rval
}

SKAT.matrixfree<-function(G,weights=function(maf) dbeta(maf,1,25), model=NULL){
  UseMethod("SKAT.matrixfree", model)
}

SKAT.matrixfree.lm<-function(G,weights=function(maf) dbeta(maf,1,25), model=NULL){
  center<-colMeans(G)
  ww<-weights(center/2)
  spG<-Matrix(G,sparse=TRUE)%*%Diagonal(x=ww)
  qr<-model$qr


rval<-list(
     mult=function(X){
	base::qr.resid(qr, as.matrix(spG%*%X))/sqrt(2)
     },
     tmult=function(X){
         crossprod(spG,qr.resid(qr,X))/sqrt(2)
      },   
    trace=NULL,
    ncol=ncol(G),
    nrow=nrow(G)
  )
  class(rval)<-"matrixfree"
  rval
}

SKAT.matrixfree.glm<-function(G,weights=function(maf) dbeta(maf,1,25), model=NULL){
  center<-colMeans(G)
  ww<-weights(center/2)
  fw <- model$family$variance(fitted(model))
  spG<-Diagonal(x=fw)%*%Matrix(G,sparse=TRUE)%*%Diagonal(x=ww)
  qr<-model$qr

rval<-list(
     mult=function(X){
	base::qr.resid(qr, as.matrix(spG%*%X))/sqrt(2)
     },
     tmult=function(X){
         crossprod(spG,qr.resid(qr,X))/sqrt(2)
      },   
    trace=NULL,
    ncol=ncol(G),
    nrow=nrow(G)
  )
  class(rval)<-"matrixfree"
  rval
}







pchisqsum<- function (x, df, a, lower.tail = FALSE, method=c("saddlepoint","integration","satterthwaite")) 
{
    satterthwaite <- function(a, df) {
        if (any(df > 1)) {
            a <- rep(a, df)
        }
        tr <- mean(a)
        tr2 <- mean(a^2)/(tr^2)
        list(scale = tr * tr2, df = length(a)/tr2)
    }
    method <- match.arg(method)
    
    ## can happen with randomised trace estimator if most remaining singular values are very small.
    ##
    if (any(bad.df <- (df<1))){
      warning("Negative/fractional df removed")
      df[bad.df]<-1
      a[bad.df]<-0
    }
    df<-round(df)
    ##
   
   
    sat <- satterthwaite(a, df)
    guess <- pchisq(x/sat$scale, sat$df, lower.tail = lower.tail)
    if (method == "satterthwaite") 
        return(guess)
    method <- match.arg(method)
    if (method == "integration" && !suppressWarnings(requireNamespace("CompQuadForm"))) {
        warning("Package 'CompQuadForm' not found, using saddlepoint approximation")
        method <- "saddlepoint"
    }
    abstol <- guess/1000
    abstol <- pmax(1e-09, abstol)
    reltol <- rep(1/1000, length(abstol))
    if (method == "integration") {
            for (i in seq(length = length(x))) {
                f <- davies(x[i], a, df, acc = 1e-09)
                if (f$ifault > 0) 
                  warning("Probable loss of accuracy ")
                guess[i] <- f$Qq
            }   
        if (lower.tail) 
            guess <- 1 - guess
        return(guess)
    }
    if (method == "saddlepoint") {
        for (i in seq(length = length(x))) {
            lambda <- rep(a, df)
            sad <- sapply(x, saddle, lambda = lambda)
            if (lower.tail) 
                sad <- 1 - sad
            guess <- ifelse(is.na(sad), guess, sad)
        }
        return(guess)
    }
}

saddle<-function (x, lambda) 
{
    d <- max(lambda)
    lambda <- lambda/d
    x <- x/d
    k0 <- function(zeta) -sum(log(1 - 2 * zeta * lambda))/2
    kprime0 <- function(zeta) sapply(zeta, function(zz) sum(lambda/(1 - 
        2 * zz * lambda)))
    kpprime0 <- function(zeta) 2 * sum(lambda^2/(1 - 2 * zeta * 
        lambda)^2)
    n <- length(lambda)
    if (any(lambda < 0)) {
        lmin <- max(1/(2 * lambda[lambda < 0])) * 0.99999
    }
    else if (x > sum(lambda)) {
        lmin <- -0.01
    }
    else {
        lmin <- -length(lambda)/(2 * x)
    }
    lmax <- min(1/(2 * lambda[lambda > 0])) * 0.99999
    hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, 
        upper = lmax, tol = 1e-08)$root
    w <- sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
    v <- hatzeta * sqrt(kpprime0(hatzeta))
    if (abs(hatzeta) < 1e-04) 
        NA
    else pnorm(w + log(v/w)/w, lower.tail = FALSE)
}
