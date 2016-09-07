pchisqsum_smf <-
function(x,Mmult,tMmult,ncolM,traceM, n=100,p=10,q=2, tr2.sample.size=300,
         method=c("saddlepoint","integration"),remainder="warn"){
    method <- match.arg(method)
    Q <- matfreeQ(Mmult, tMmult, n + p, ncolM, q = q)
    B <- tMmult(Q)
	
    ## singular values from projected matrix
    s <- svd(B,nu=0,nv=n)
    ee<-s$d[1:n]^2
	
    ## Q1 spans the range of projected matrix
    Q1<-Q%*%s$v[,1:n]
    k<-tr2.sample.size; trace.full=traceM
    
    Omega <- matrix(rnorm(k * ncolM), ncol = k)
    AOmega0<-Mmult(Omega)
    ## project AOmega0 orthogonal to range of projected matrix
    
    AOmega <- AOmega0-Q1%*%crossprod(Q1,AOmega0)
    AAOmega <- tMmult(AOmega)
    tr <- sum(rowSums(AOmega * AOmega))/k
    trsquared <- sum(rowSums(AAOmega * AAOmega))/k
    
    if (is.null(trace.full) || (sum(ee)>trace.full)){
    	tr.small<-tr
    	tr2.small<-trsquared    	
    } else {
        ## subtract from known trace
        tr.small<-trace.full-sum(ee)
        ## ratio estimator
        tr2.small = trsquared * (tr.small/tr)^2
    }
    
    scale <- tr2.small/tr.small
    nu <- (tr.small^2)/tr2.small
    pchisqsum(x, c(rep(1, n), nu), c(ee, scale), method = method, 
              lower.tail = FALSE, remainder = remainder)
}
