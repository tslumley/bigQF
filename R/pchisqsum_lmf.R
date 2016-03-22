pchisqsum_lmf <-
function(x,Mmult,tMmult,ncolM,nrowM, traceM, n=100, tr2.sample.size=300, method=c("saddlepoint","integration")){
	method<-match.arg(method)
	extM<-extmat(function(x) as.numeric(Mmult(x)), function(x) as.numeric(tMmult(x)), nrow=nrowM, ncol= ncolM)
	ee<-svd::trlan.svd(extM, neig=n)$d^2
	tr <- traceM
	tr2<-trace_mf(Mmult,tMmult,k=tr2.sample.size,ncolM,trace.full=traceM)
	tr.small<-tr-sum(ee)
	tr2.small<-tr2-sum(ee^2)
	scale<-tr2.small/tr.small
	nu<-(tr.small^2)/tr2.small
    pchisqsum(x, c(rep(1,n), ceiling(nu)), c(ee, scale), method=method,lower.tail=FALSE)
}
