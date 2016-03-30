pchisqsum_smf <-
function(x,Mmult,tMmult,ncolM,traceM, n=100,p=10,q=2, tr2.sample.size=300, method=c("saddlepoint","integration")){
	method<-match.arg(method)
	Q<-matfreeQ(Mmult,tMmult, n+p,ncolM,q=q)
	B<-tMmult(Q)
	ee<-svd(B,nu=0,nv=0)$d[1:n]^2

	## estimates tr2 always, tr if trace.full is NULL
	traces <-trace_mf(Mmult,tMmult,k=tr2.sample.size,ncolM,trace.full=traceM)
	tr2 <- traces[1]
	tr <- traces[2]
	
	tr.small<-tr-sum(ee)
	tr2.small<-tr2-sum(ee^2)
	scale<-tr2.small/tr.small
	nu<-(tr.small^2)/tr2.small
    pchisqsum(x, c(rep(1,n), ceiling(nu)), c(ee, scale), method=method,lower.tail=FALSE)
}
