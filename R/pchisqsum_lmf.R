pchisqsum_lmf <-
function(x,Mmult,tMmult,ncolM,nrowM, traceM, n=100, tr2.sample.size=300, method=c("saddlepoint","integration"),remainder="warn"){
	method<-match.arg(method)
	extM<-extmat(function(x) as.numeric(Mmult(x)), function(x) as.numeric(tMmult(x)), nrow=nrowM, ncol= ncolM)
	ee<-svd::trlan.svd(extM, neig=n)$d^2

	## estimates tr2 always, tr if trace.full is NULL
	traces <-trace_mf(Mmult,tMmult,k=tr2.sample.size,ncolM,trace.full=traceM)
	tr2 <- traces[1]
	tr <- traces[2]

	tr.small<-tr-sum(ee)
	tr2.small<-tr2-sum(ee^2)
        if (tr.small==0) {
            scale<-0
            nu<-1
        } else {
            scale<-tr2.small/tr.small
            nu<-(tr.small^2)/tr2.small
        }
        pchisqsum(x, c(rep(1,n), nu), c(ee, scale), method=method,lower.tail=FALSE,remainder=remainder)
}
