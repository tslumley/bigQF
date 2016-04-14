pchisqsum_partial <-
function(x,M,n=100,method=c("saddlepoint","integration"),tol=1e-3,remainder="warn"){
	method<-match.arg(method)
	if (n>0){
		ee<-svd::trlan.eigen(M,n,opts=list(tol=tol))
	} else {
		ee<-list(values=numeric(0))
	}
	tr<-sum(diag(M))
	tr2<-sum(M^2)
	tr.small<-tr-sum(ee$d)
	tr2.small<-tr2-sum(ee$d^2)
	scale<-tr2.small/tr.small
	nu<-(tr.small^2)/tr2.small
	if (n>0)
	    pchisqsum(x, c(rep(1,n),ceiling(nu)), c(ee$d, scale), method=method,lower.tail=FALSE,remainder=remainder)
	else
	    pchisq(x/scale, nu,lower.tail=FALSE)
}
