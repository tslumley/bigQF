pchisqsum_ssvd <-
function(x,M,n=100,p=10,q=0, method=c("saddlepoint","integration"),remainder="warn"){
	method<-match.arg(method)
	Q<-srfht(M,n+p,q=q)
	B<-t(Q)%*%M%*%Q
	ee<-eigen(B,symmetric=TRUE,only.values=TRUE)$values[1:n]
	tr<-sum(diag(M))
	tr2<-sum(M^2)
	tr.small<-tr-sum(ee)
	tr2.small<-tr2-sum(ee^2)
	scale<-tr2.small/tr.small
	nu<-(tr.small^2)/tr2.small
    pchisqsum(x, c(rep(1,n),ceiling(nu)), c(ee, scale), method=method,lower.tail=FALSE,remainder=remainder)
}
