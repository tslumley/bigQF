pchisqsum_rsvd <-
function(x,M,n=100,p=10,q=2, tr2.sample.size=100, method=c("saddlepoint","integration"),remainder="warn"){
	method<-match.arg(method)
	Q<-srfht2(M,n+p,q=q)
	B<-M%*%Q
	ee<-svd(B,nu=0,nv=0)$d[1:n]^2
	diags <- colSums(M^2)
	tr<-sum(diags)
	if (tr2.sample.size>0){
		tr2<-tracefht(M,k=tr2.sample.size,trace.full=tr)
	} else {
		Ms<-crossprod(M)
		tr2<- sum(Ms^2)
	}	
	tr.small<-tr-sum(ee)
	tr2.small<-tr2-sum(ee^2)
	scale<-tr2.small/tr.small
	nu<-(tr.small^2)/tr2.small
    pchisqsum(x, c(rep(1,n), nu), c(ee, scale), method=method,lower.tail=FALSE,remainder=remainder)
}
