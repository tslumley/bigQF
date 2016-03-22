pchisqsum_Gsatt <-
function (x, M, tr2.sample.size=100){
	  tr<-sum(M^2)
	  tr2<-tracefht(M,k=tr2.sample.size,trace.full=tr)
	  scale<-tr2/tr
	  nu<-(tr^2)/tr2
      pchisq(x/scale, nu,lower.tail=FALSE)
	}
	
pchisqsum_spsatt <-
function (x, Mmult,tMmult,ncolM,traceM, tr2.sample.size=100){
	tr <- traceM
	tr2<-trace_mf(Mmult,tMmult,k=tr2.sample.size,ncolM,trace.full=traceM)
	  scale<-tr2/tr
	  nu<-(tr^2)/tr2
      pchisq(x/scale, nu,lower.tail=FALSE)
	}
