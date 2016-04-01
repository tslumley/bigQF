trace_mf <-
function(Mmult,tMmult,k,cols,trace.full=NULL){
        Omega<-matrix(rnorm(k*cols), ncol=k)
	AOmega<-Mmult(Omega)
	AAOmega<-tMmult(AOmega)
	
	tr<-sum(rowSums(AOmega*AOmega))/k
	trsquared<-sum(rowSums(AAOmega*AAOmega))/k
	
	if (is.null(trace.full))
		c(tr2=trsquared, tr=tr)
	else
		c(tr2=trsquared*(trace.full/tr)^2, tr=trace.full)
}
