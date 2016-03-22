tracefht <-
function(A,k,trace.full=NULL){
	m<-NROW(A)
	mbig<-2^ceiling(log2(m))
	R<-sample(c(-1,1),m,replace=TRUE)
	Astar<-A*R
	idx<-sample(mbig,k)
	AOmega<-fht(Astar)[idx,]
	AAOmega<-tcrossprod(AOmega,A)
	if (!is.null(trace.full))
		tr<-sum(rowSums(AOmega*AOmega))/k
	trsquared<-sum(rowSums(AAOmega*AAOmega))/k
	if (is.null(trace.full))
		trsquared
	else
		trsquared*(trace.full/tr)^2 
}
