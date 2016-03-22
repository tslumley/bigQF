tracefht_liu <-
function(A,k,trace.full=sum(A^2)){
	m<-NROW(A)
	mbig<-2^ceiling(log2(m))
	R<-sample(c(-1,1),m,replace=TRUE)
	Astar<-A*R
	idx<-sample(mbig,k)
	AOmega<-fht(Astar)[idx,]
	AAOmega<-tcrossprod(AOmega,A)
	A3Omega<-AAOmega%*%A
	A4Omega<-tcrossprod(A3Omega,A)
	
	tr<-sum(rowSums(AOmega*AOmega))/k
	trAA2<-sum(rowSums(AAOmega*AAOmega))/k
	trAA3<-sum(rowSums(A3Omega*A3Omega))/k
	trAA4<-sum(rowSums(A4Omega*A4Omega))/k

	c(trAA2*(trace.full/tr)^2, trAA3*(trace.full/tr)^3,trAA4*(trace.full/tr)^4)
}
