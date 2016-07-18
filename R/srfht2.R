srfht2 <-
function(A,k,q=0){
	m<-NROW(A)
	mbig<-2^ceiling(log2(m))
	R<-sample(c(-1,1),m,replace=TRUE)
	Astar<-A*R
	idx<-sample(mbig,k)
	AOmega<-fht(Astar)[idx,]/sqrt(k)
	Q<-qr.Q(qr(t(AOmega)))
	for(i in seq_len(q)){
		tildeQ<-qr.Q(qr(A%*%Q))
		Q<-qr.Q(qr(crossprod(A,tildeQ)))
	}
	Q
	}
