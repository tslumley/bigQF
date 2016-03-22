matfreeQ <-
function(Mmult,tMmult,m,cols,q=0){
	Omega<-matrix(rnorm(m*cols), ncol=m)
	AOmega<-Mmult(Omega)
	Q<-qr.Q(qr(AOmega))
	for(i in seq_len(q)){
		tildeQ<-qr.Q(qr(tMmult(Q)))
		Q<-qr.Q(qr(Mmult(tildeQ)))
	}
	Q
	}
