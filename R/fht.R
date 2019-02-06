small_fht <-
function(A){
	m<-NROW(A)
	mbig<-2^ceiling(log2(m))
	Abig<-matrix(0,nrow=mbig,ncol=NCOL(A))
	Abig[1:m,]<-A
	.C("mfwht", Abig, as.integer(mbig), as.integer(NCOL(Abig)))[[1]]
}


fht<-function (A,big=TRUE) 
{
    m <- NROW(A)
    mbig <- 2^ceiling(log2(m))
    Abig <- matrix(0, nrow = mbig, ncol = NCOL(A))
    Abig[1:m, ] <- A
    if (big)
	    .Call("big_mfwht",Abig)
    else 
	    .C("mfwht", Abig, as.integer(mbig), as.integer(NCOL(Abig)))[[1]]
}
