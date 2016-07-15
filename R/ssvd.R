ssvd<-function(M,n,U=FALSE, V=FALSE,q=3, p=10){
   if(inherits(M,"matrixfree")){
        Q<-matfreeQ(M$mult,M$tmult, n+p, M$ncol,q=q)
	B<-M$tmult(Q)
   } else {
	Q<-srfht2(M,n+p,q=q)
	B<-M%*%Q
   }
   s<-svd(B,nu=if (U) n else 0, nv=if(V) n else 0)
   s$d<-s$d[1:n]
   if (V) s$v<-Q%*%s$v
   s
}


seigen<-function(M, n, only.values=TRUE, q=3, symmetric=FALSE,spd=FALSE,p=10){
	Q<-srfht(M,n+p,q=q)
	B1<-M%*%Q
	B2<-crossprod(Q,B1)
	if (spd){
	  C<-chol(B2)
	  F<-backsolve(C,t(B1))
	  s<-svd(F,nu=0, nv= if(only.values) 0 else n)
	  ee<-list(vectors=s$v,values=s$d[1:n]^2)
	} else {
   	  ee <- eigen(B2, symmetric = symmetric, only.values = only.values)
    	  ee$values <- ee$values[1:n]
    	  if (!only.values) 
             ee$vectors <- Q %*% ee$vectors[, 1:n]
	}
	ee
}
