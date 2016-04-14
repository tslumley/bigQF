
SKAT.matrixfree.lmekin<-function(G,  weights = function(maf) dbeta(maf, 1, 25), model, kinship){
  famSKAT(G, model, kinship, weights)
}

famSKAT<-function(G,  model, kinship,  weights = function(maf) dbeta(maf, 1, 25)){
	center <- colMeans(G)
    ww <- weights(center/2)
    spG <- Matrix(G, sparse = TRUE) %*% Diagonal(x = ww)

	
    model$theta <- c(model$vcoef$id, model$sigma^2)
    SIGMA <- model$theta[1] * 2 * kins + model$theta[2] * Diagonal(nrow(kins))
    CholSigma<-t(chol(SIGMA)) 
    if (is.null(model$x)) stop("null model must be fitted with x=TRUE")
    Z<-model$x
    qr<-qr(as.matrix(solve(CholSigma,Z)))
	
	#debugging
	# sinvZ<-solve(CholSigma,Z)
	# Pi0<- Diagonal(600)-sinvZ%*%solve(crossprod(sinvZ))%*%t(sinvZ)
	# tildeG<-t(spG)%*%solve(SIGMA)%*%CholSigma%*%Pi0
	# H<-tcrossprod(tildeG)


    rval <- list(mult = function(X) {
    		base::qr.resid(qr,as.matrix(solve(CholSigma,(spG %*% X))))
    }, tmult = function(X) {
        crossprod(spG, solve(t(CholSigma), base::qr.resid(qr,X)))
    }, trace = NULL, ncol = ncol(G), nrow = nrow(G))
    class(rval) <- "matrixfree"
    rval
}