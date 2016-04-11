
famSKAT<-function(G,  formula, kinship, data=NULL, weights = function(maf) dbeta(maf, 1, 25)){
	center <- colMeans(G)
    ww <- weights(center/2)
    spG <- Matrix(G, sparse = TRUE) %*% Diagonal(x = ww)

	nullmodel <- coxme::lmekin(formula = update(formula, 
            "~.+ (1|id)"), data = data, varlist = 2 * kins, method = "REML",x=TRUE)
    nullmodel$theta <- c(nullmodel$vcoef$id, nullmodel$sigma^2)
    SIGMA <- nullmodel$theta[1] * 2 * kins + nullmodel$theta[2] * Diagonal(nrow(kins))
	CholSigma<-t(chol(SIGMA))
	Z<-nullmodel$x
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
    }, trace = NULL, ncol = ncol(G), nrow = nrow(G), debugging=list(Pi0=Pi0,tildeG=tildeG,H=H))
    class(rval) <- "matrixfree"
    rval
}
