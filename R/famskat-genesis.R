##
## Not done yet
##

famSKAT.mattmodel<-function (G, model, weights = function(maf) dbeta(maf, 
    1, 25)) 
{
    center <- colMeans(G)
    ww <- weights(center/2)
    spG <- Matrix(G, sparse = TRUE) %*% Diagonal(x = ww)
    CholSigmaInv <- Matrix(model$cholSigmaInv)
    SIGMAinv <- tcrossprod(CholSigmaInv)
    Z <- model$model.matrix
    qr <- qr(as.matrix(CholSigmaInv%*%Z))
    mu <- Z %*% model$fixef[,1]
    res <- model$workingY - mu
    rval <- list(mult = function(X) {
        base::qr.resid(qr, as.matrix(CholSigmaInv%*% (spG %*% 
            X)))
    }, tmult = function(X) {
        crossprod(spG, crossprod(CholSigmaInv, base::qr.resid(qr, 
            X)))
    }, trace = NULL, ncol = ncol(G), nrow = nrow(G), Q = function() {
        stdres <- SIGMAinv%*%res
        s = crossprod(spG, stdres)
        sum(s^2)
    })
    class(rval) <- "matrixfree"
    rval
}
