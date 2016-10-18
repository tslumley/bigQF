##
## Needs optimising
##


famSKAT.GENESIS.nullMixedModel <-function (G, model,threshold=1e-10, weights = function(maf) dbeta(maf, 
    1, 25),...) 
{
    center <- colMeans(G)
    ww <- weights(center/2)
    spG <- Matrix(G, sparse = TRUE) %*% Diagonal(x = ww)
    CholSigmaInv <- Matrix(model$cholSigmaInv)
    SIGMA<-chol2inv(t(model$cholSigmaInv))
    SIGMA[abs(SIGMA)<threshold]<-0
    SIGMA<-Matrix(SIGMA)
    CholSigma<-t(chol(SIGMA))
    SIGMAinv <- tcrossprod(CholSigmaInv)
    Z <- model$model.matrix
    qr<-qr(as.matrix(solve(CholSigma,Z)))
    mu <- Z %*% model$fixef[,1]
    res <- model$workingY - mu


    rval <- list(mult = function(X) {
        base::qr.resid(qr,as.matrix(solve(CholSigma,(spG %*% X))))
    }, tmult = function(X) {
        crossprod(spG, solve(t(CholSigma), base::qr.resid(qr,X)))
    },
    trace = NULL,
    ncol = ncol(G),
    nrow = nrow(G),
     Q = function() {
        stdres <- SIGMAinv%*%res
        s = crossprod(spG, stdres)
        sum(s^2)
    })
    

    class(rval) <- c("famSKAT_genesis","famSKAT","matrixfree")
    rval
}


update.famSKAT_genesis<-function(object, G,...){
    center <- colMeans(G)

    
    ## Get the computationally expensive bits from the existing object
    e<-environment(object$mult)
    weights<-get("weights",e)
    
    ww <- weights(center/2)
    spG <- Matrix(G, sparse = TRUE) %*% Diagonal(x = ww)

    CholSigma<-get("CholSigma",e)
    SIGMA<-get("SIGMA",e)
    SIGMAinv<-get("SIGMAinv",e)
    qr<-get("qr",e)
    res<-get("res",e)

    ## return new object
    rval <- list(mult = function(X) {
        base::qr.resid(qr,as.matrix(solve(CholSigma,(spG %*% X))))
    }, tmult = function(X) {
        crossprod(spG, solve(t(CholSigma), base::qr.resid(qr,X)))
    },
    trace = NULL,
    ncol = ncol(G),
    nrow = nrow(G),
    Q=function(){
        stdres<-solve(SIGMA,res)
        s=crossprod(spG, stdres)
        sum(s^2)
    }
    )
    class(rval) <- c("famSKAT_lmekin","famSKAT","matrixfree")
    rval
    
}
