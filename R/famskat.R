
SKAT.matrixfree.lmekin<-function(G,  weights = function(maf) dbeta(maf, 1, 25), model=NULL, kinship,...){
    famSKAT(G, model, kinship, weights)
}


famSKAT<-function(G, model, ...) UseMethod("famSKAT",model)

famSKAT.lmekin<-function(G,  model, kinship,  weights = function(maf) dbeta(maf, 1, 25),...){
    center <- colMeans(G)
    ww <- weights(center/2)
    spG <- Matrix(G, sparse = TRUE) %*% Diagonal(x = ww)


    model$theta <- c(model$vcoef$id, model$sigma^2)
    SIGMA <- model$theta[1] * 2 * kinship + model$theta[2] * Diagonal(nrow(kinship))
    CholSigma<-t(chol(SIGMA))
    if (is.null(model$x) || is.null(model$y)) stop("null model must be fitted with x=TRUE, y=TRUE")
    Z<-model$x
    qr<-qr(as.matrix(solve(CholSigma,Z)))
    mu <- model$x%*%model$coefficients$fixed
    res<-model$y-mu

    ##debugging
    ## sinvZ<-solve(CholSigma,Z)
    ## Pi0<- Diagonal(600)-sinvZ%*%solve(crossprod(sinvZ))%*%t(sinvZ)
    ## tildeG<-t(spG)%*%solve(SIGMA)%*%CholSigma%*%Pi0
    ## H<-tcrossprod(tildeG)


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

update.famSKAT_lmekin<-function(object, G,...){
    center <- colMeans(G)

    
    ## Get the computationally expensive bits from the existing object
    e<-environment(object$mult)
    weights<-get("weights",e)
    
    ww <- weights(center/2)
    spG <- Matrix(G, sparse = TRUE) %*% Diagonal(x = ww)

    CholSigma<-get("CholSigma",e)
    SIGMA<-get("SIGMA",e)
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
