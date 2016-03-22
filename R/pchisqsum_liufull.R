pchisqsum_liufull <-
function (q, H) 
{
	c1<-sum(diag(H))
	H2<-crossprod(H)
	c2 <- sum(diag(H2))
	H3<-crossprod(H,H2)
    c3 <- sum(diag(H3))
    c4 <- sum(H3^2)
    s1 <- c3/(c2^(3/2))
    s2 <- c4/c2^2
    muQ <- c1
    sigmaQ <- sqrt(2 * c2)
    tstar <- (q - muQ)/sigmaQ
    if (s1^2 > s2) {
        a <- 1/(s1 - sqrt(s1^2 - s2))
        delta <- s1 * a^3 - a^2
        l <- a^2 - 2 * delta
    }
    else {
        a <- 1/s1
        delta <- 0
        l <- c2^3/c3^2
    }
    muX <- l + delta
    sigmaX <- sqrt(2) * a
    Qq <- pchisq(tstar * sigmaX + muX, df = l, ncp = delta, lower.tail = FALSE)
    return(Qq)
}
