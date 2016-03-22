pchisqsum_liu <-
function (q, M, tr2.sample.size=300) 
{
	c1<-sum(M^2)
	traces<-tracefht_liu(M,k=tr2.sample.size,trace.full=c1)
	c2 <- traces[1]
    c3 <- traces[2]
    c4 <- traces[3]
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
