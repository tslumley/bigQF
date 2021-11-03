\name{SKAT.example}
\alias{SKAT.example}
\docType{data}
\title{
Data example from SKAT package
}
\description{
These data (probably synthetic) come from the SKAT package. The data set is too small for the leading-eigenvalue approximation to really make sense, but it provides some numerical comparison.  The SKAT Q statistic should match exactly, the p-values should be fairly close.
}
\usage{data("SKAT.example")}
\format{
  SKAT.example contains the following objects:
  \describe{
    \item{Z}{a numeric genotype matrix of 2000 individuals and 67 SNPs. 
	Each row represents a different individual, and each column represents a different SNP marker. }
    \item{X}{a numeric matrix of 2 covariates. }
    \item{y.c}{a numeric vector of continuous phenotypes.}
    \item{y.b}{a numeric vector of binary phenotypes.}
  }
}

\source{
\url{https://www.hsph.harvard.edu/skat/}
}

\examples{
data(SKAT.example)

skat1mf <- SKAT.matrixfree(SKAT.example$Z)
Q<-skat1mf$Q(SKAT.example$y.c)
all.equal(as.numeric(Q), 234803.786)
## correct value is 0.01874576
pQF(Q, skat1mf, neig=4, convolution.method="integration")

skat2mf <- SKAT.matrixfree(SKAT.example$Z, model=lm(y.c~1, data=SKAT.example))
Q<-skat2mf$Q()
all.equal(Q, 234803.786)
## correct value is 0.01874576
pQF(Q, skat2mf, neig=4, convolution.method="integration")

skat3mf <- SKAT.matrixfree(SKAT.example$Z, model=lm(y.c~X, data=SKAT.example))
Q<-skat3mf$Q()
all.equal(Q, 298041.542)
## correct value is 0.002877041
pQF(Q, skat3mf, neig=4, convolution.method="integration")


}
\keyword{datasets}
