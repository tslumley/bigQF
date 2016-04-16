\name{seqMetaExample}
\alias{seqMetaExample}
\alias{Z1}
\alias{Z2}
\alias{pheno1}
\alias{pheno2}
\alias{kins}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Example data, from seqMeta package
}
\description{
Contains simulated data for two cohorts taken from the seqMeta package. The individual genes are too small for this to be a good use of the leading-eigenvalue approximation, but the data at least allow basic numerical comparisons.
}
\usage{data(seqMetaExample)}
\format{
	This contains simulated data for two cohorts
  \describe{
    \item{\code{Z1},\code{Z2}}{
    	Genotype matrices for cohorts 1 and 2 respectively
    }
    \item{\code{pheno1},\code{pheno2}}{
    	phenotype matrices for cohorts 1 and 2 respectively
    }
    \item{kins}{
    	The kinship matrix for cohort 2
    }
  }
}
\source{
\url{https://github.com/DavisBrian/seqMeta}
}
\examples{
data(seqMetaExample)

m<-coxme::lmekin(y~sex+bmi+(1|id),data=pheno2,varlist=2*kins,
	x=TRUE,y=TRUE,method="REML")

#first gene
g1snps<-c("1000001", "1000002", "1000003", "1000004", "1000005", "1000006", 
"1000007", "1000008", "1000009", "1000010", "1000012", "1000013", 
"1000014", "1000015")
Z2gene1<-Z2[,g1snps]

f<-famSKAT(Z2gene1, m, kins)
Q<-f$Q()
all.equal(Q, 56681.209)
## correct p is 0.742756401
pQF(Q,f,neig=4)


}
\keyword{datasets}
