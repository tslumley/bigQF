\name{SKAT.matrixfree}
\alias{SKAT.matrixfree}
\alias{SKAT.matrixfree.lm}
\alias{SKAT.matrixfree.glm}
\alias{SKAT.matrixfree.lmekin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make 'matrix-free' object for SKAT test
}
\description{
'Matrix-free' or 'implicit' linear algebra uses a matrix only through the linear operations of multiplying by the matrix or its transpose. It's suitable for sparse matrices, and also for structured matrices that are not sparse but still have fast algorithms for multiplication.  The Sequence Kernel Association Test is typically performed on sparse genotype data, but the matrix involved in computations has been centered and is no longer sparse.
}
\usage{
SKAT.matrixfree(G,weights=function(maf) dbeta(maf,1,25), model=NULL,...)
\method{SKAT.matrixfree}{lm}(G,weights=function(maf) dbeta(maf,1,25), model=NULL,...)
\method{SKAT.matrixfree}{glm}(G,weights=function(maf) dbeta(maf,1,25), model=NULL,...)
\method{SKAT.matrixfree}{lmekin}(G,weights=function(maf) dbeta(maf,1,25), model=NULL, kinship,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
A 0/1/2 matrix whose columns are markers and whose rows are samples. Should be mostly zero.
}
  \item{weights}{
A weight function used in SKAT: the default is the standard one.
}
\item{model}{
A  linear model, generalised linear model, or \code{lmekin} object used for adjustment, or \code{NULL}
}
\item{kinship}{A (sparse) kinship matrix. The \code{model} will have used \code{2*kinship} in its
\code{varlist} argument.
}
\item{\dots}{to keep CMD check happy}
}
\value{
An object of class \code{matrixfree}
}
\details{
If the adjustment model is \code{NULL} the object contains the trace of the underlying quadratic form, if the adjustment model is not \code{NULL} the trace will be estimated using Hutchinson's randomised estimator inside \code{\link{pQF}}.  The \code{lmekin} method calls \code{\link{famSKAT}}.

}

\references{
Lee, S., with contributions from Larisa Miropolsky, and Wu, M. (2015). SKAT: SNP-Set (Sequence) Kernel Association Test. R package version 1.1.2.

Lee, S., Wu, M. C., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011). Rare-variant association testing for sequencing data with the sequence kernel association test. American Journal of Human Genetics, 89:82-93.

Wu, M. C., Kraft, P., Epstein, M. P., Taylor, D. M., Channock, S. J., Hunter, D. J., and Lin, X. (2010). Powerful SNP set analysis for case-control genome-wide association studies. American Journal of Human Genetics, 86:929-942.
}
\author{
Thomas Lumley
}



\seealso{
\code{\link{pQF}}
\code{\link{sparse.matrixfree}}
}
\examples{
data(sequence)
skat<-SKAT.matrixfree(sequence)
skat$trace
pQF(c(33782471,7e7,1e8), skat,n=100, tr2.sample.size=500)

data(SKAT.example)

skat1mf <- SKAT.matrixfree(SKAT.example$Z)
(Q<-skat1mf$Q(SKAT.example$y.c))
all.equal(as.numeric(Q), 234803.786)

## correct value is 0.01874576
pQF(Q, skat1mf, neig=4, convolution.method="integration")


skat3mf <- SKAT.matrixfree(SKAT.example$Z, model=lm(y.c~X, data=SKAT.example))
(Q<-skat3mf$Q())
all.equal(Q, 298041.542)

## correct value is 0.002877041
pQF(Q, skat3mf, neig=4, convolution.method="integration")


}