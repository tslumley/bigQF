\name{SKAT.matrixfree}
\alias{SKAT.matrixfree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make 'matrix-free' object for SKAT test
}
\description{
'Matrix-free' or 'implicit' linear algebra uses a matrix only through the linear operations of multiplying by the matrix or its transpose. It's suitable for sparse matrices, and also for structured matrices that are not sparse but still have fast algorithms for multiplication.  The Sequence Kernel Association Test is typically performed on sparse genotype data, but the matrix involved in computations has been centered and is no longer sparse. 
}
\usage{
SKAT.matrixfree(G,weights=function(maf) dbeta(maf,1,25))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
A 0/1/2 matrix whose rows are markers and whose columns are samples. Should be mostly zero.
}
  \item{weights}{
A weight function used in SKAT: the default is the standard one.
}
}

\value{
An object of class \code{matrixfree}

}
\references{
Lee, S., with contributions from Larisa Miropolsky, and Wu, M. (2015). SKAT: SNP-Set (Sequence) Kernel Association Test. R package version 1.1.2.

Lee, S., Wu, M. C., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011). Rare-variant association testing for sequencing data with the sequence kernel association test. American Journal of Human Genetics, 89:82-93.
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
pQF(c(1e10,1.2e10), skat,n=100)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
