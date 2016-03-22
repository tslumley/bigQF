\name{sparse.matrixfree}
\alias{sparse.matrixfree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make 'matrix-free' object from (sparse) Matrix
}
\description{
Packages a matrix (which will typically be a sparse Matrix) for 'matrix-free' or 'implicit' stochastic SVD. 
}
\usage{
sparse.matrixfree(M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
A matrix or Matrix
}
}


\value{
Object of class \code{'matrixfree'}, with components
\item{mult}{Function to multiply by \code{M}}
\item{tmult}{Function to multiply by \code{t(M)}}
\item{trace}{trace of \code{t(M)\%*\%M}, needed for \code{\link{pQF}}}
\item{ncol}{dimensions of \code{M}}
\item{nrow}{dimensions of \code{M}}
}

\author{
Thomas Lumley
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pQF}},\code{link{seigen}},\code{\link{ssvd}}, \code{\link{SKAT.matrixfree}}
}
\examples{
data(sequence)
Msp<-sparse.matrixfree(sequence)
ssvd(Msp,n=10)

\dontrun{
## this is slow
svd(sequence,nu=0,nv=0)$d[1:10]
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
