\name{pQF}
\alias{pQF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tail probabilities for quadratic forms
}
\description{
Computes the upper tail probability for quadratic forms in standard Normal variables, using a leading-eigenvalue approximation that is feasible even for large matrices. Can take advantage of sparse matrices. 
}
\usage{
pQF(x, M, method = c("ssvd", "lanczos", "satterthwaite"), neig = 100,
  tr2.sample.size = 300, q = NULL,
  convolution.method = c("saddlepoint", "integration"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of quantiles to compute tail probabilities
}
  \item{M}{
If \code{M} is square, it is the matrix in the quadratic form and is assumed to be symmetric. If it is not square, \code{crossprod(M)} is the matrix in the quadratic form, although this matrix is never formed explicitly. It can also be an object of class \code{"matrixfree"} to allow matrix-free multiplication for, eg, sparse matrices; see \code{\link{SKAT.matrixfree}} for an example.
}
  \item{method}{ Use stochastic SVD ("ssvd") or a thick-restarted Lanczos algorithm ("lanczos") to extract the leading eigenvalues, or use the naive Satterthwaite approximation ("satterthwaite")}

  \item{neig}{
Number of leading eigenvalues to use
}
  \item{tr2.sample.size}{
When \code{M} is not square, a randomised estimator for the trace of $crossprod(M)^2$ is used. This is the sample size.}
  \item{q}{
Power iteration parameter for the stochastic SVD (see the Tropp et al reference)
}
  \item{convolution.method}{
For either the "ssvd" or "lanczos" methods, how the convolution of the leading-eigenvalue terms is performed: by Kuonen's saddlepoint approximation or Davies' algorithm inverting the characteristic function
}
}
\details{
Increasing \code{neig} or \code{q} will improve the accuracy of approximation. In simulated human sequence data, \code{neig=100} is satisfactory for 5000 samples and markers. Increasing \code{q} will help when the singular values of \code{M} decrease slowly.  The default sample size for the randomized trace estimator seems to be large enough.
}
\value{
Vector of upper tail probabilities
}
\references{
Lumley et al. "Sequence kernel association tests for large sets of markers: tail probabilities for large quadratic forms"

Nathan Halko, Per-Gunnar Martinsson, Joel A. Tropp (2010) "Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions" \url{http://arxiv.org/abs/0909.4061}.
}
\author{
Thomas Lumley
}

\seealso{
\code{\link{ssvd}},\code{\link{SKAT.matrixfree}}
}
\examples{
data(sequence)
dim(sequence)

skat<-SKAT.matrixfree(sequence)
pQF(c(1e10,1.2e10),skat, n=100)

\dontrun{
# these take a few minutes
scG<-scale(sequence)
wuweights<-function(maf) dbeta(maf,1,25)
tildeGt<-t(wuweights(attr(scG,"scaled:center")/2)*t(scG))
pQF(c(1e10,1.2e10),tildeGt, n=100)

H<-crossprod(tildeGt)
pQF(c(1e10,1.2e10), H, n=100)
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
