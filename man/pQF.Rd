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
  tr2.sample.size = 500, q = NULL,
  convolution.method = c("saddlepoint", "integration"),
  remainder.underflow=c("warn","missing","error"))
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
When \code{M} is not square, a randomised estimator for the trace of \code{crossprod(M)^2} is used. This is the sample size. When \code{M} is of \code{"matrixfree"} and does not include the trace of {crossprod(M)} as a component, this sample size will also be used to compute that. See Hutchinson reference. }
  \item{q}{
Power iteration parameter for the stochastic SVD (see the Halko et al reference)
}
  \item{convolution.method}{
For either the "ssvd" or "lanczos" methods, how the convolution of the leading-eigenvalue terms is performed: by Kuonen's saddlepoint approximation or Davies' algorithm inverting the characteristic function
}
\item{remainder.underflow}{
How should underflow of the remainder term (see Details) be handled? The default is a warning, with \code{"error"} an error is thrown, and with \code{"miss"} the return value will be \code{NaN}.
}
}
\details{
Increasing \code{neig} or \code{q} will improve the accuracy of approximation. In simulated human sequence data, \code{neig=100} is satisfactory for 5000 samples and markers. Increasing \code{q} should help when the singular values of \code{M} decrease slowly.  The default sample size for the randomized trace estimator seems to be large enough.

If the remainder term in the approximation has less than 1 degree of freedom it will be dropped, and \code{remainder.underflow} controls how this is handled. The approximation will then be anticonservative, but usually not seriously so. 

In earlier versions, \code{method="satterthwaite"} rounded the number of
degrees of freedom up to the next integer; it now does not round.

In the current version when \code{M} is of class \code{"matrix-free"}
and \code{method="ssvd"}, an improved estimator of the remainder term is
used. Instead of using Hutchinson's randomised trace estimator and
subtracting the known eigenvalues, we apply the randomised trace
estimator after projecting orthogonal to the known eigenvectors. After
more evaluation, the improvement is likely to be extended to the other
methods for rectangular matrices.

By default, Davies's algorithm is run with a tolerance of \code{1e-9}. This can be changed by setting, eg,  \code{options(bigQF.davies.threshold=1e-12)}


}
\value{
Vector of upper tail probabilities
}
\references{
Lumley et al. (2018) "Sequence kernel association tests for large sets of markers: tail probabilities for large quadratic forms" Genet Epidemiol. 2018 Sep;42(6):516-527. doi: 10.1002/gepi.22136

Tong Chen, Thomas Lumley (2019) Numerical evaluation of methods approximating the distribution of a large quadratic form in normal variables. Computational Statistics & Data Analysis.  139: 75-81,

Thomas Lumley (2017) "How to add chi-squareds" \url{https://notstatschat.rbind.io/2017/12/06/how-to-add-chi-squareds/}

Thomas Lumley (2016) "Large quadratic forms" \url{https://notstatschat.rbind.io/2016/09/27/large-quadratic-forms/}

Nathan Halko, Per-Gunnar Martinsson, Joel A. Tropp (2010) "Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions" \url{https://arxiv.org/abs/0909.4061}.

Hutchinson, M. F. (1990). A stochastic estimator of the trace of the influence matrix for Laplacian smoothing splines. Communications in Statistics - Simulation and Computation, 19(2):433-450.
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
skat$trace
pQF(c(33782471,7e7,1e8),skat, n=100)

\donttest{
# Don't run these; they take a few minutes
G<-sequence
wuweights<-function(maf) dbeta(maf,1,25)
tmp<-wuweights(colMeans(G)/2)*t(G)
tildeGt<-t(tmp-rowMeans(tmp))/sqrt(2)
sum(tildeGt^2)

pQF(c(33782471,7e7,1e8), tildeGt, n=100)

H<-crossprod(tildeGt)
pQF(c(33782471,7e7,1e8), H, n=100)
}


}