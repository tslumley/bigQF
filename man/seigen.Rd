\name{seigen}
\alias{seigen}
\alias{ssvd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stochastic singular value decomposition
}
\description{
Extract the leading eigenvalues of a large matrix and their eigenvectors, using random projections. 
}
\usage{
ssvd(M, n,U=FALSE, V=FALSE,q=3, p=10) 
seigen(M, n, only.values = TRUE, q = 3, symmetric = FALSE, spd = FALSE, p = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
A square matrix
}
  \item{n}{
Number of eigenvalues to extract
}
\item{U,V}{If \code{TRUE} return the left (respectively, right) singular vectors as well as the singular values}
  \item{only.values}{
If \code{TRUE}, only extract the eigenvalues, otherwise also extract corresponding eigenvectors
}
  \item{q}{
Number of power iterations to use in constructing the projection basis (zero or more)
}
  \item{symmetric}{
If \code{TRUE}, assume the matrix is symmetric
}
  \item{spd}{
If \code{TRUE}, assume the matrix is positive definite and use the Nystrom method to improve estimation. 
}
  \item{p}{
The oversampling parameter: number of extra dimensions above n for the random projection
}
}
\details{
The parameters p and q are as in the reference. Both functions use Algorithm 4.3 to construct a projection; \code{ssvd} then uses Algorithm 5.1. 
With \code{spd=TRUE}, \code{seigen} uses Algorithm  5.5, otherwise Algorithm 5.3
}
\value{
A list with components
\item{values}{eigenvalues}
\item{vectors}{matrix whose columns are the corresponding eigenvectors}
}
\references{
Nathan Halko, Per-Gunnar Martinsson, Joel A. Tropp (2010) "Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions" \url{https://arxiv.org/abs/0909.4061}.
}
\author{
Thomas Lumley
}
\note{
Unlike the Lanczos-type algorithms, this is accurate only for large matrices.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ssvd}}, \code{\link[base]{eigen}}
}
\examples{
data(sequence)
G<-sequence[1:1000,]
H<-tcrossprod(G)
seigen(H,n=10,spd=TRUE,q=5)


eigen(H, symmetric=TRUE,only.values=TRUE)$values[1:10]

  
}