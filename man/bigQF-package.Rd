\name{bigQF-package}
\alias{bigQF-package}
\alias{bigQF}
\docType{package}
\title{
\packageTitle{bigQF}
}
\description{
\packageDescription{bigQF}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{bigQF}
\packageIndices{bigQF}
This package computes tail probabilities for large quadratic forms, with
the motivation being the SKAT test used in DNA sequence association
studies.

The true distribution is a linear combination of 1-df chi-squared
distributions, where the coefficients are the non-zero eigenvalues of
the matrix \code{A} defining the quadratic form \eqn{z^TAz}.  The package uses an
approximation to the distribution consisting of the largest \code{neig} terms in the
linear combination plus the Satterthwaite approximation to the rest of
the linear combination. 

The main function is \code{\link{pQF}}, which has options for how to
compute the leading eigenvalues (Lanczos-type algorithm or stochastic
SVD) and how to compute the linear combination (inverting the
characteristic function or a saddlepoint approximation).  The Lanczos
algorithm is from the \code{svd} package; the stochastic SVD can be
called directly via \code{\link{ssvd}} or \code{\link{seigen}}

Given a square matrix, \code{pQF} uses it as \code{A}. If the input is a
non-square matrix \code{M}, then \code{A} is \code{crossprod(M)}. The
function can also be used matrix-free, given an object containing
functions to compute the product and transpose-product by \code{M}. This last option
is described in the \code{"matrix-free"} vignette. The matrix-free
algorithm also uses a randomised estimator to estimate
the trace of \code{crossprod(A)}. The function \code{sparse.matrixfree} constructs a object for
matrix-free use of \code{pQF} from a sparse Matrix object. The
algorithms are described in the Lumley et al (2018) reference.

Finally, there are functions specifically for the SKAT family of genomic
tests. These take a genotype matrix and an adjustment model as arguments
and produce an object that contains the test statistic in its
\code{Q} component and which can be used as an argument to \code{pQF} to
extract p-values: \code{\link{SKAT.matrixfree}} and \code{famSKAT}. The
vignette \code{"Checking pQF vs SKAT"} compares \code{SKAT.matrixfree}
to the \code{SKAT} package and illustrates how it can be used

}
\author{
\packageAuthor{bigQF}

Maintainer: \packageMaintainer{bigQF}
}
\references{
  
Tong Chen, Thomas Lumley (2019) Numerical evaluation of methods approximating the distribution of a large quadratic form in normal variables. Computational Statistics & Data Analysis.  139: 75-81,

Lumley et al. (2018) Sequence kernel association tests for large sets of markers: tail probabilities for large quadratic forms. Genet Epidemiol
. 2018 Sep;42(6):516-527. doi: 10.1002/gepi.22136

Nathan Halko, Per-Gunnar Martinsson, Joel A. Tropp (2010) Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions. \url{https://arxiv.org/abs/0909.4061}.

Lee, S., with contributions from Larisa Miropolsky, and Wu, M. (2015). SKAT: SNP-Set (Sequence) Kernel Association Test. R package version 1.1.2.

Lee, S., Wu, M. C., Cai, T., Li, Y., Boehnke, M., and Lin,
X. (2011). Rare-variant association testing for sequencing data with the
sequence kernel association test. American Journal of Human Genetics,
89:82-93.


}

