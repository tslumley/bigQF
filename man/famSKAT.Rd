\name{famSKAT}
\alias{famSKAT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implicit matrix for family-based SKAT test
}
\description{
Like \code{link{SKAT.matrixfree}} but for the family-based test of Chen and colleagues
}
\usage{
famSKAT(G, formula, kinship, data = NULL, weights = function(maf) dbeta(maf, 1, 25))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
A 0/1/2 matrix whose columns are markers and whose rows are samples. Should be mostly zero.

}
  \item{formula}{
A model formula for the null (adjustment) model
}
  \item{kinship}{
Sparse kinship matrix (ie, diagonal = 0.5, not 1)
}
  \item{data}{
Data frame for the variables in \code{formula}
}
  \item{weights}{
A weight function used in SKAT: the default is the standard one.

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
An object of class \code{matrixfree}
}
\references{
Chen H, Meigs JM, Dupuis J (2013) Sequence Kernel Association Test for Quantitative Traits in Family Samples. Genet Epidemiol. 37(2): 196–204.
}

\note{
The matrix and test statistic both differ by a factor of 1/2 from \code{SKAT.matrixfree} when used with unrelated individuals.
}


\seealso{
\code{\link{SKAT.matrixfree}}, \code{\link{pQF}}
}
\examples{
}
