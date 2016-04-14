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
famSKAT(G, model, kinship, data = NULL, weights = function(maf) dbeta(maf, 1, 25))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
A 0/1/2 matrix whose columns are markers and whose rows are samples. Should be mostly zero.

}
  \item{model}{
Object of class \code{lmekin} giving the adjustment model, which must have been fitted with \code{x=TRUE}
}
  \item{kinship}{
The sparse kinship matrix: the \code{model} will have used \code{2*kinship} in its \code{varlist} argument.
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
The matrix and test statistic both differ by a factor of \code{var(y)/2} from \code{SKAT.matrixfree} when used with unrelated individuals.
}


\seealso{
\code{\link{SKAT.matrixfree}}, \code{\link{pQF}}
}
\examples{
}
