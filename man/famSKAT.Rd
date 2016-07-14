\name{famSKAT}
\alias{famSKAT}
\alias{update.famSKAT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implicit matrix for family-based SKAT test
}
\description{
Like \code{link{SKAT.matrixfree}} but for the family-based test of Chen and colleagues
}
\usage{
famSKAT(G, model, kinship,  weights = function(maf) dbeta(maf, 1, 25))
\S3method{update}{famSKAT}(object,G)
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
  \item{weights}{
A weight function used in SKAT: the default is the standard one.

}
\item{object}{An existing \code{famSKAT} object to be updated with a new set of genotypes (eg for a new gene or genomic window), keeping the same phenotype and kinship structure.  Avoids recomputing the Cholesky square root of the phenotype variance matrix, which will often be a computational bottleneck.
}
}


\value{
An object of class \code{c("famSKAT","matrixfree")}
}
\references{
Chen H, Meigs JM, Dupuis J (2013) Sequence Kernel Association Test for Quantitative Traits in Family Samples. Genet Epidemiol. 37(2): 196-204.
}

\note{
The matrix and test statistic both differ by a factor of \code{var(y)/2} from \code{SKAT.matrixfree} when used with unrelated individuals (because the Chen et al reference differs from the original SKAT paper by the same factor)
}


\seealso{
\code{\link{SKAT.matrixfree}}, \code{\link{pQF}}
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
