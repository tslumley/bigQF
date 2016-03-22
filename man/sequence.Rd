\name{sequence}
\alias{sequence}
\docType{data}
\title{
Simulated human DNA variant sequence
}
\description{
A matrix with number of copies of the minor allele for 4028 variants on 5000 people, simulated using the Markov Coalescent Simulator of Chen and coworkers.
}
\usage{data("sequence")}
\format{
  The format is:
 num [1:5000, 1:4028] 0 0 0 0 0 0 0 0 0 0 ...
}

\source{
\url{https://github.com/gchen98/macs}

}
\references{
Gary K Chen, Paul Marjoram, and Jeffrey D. Wall (2009) Fast and flexible simulation of DNA sequence data Genome Research 19:136-142. \url{http://genome.cshlp.org/content/19/1/136}
}
\examples{
data(sequence)
summary(colMeans(sequence))
}
\keyword{datasets}
