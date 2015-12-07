\name{q012}
\alias{q012}
\title{
Probabilities for pairwise relationships
}
\description{
Calculates the probability distribution for a pair of
individuals conditionally on 0,1, and 2 IBD alleles.
}
\usage{
q012(p = c(0.5, 0.5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
A numerical vector with allele frequencies
}
}
\details{
The function calls \code{\link{oneMarkerDistribution}} for IBD=0,1 and 2.
}
\value{
\item{q0}{Joint distribution given IBD=0}
\item{q1}{Joint distribution given IBD=1}
\item{q2}{Joint distribution given IBD=2}
}
\references{
None
}
\author{
Thore Egeland <Thore.Egeland@nmbu.no>
}


\examples{
require(paramlink)
q012()
}


