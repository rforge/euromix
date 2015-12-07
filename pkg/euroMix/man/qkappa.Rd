\name{qkappa}
\alias{qkappa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates joint distribution for a pair of individals given IBD probabilities
}
\description{
Based on conditional distribution given IBD from \code{\link{q012}},
the joint probability distribution for two individuals are given for
specified IBD probabilities
}
\usage{
qkappa(kappa = c(0, 1, 0), q = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kappa}{
Three reals summing to 1 giving IBD (0,1,2) probabilities
}
  \item{q}{
The joint probability distribution for two individuals 
}
}

\value{
%%  ~Describe the value returned
A matrix giving the joint distribution
}
\references{
To appear
}
\author{
Thore Egeland <Thore.Egeland@nmbu.no>
}

\examples{
require(paramlink)
#Sibs. One SNP marker with 
qkappa(kappa=c(0.25,0.5,0.25),q012(p=c(0.2,0.8)))

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
