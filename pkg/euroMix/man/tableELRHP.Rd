\name{tableELRHP}
\alias{tableELRHP}
\title{
Calculates \code{E(LR(HP))}
}
\description{
Calculates \code{E(LR(HP)), SD(LR(HP) and SD(LR(HD))} exactly
Answers is independent of allele frequencies except for \code{SD(LR(HD))}
}
\usage{
tableELRHP(L = 4,p=rep(1/L,L))
}
\arguments{
  \item{L}{
Integer, at least 1.
}
  \item{p}{
Allele frequencies, vector of length $L$
}
}

\value{
A table for a set of pairwise relationships.
}
\references{
Slooten and Egeland (to appear)
}
\author{
Thore Egeland Thore.Egeland@nmbu.no
}
\note{
Other pairwise relationships requires simple changes in code.
}



\examples{
tableELRHP(L = 2)
}

