\name{tableELRHP}
\alias{tableELRHP}
\title{
Calculates \code{E(LR(HP))}
}
\description{
Calculates \code{E(LR(HP))} exactly, by formula and numerically for checking for
given number of alleles; answers is independent of allele frequencies
}
\usage{
tableELRHP(L = 4)
}
\arguments{
  \item{L}{
Integer, at least 1.
}
}

\value{
A table for a set of pairwise relationships.
}
\references{
Slooten and Egeland (to appear)
}
\author{
Thore Egeland Thore.Egeland@umb.no
}
\note{
Other pairwise relationships requires simple changes in code.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
tableELRHP(L = 2)
#Table in paper
#cbind(tableELRHP(L=2)[,-5],tableELRHP(L=4)[,4],tableELRHP(L=10)[,4])
}

