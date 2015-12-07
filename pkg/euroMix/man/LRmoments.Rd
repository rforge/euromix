\name{LRmoments}
\alias{LRmoments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates expectation, standard deviation and skewness of LR under HP and HD
}
\description{
Exact numerical calculation
}
\usage{
LRmoments(p = c(0.5, 0.5), kappaP = c(0, 1, 0),kappaD = c(1, 0, 0),log10=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
Allele frequencies
}
  \item{kappaP}{
Probabilities of 0,1 and 2 alleles IBD corresponding to pedigree for HP
}
  \item{kappaD}{
Probabilities of 0,1 and 2 alleles IBD corresponding to pedigree for HD
}
  \item{log10}{
If \code{TRUE, LR} is log10 transformed.
}
}

\value{
\item{moments }{expectation, standard deviation and skewness of LR under HP and HD}
\item{LRtable}{Distribution of LR under HP and HD}

}
\references{
Slooten and Egeland (2013, submitted)
}
\author{
Thore Egeland Thore.Egeland@nmbu.com
}

\examples{
LRmoments(kappaP=c(0,1,0)) #Motivating example
LRmoments(kappaP=c(0,0.25,0.75)) #skew(LR(HP))<0
#Appendix of Slooten and Egeland (2013, submitted)
\dontrun{data(db)
p=db[db$Marker=="VWA",]$Freq
LRmoments(p=p,kappaP=c(0,1,0))}


}
