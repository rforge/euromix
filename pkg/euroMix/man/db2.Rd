\name{db2}
\alias{db2}
\docType{data}
\title{
Allele database.
}
\description{
Norwegian database for 10 SGM Plus markers.
}
\usage{data(db2)}
\format{
  A data frame with 119 observations on the following 3 variables.
  \describe{
    \item{\code{Marker}}{a factor with levels corresponding to name of markers}
    \item{\code{Allele}}{a numeric vector denoting allele}
    \item{\code{Frequency}}{a numeric vector in (0,1)}
  }
}
\details{
The format is convenient for R.
}
\source{
Andreassen et al. (2007).
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(db2)
#Checks that frequencies add to 1
lapply(split(db2$Frequency,db2$Marker),sum)
#Finds number of alleles for all markers
unlist(lapply(split(db2$Frequency,db2$Marker),length))
#A closer look at the marker TH01
TH01=db2[db2$Marker=="TH01",]
barplot(TH01$Frequency)
}
\keyword{datasets}
