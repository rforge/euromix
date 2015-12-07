\name{simMixParamlink}
\alias{simMixParamlink}

\title{
Generates DNA mixtures
}
\description{
A DNA mixture is generated from individual genotypes using the R package paramlink
}
\usage{
simMixParamlink(y, alleles)
}

\arguments{
  \item{y}{
\code{linkdat} object from paramlink 
}
\item{alleles}{
Alleles in original form.
}
}
\details{
The alleles are internally represented as consecutive integers \code{1,2,...}, mixtures are generated
and transferred back to original allele values
}
\value{
A list of length equal to the number of markers (or simulations) each giving the mixture
}

\author{
Thore Egeland <Thore.Egeland@nmbu.no>
}

\examples{
#Example 1
require(paramlink)
x=cousinPed(1)
x=swapSex(addOffspring(x,father=7,mother=8,noff=2),ids=9)
plot(x)
data(db)
locus="FGA"
afreq1=db[db$Marker==locus,3]
alleles=db[db$Marker==locus,2]
m1=marker(x,alleles=alleles,afreq=afreq1)
y=markerSim(x,N=3,available=c(9,10),partialmarker=m1,verbose=FALSE,loop=7,seed=2)
res=simMixParamlink(y,alleles)
plot(y,marker=1:3)

#Example 2  With conditioning
x=halfCousinPed(2)
data(db)
locus="FGA"
afreq1=db[db$Marker==locus,3]
alleles=db[db$Marker==locus,2]
g.13=c(18,19.2)
m1=marker(x,13,g.13,alleles=alleles,afreq=afreq1)
y=markerSim(x,N=2,available=c(8,9,12),partialmarker=m1,verbose=FALSE,seed=2)
res=simMixParamlink(y,alleles)
plot(y,marker=1:2,cex=0.7,starred=13)

}



