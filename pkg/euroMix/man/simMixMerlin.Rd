\name{simMixMerlin}
\alias{simMixMerlin}
\title{
A DNA mixture is generated from individual genotypes using \code{paramlink} and \code{Merlin}
}
\description{
A \code{\link{linkdat}} object is created. MERLIN files can be generated or mixtures can be generated based on existing files.
This function requires MERLIN to be installed and correctly pointed to in the PATH environment variable. 
}
\usage{
simMixMerlin(x, aa, afreq, options=NULL, seed = 12345, generate = FALSE)
}

\arguments{
  \item{x}{
\code{\link{linkdat}} object
}
  \item{aa}{
allele list. \code{aa[[1]]} contains alleles for marker 1. MERLIN has an upper limit on the number of alleles
therefore has problems with one marker, SE33, in \code{\link{db}}
}
  \item{afreq}{
Frequency list
}
  \item{options}{
A character with additional options to pass on to Merlin
}
  \item{seed}{
  Random seed to pass on to Merlin. If not set, Merlin will return the same simulated data each time 
  }
  \item{generate}{
If TRUE, Merlin files are generated
}
}

\value{
\item{y }{\code{\link{linkdat}} object}
\item{comp2 }{list of mixtures}
}

\author{
Thore Egeland <Thore.Egeland@nmbu.no>
}


\examples{
\dontrun{
#Example 1
require(paramlink)
data(db)
x=cousinPed(1)
x=swapSex(addOffspring(x,father=7,mother=8,noff=2),ids=10)
db2=split(db,db$Marker)
Nmarkers=5
aa=vector("list",Nmarkers)
afreq=vector("list",Nmarkers)
for (i in 1:Nmarkers){
aa[[i]]=db2[[i]]$Allel
afreq[[i]]=db2[[i]]$Frequency
m=marker(x,9,c(1,1),10,c(1,1),alleles=1:length(aa[[i]]),afreq=afreq[[i]])
x=addMarker(x,m)
}
res=simMixMerlin(x,aa,afreq,generate=TRUE)
#The map file generated default above leads to tightly linked markers. The map file
#can be edited and simMixMerlin rerun with generate=FALSE. 
#Example 2 
#Next we consider an example #with two markers (for simplicity), 
#D12 and VWA, markers 3 and 23 in db and illustrate how the map 
#file is edited to account for linkage
x=cousinPed(1)
x=swapSex(addOffspring(x,father=7,mother=8,noff=2),ids=10)
Nmarkers=2
aa=vector("list",Nmarkers)
afreq=vector("list",Nmarkers)
i=0
for (j in c(3,23)){
i=i+1
aa[[i]]=db2[[i]]$Allel
afreq[[i]]=db2[[i]]$Frequency
m=marker(x,9,c(1,1),10,c(1,1),alleles=1:length(aa[[i]]),
afreq=afreq[[i]])
x=addMarker(x,m)
}
res=simMixMerlin(x,aa,afreq,generate=TRUE)
#Next edit map file, normally this is done  
#simpler than below
map=read.table("merlin.map",header=FALSE)
map[,1]=c(12,13)
map[,3]=c(0.5,0.5)
write.table(map,"merlin.map",col.names=FALSE,quote=FALSE,
row.names=FALSE)
res=simMixMerlin(x,aa,afreq,generate=FALSE)
}
}

