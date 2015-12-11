\name{LRstat}
\alias{LRstat}

\title{
Distribution of LR(HP) and LR(HD)
}
\description{
Distribution of LR(HP) and LR(HD) are calculated as well as some summary statistics. 
}
\usage{
LRstat(ped_claim, ped_true, ids, alleles, afreq = NULL, 
known_genotypes = list(), loop_breakers = NULL, Xchrom = F, plot = T)
}
\arguments{
  \item{ped_claim}{
a linkdat object, or a list of several linkdat and/or singleton objects, describing
the claimed relationship. If a list, the sets of ID labels must be disjoint, that is, 
all ID labels must be unique.
}
  \item{ped_true}{
a linkdat object, or a list of several linkdat and/or singleton objects,
describing the true relationship. 
ID labels must be consistent with ped_claim.
individuals available for genotyping.
}
  \item{ids}{
individuals available for genotyping.
}
  \item{alleles}{
a numeric or character vector containing marker alleles names
}
  \item{afreq}{
a numerical vector with allele frequencies. 
An error is given if they don't sum to 1 (rounded to 3 decimals).
}
  \item{known_genotypes}{
list of triplets (a, b, c), indicating that individual a has genotype b/c.
}
  \item{loop_breakers}{
a numeric containing IDs of individuals to be used as loop breakers.
Relevant only if any of the pedigrees has loops. See breakLoops.
}
  \item{Xchrom}{
a logical: Is the marker on the X chromosome?
}
  \item{plot}{
either a logical or the character "plot_only", controlling if a plot should be produced. 
If "plot_only", a plot is drawn, but no further computations are done 
(useful for reproducing the plot in computer-intensive applications)
}
}
\details{
Connected to joint work with Klaas Slooten
}
\value{
\item{main }{Expected values, variances of LR(HP) and LR(HD). RMNE}
\item{extra }{P(data|ped_claim) and P(data|ped_true)}
\item{LRdist }{Distribution of LR}
}

\author{
Thore Egeland Thore.Egeland@nmbu.no
}

\seealso{
See Also \code{\link{exclusionPower}}
}
\examples{
HP = nuclearPed(noffs=1, sex=2) # Specifies individual 1 as the father of 3
HD= list(singleton(id=1,sex=1), singleton(id=3, sex=2)) # Specifies 1 and 3 as unrelated
p=c(0.2,0.3,0.5);L=length(p)
available = c(1, 3)
res=LRstat(HP, HD, available, alleles = 1:L, afreq=p)
E.LR.HP=res$main[1]
stopifnot(abs(E.LR.HP-(L+3)/4)<1e-06)
res$LRdist #Distribution of LR
}
