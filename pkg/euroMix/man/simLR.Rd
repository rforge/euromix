\name{simLR}
\alias{simLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood for mixtures that may have related contributors and drop-in and drop-out of alleles
}
\description{
Likelihood for mixtures that may have related contributors and drop-in and drop-out of alleles. For a general description of the problem, see see \code{\link{paraMix}}. As opposed to \code{\link{paraMix}}, drop-in and drop-out of alleles are allowed. The likelihood is based on simulations from an urn model. Possible mixtures are simulated by applying drop-in and drop-out to genotypes for the assumed contributors. Genotypes for unknown contributors are simulated conditioned on the pedigree.
}
\usage{
simLR(R, x, alleles, afreq, pDO, pDI, N, known_genotypes = NULL, 
ped = NULL, id.U = NULL, id.V = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{
  Integers, mixture
}
  \item{x}{
  Number of unknown contributors
}
  \item{alleles}{
  Integers indicating alleles for marker 
}
  \item{afreq}{
  A numerical vector with allele frequencies
}
  \item{pDO}{
  Probability of drop-out applied per allele
}
  \item{pDI}{
  Probability of drop-in per locus
}
  \item{N}{
  Number of simulations
}
  \item{known_genotypes}{
  List of known genotypes. If a pedigree is specified, each element must a triplet of integers corresponding to (id,allele1,allele2).
  If no pedigree is specified, the id can be omitted.
}
  \item{ped}{
  \code{\link{linkdat}} object, or a list of such (if disconnected), describing the claimed relationship.
}
  \item{id.U}{
  Integers indicating untyped contributors (e.g.,suspect(s)). Only relevant if a pedigree is specified.
}
  \item{id.V}{
  Integers indicating typed non-contributors. Only relevant if a pedigree is specified.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
p.R: the likelihood of the mixture R
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Guro Dorum and Thore Egeland <guro.dorum@nmbu.no>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{paraMix}}.
}
\examples{
require(paramlink)
alleles <-1:4
p <- c(0.044, 0.166, 0.11, 0.68)
names(p) <- alleles
R <- 1:3
known <- list(c(6,1,2),c(4,3,3))
x <- halfCousinPed(0)
y <- singleton(6,sex=2)
pDO <- 0.1
pDI <- 0.05
N <- 20000
lp <- simLR(R=R, x=0, alleles=alleles, afreq=p, 
pDO, pDI, N, known_genotypes=known, ped=list(x,y))
ld <- simLR(R=R, x=1, alleles=alleles, afreq=p, pDO, pDI, 
N, known_genotypes=known, ped=list(x,y), id.U=5,id.V=4)
lp/ld
}
