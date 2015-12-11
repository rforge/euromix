\name{paraMix}
\alias{paraMix}
\alias{checkInput} 

\title{
Likelihood for mixtures with related contributors based on \code{\link{paramlink}}
}
\description{
A DNA mixture (R) has been observed and some individuals may have been typed.
Some of these typed individuals are known contributors to the mixture, some are known non-contributors.
In addition, there may be specified untyped individuals that have contributed to the mixture.
Individuals can be specified as members of a pedigree defined by a  \code{\link{linkdat}} object \code{x}
corresponding to a hypothesis \code{H}. Relevant individuals unrelated to all others, are defined using
\code{\link{singleton}}.The likelihood 

\code{Pr(mixture,Typed contributors,Typed non-contributors|H)=P(R,T,V|H)}

is calculated; the notation on the right hand side corresponds to that of Curran, Gill and Bill (2005).
A plot is also produced summarising the essential information.
Compared to previous literature and methods, including a series of papers by Fung and Hu, we generalise calculations to allow for
general, possibly inbred, pedigrees. 
Typically calculations are performed for competing hypotheses and the ratio of likelihoods, 
the likelihood ratio \code{LR} is calculated and reported. Previous methods have assumed the relationships 
between typed contributors to be same for the competing hypotheses. This restriction does not apply for
our approach. The calculation may also be used for identification cases where a mixture and reference samples are available.
Likelihood calculations are performed using the \code{\link{likelihood}} of \code{\link{paramlink}}.
The function \code{checkInput} checks input to \code{paraMix}.
}
\usage{
paraMix(x, R, id.U, id.V = NULL, alleles, afreq = NULL, 
Xchrom= FALSE, known_genotypes = list(), loop_breakers =NULL, 
eliminate = 0, check = TRUE, plot = TRUE, title= NULL)
checkInput(x, R, id.U, id.V, alleles, all_typed, K, R_not_masked)
}

\arguments{
  \item{x}{
\code{\link{linkdat}} object, or a list of such (if disconnected), describing the claimed relationship.
}
  \item{R}{
Integers, mixture.
}
  \item{id.U}{
Integers indicating untyped contributors (e.g.,suspect(s)).
}
  \item{id.V}{
Integers indicating typed non-contributors.
}
 \item{alleles}{
Integers indicating alleles for marker. 
}
  \item{afreq}{
A numerical vector with allele frequencies. An error is given if they don't sum to 1 (rounded to 3 decimals).
}
  \item{Xchrom}{
Logical, \code{FALSE} for autosomal marker.
}
  \item{known_genotypes}{
List, each element a triplet of integers corresponding to (id,allele1,allele2) 
}
  \item{loop_breakers}{
A numeric containing IDs of individuals to be used as loop breakers. Relevant only if the pedigree has loops. See \code{\link{breakLoops}}.
}

 \item{eliminate}{
A non-negative integer, indicating the number of iterations in the internal genotype-compatibility algorithm. 
Positive values can save time if partialmarker is non-empty and the number of alleles is large.
}
 \item{check}{
If TRUE check of input is performed and calculations stop if they are likely to take too much time
}
 \item{plot}{
If TRUE a plot is produced
}
 \item{title}{
 Title of the plot
 }
 \item{all_typed}{
An integer vector identifying typed individuals
}
\item{K}{
Known alleles in contrib_typed
}
\item{R_not_masked}{
Unexplained alleles
}
}
\value{
\item{likelihod }{The likelood Pr(R,T,V|H)}
\item{allLikelihoods }{Terms adding to above Pr(R,T,V|H)}
}

\author{
Magnus Dehli Vigeland and Thore Egeland <Thore.Egeland@nmbu.no>
}

\seealso{
\code{\link{famMix}}
}
\details{
The required likelihood 
\code{Pr(R,T,V|H)=Pr(R|T,V,H)Pr(T,V|H)= Pr(T,V|H)sum_u Pr(U=u,T,V|H)}
where the sum extends over u among persons specified by id.U so that the union of u,T, V is R.
The likelihoohd for each u and the sum is returned.
Assumes alleles to be numbered 1,2,...
}

\examples{
#Example 1: Motivating example Egeland et al. (2013)
require(paramlink)
y1=swapSex(nuclearPed(3),c(3,4))
p=c(0.1,0.2,0.3,0.4)
alleles=1:length(p)
T1=c(1,1)
T2=c(2,2)
R=1:2
known=list(c(3,T1),c(4,T2))
l1=paraMix(y1,R,id.U=5,alleles=alleles,afreq=p,known_genotypes=known)
y2=swapSex(nuclearPed(1),3)
y2=addOffspring(y2,mother=2,noff=1,sex=2)
y2=relabel(y2,c(1:3,6,4),1:5)
l2=paraMix(y2,R,id.U=6,alleles=alleles,afreq=p,known_genotypes=known)
LR1=l1$lik/l2$lik
exact=1/(2*(p[1]+p[2]))
stopifnot(abs(LR1-exact)<10^(-6)) 

#Example 2. Example 1 in Egeland et al. (2013) based on Fung and Hu (2008)
#Data:
#Mixture 1/2/3
#Suspect=4, genotype 3/3
#Victim=10, genotype 1/2
#H1: Contributors were the suspect  and victim (unrelated)
#H2: Contributors were the father of suspect  and victim (unrelated)
#H3: Contributors were the brother of suspect  and victim (unrelated)
afreq=c(0.044,0.166,0.110,0.680)
alleles=1:length(afreq)
R=1:3 #Mixture
man_ped=nuclearPed(2)
victim = singleton(id=10, sex=2)
known = list(c(4,3,3),c(10,1,2)) #individual 4 is 3/3, and 10 (the victim) is 1/2.
#The likelihoods corresponding to H1,H2 and H3
l1=paraMix(list(man_ped, victim), R, id.U=NULL, id.V=NULL, 
alleles=alleles, afreq=afreq, known_genotypes=known)$lik
l2=paraMix(list(man_ped, victim), R, id.U=1, id.V=4, 
alleles=alleles, afreq=afreq, known_genotypes=known)$lik
l3=paraMix(list(man_ped, victim), R, id.U=3, id.V=4, 
alleles=alleles, afreq=afreq, known_genotypes=known)$lik
LR12=l1/l2
stopifnot(abs(LR12-3.125)<10^(-6)) 
LR13=l1/l3
stopifnot(abs(LR13- 2.355296)<10^(-6)) 
}

