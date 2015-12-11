\name{famMix}
\alias{famMix}
\title{
Likelihood for mixtures with related contributors based on \code{\link{Familias}}
}
\description{
Likelihood for mixtures with related contributors based on \code{\link{Familias}}. For a general description
of the problem, see \code{\link{paraMix}}. As opposed to \code{\link{paraMix}} this function uses the R version
of \code{\link{Familias}} for likelihood calculation and therefore theta-correction, mutation models and silent allele frequencies
(but not X-chromosomes or simulation)
are accomodated. 
}
\usage{
famMix(x, R, id.U, id.V = NULL, partialmarker = NULL, theta = 0, 
mutationRateFemale = 0, mutationRateMale = 0, 
mutationModelFemale = "stable", mutationModelMale = "stable", 
mutationRangeFemale = 0.1, mutationRangeMale = 0.1, 
silentFrequency = 0,check=TRUE)
}
\arguments{
  \item{x}{
\code{\link{linkdat}} object.
}
  \item{R}{
Integers, mixture.
}
  \item{id.U}{
List of unknown contributors (e.g.,suspect(s)).
}
  \item{id.V}{
Integers indicating typed non-contributors.
}
  \item{partialmarker}{
A \code{\link{marker}} object. 
}
  \item{theta}{
Real in [0,1]
}
  \item{mutationRateFemale}{
See \code{\link{FamiliasLocus}}.
}
  \item{mutationRateMale}{
See \code{\link{FamiliasLocus}}.
}
  \item{mutationModelFemale}{
See \code{\link{FamiliasLocus}}.
}
  \item{mutationModelMale}{
See \code{\link{FamiliasLocus}}.
}
  \item{mutationRangeFemale}{
See \code{\link{FamiliasLocus}}.
}
  \item{mutationRangeMale}{
See \code{\link{FamiliasLocus}}.
}
  \item{silentFrequency}{
Real in [0,1].
}
 \item{check}{
If TRUE check of input is performed and calculations stop if they are likely to take too much time.
}
}

\details{
See \code{\link{paraMix}}.
}

\value{
\item{x}{ \code{\link{linkdat}} object updated with genotypes of missing individuals specified by \code{id.U}}
\item{likelihod }{The likelood Pr(R,T,V|H)}
\item{allLikelihoods }{Terms adding to above Pr(R,T,V|H)}
}

\references{
Egeland et al (2013)
}
\author{
Thore Egeland <Thore.Egeland@nmbu.no>
}

\seealso{
\code{\link{paraMix}}
}
\examples{
#Example
require(paramlink)
require(Familias)

% #Example 1. Motivating example in Egeland, Dorum and Sheehan (2013). See
% #paraMix for another implementation
% theta=0.03
% y1=swapSex(nuclearPed(3),c(3,4))
% afreq=c(0.1,0.2,0.3,0.4)
% afreq=afreq/sum(afreq)
% alleles=1:length(afreq)
% T1=c(1,1)
% T2=c(2,2)
% R=1:2
% m1=marker(y1,3,T1,4,T2,alleles=alleles,afreq=afreq)
% par(mfcol=c(1,2))
% plot(addMarker(y1,m1),marker=1,title="Mixture=1/2;id.U=5;id.V=None")
% l1=famMix(y1,R,5,partialmarker=m1,theta=theta)
% y2=swapSex(nuclearPed(1),3)
% y2=addOffspring(y2,mother=2,noff=1,sex=2)
% m2=marker(y2,3,T1,5,T2,alleles=alleles,afreq=afreq)
% plot(addMarker(y2,m2),marker=1,title="Mixture=1/2;id.U=4;id.V=None")
% par(mfcol=c(1,1))
% l2=famMix(y2,R,4,partialmarker=m2,theta=theta)
% LR1=l1$lik/l2$lik
% pp=sum(afreq[1:2])
% exact=(1+3*theta)/(4*theta+pp*(1-theta))*(1/2)
% stopifnot(abs(LR1-exact)<10^(-6))


% #Example 2. 
% #The hypotheses are
% #Hp: contributors were the victim and the suspect ((unrelated)):
% #Hd: contributors were the victim and a relative of the suspect
% #We consider the marker FGA with alleles (original denomination  and
% #frequencies in paranthesis)
% # 1 (20, 0.044), 2 (24,0.166), 3 (25,0.11), 4( rest,0.680)
% #The victim is 1/2, the suspect 3/3.
% #The mixture is 1/2/3
% #First the pedigree is defined to include to unrelated persons
% #(which happens to appear as married; the points is that they are unrelated)

% yp=nuclearPed(1)
% afreq=c(0.044,0.166,0.110,0.680)
% alleles=1:length(afreq)
% R=1:3
% m1=marker(yp,2,c(1,2),1,c(3,3),alleles=alleles,afreq=afreq)
% par(mfcol=c(1,2))
% plot(addMarker(yp,m1),marker=1,id.labels=c("Susp(1)","Vict(2)","Extra (3)"),title="Hp. Fung&Hu, Table 7.3, FGA")
% lp=famMix(yp,R,id.U=NULL,partialmarker=m1) 

% #Consider next Hd with the contributor being the father of the suspect
% yd=nuclearPed(1)
% yd=swapSex(addOffspring(yd,father=1,noff=1),5)
% m2=marker(yd,4,c(1,2),3,c(3,3),alleles=alleles,afreq=afreq)
% plot(addMarker(yd,m2),marker=1,title="Hd. Relative:Father of suspect.Fung&Hu, Table 7.3, FGA")
% par(mfcol=c(1,1))
% ld1=famMix(yd,R,id.U=1,id.V=3,partialmarker=m2)
% LR1=lp$likelihood/ld1$likelihood 
% stopifnot(abs(LR1-3.125)<10^(-6))
% # Next example is extended to allow for theta correction confirming Table 7.9 Fung&Hu (2008)
% theta=0.01
% (LR1=famMix(yp,R,id.U=NULL,partialmarker=m1,theta=theta)$lik/famMix(yd,R,id.U=1,id.V=3,partialmarker=m2,theta=theta)$lik)
% stopifnot(abs(LR1-2.886771)<10^(-6))
% # Next example is extended to allow for silent allele
% silentFrequency=0.02
% l1=famMix(yp,R,id.U=NULL,partialmarker=m1,theta=theta,silentFrequency=silentFrequency)$lik
% l2=famMix(yd,R,id.U=1,id.V=3,partialmarker=m2,theta=theta,silentFrequency=silentFrequency)$lik
% LR1=l1/l2
% stopifnot(abs(LR1-3.045986)<10^(-6))
% # Next example is extended to allow for mutation. Default model with mutation rate 0.004.
% mutationRateFemale=mutationRateMale=Rmut=0.004
% l1=famMix(yp,R,id.U=NULL,partialmarker=m1,theta=theta,mutationRateFemale=Rmut,mutationRateMale=Rmut)$lik
% l2=famMix(yd,R,id.U=1,id.V=3,partialmarker=m2,theta=theta,mutationRateFemale=Rmut,mutationRateMale=Rmut)$lik
% LR1=l1/l2
% stopifnot(abs(LR1-2.917555)<10^(-6))

% #Example 3. Mutation
% #The mother and child, individuals 2 and 3 are knownto be contributors.
% #The question is wheteher the father (1 below also is a contributor).
% #In this case, mutation does not matter as the familiy relationships is the same for both hypotheses
% xp=swapSex(nuclearPed(1),3)
% alleles=1:4
% afreq=c(0.1,0.2,0.3)
% p=c(afreq,1-sum(afreq))
% T1=c(1,1);T2=c(3,3);T3=c(2,3)
% mp=marker(xp,1,T1,2,T2,3,T3,alleles=alleles,afreq=p)
% par(mfcol=c(1,2))
% plot(addMarker(xp,mp),marker=1)
% R=1:3
% mu=0.004;mod="stable"
% lik.HP=famMix(xp,R,id.U=NULL,partialmarker=mp,mutationRateFemale =mu, mutationRateMale = mu,mutationModelFemal=mod,mutationModelMale=mod)
% xd=addOffspring(xp,mother=2,noff=1)
% md=marker(xd,1,T1,2,T2,3,T3,alleles=alleles,afreq=p)
% plot(addMarker(xd,md),marker=1)
% lik.HD=famMix(xd,R,id.U=4,id.V=1,partialmarker=md,mutationRateFemale =mu, mutationRateMale = mu,mutationModelFemal=mod,mutationModelMale=mod)
% LR.D18=lik.HP$likelihood/lik.HD$likelihood
% stopifnot(abs(LR.D18-1/(p[1]^2+2*p[1]*p[2]+2*p[1]*p[3]))<10^(-6))
% par(mfcol=c(1,1))
}

