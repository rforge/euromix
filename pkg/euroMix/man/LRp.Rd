\name{LRp}
\alias{LRp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the p-value corresponding to a likelihood ratio.
}
\description{
Computes the likelihood ratio for the given hypotheses and finally calculates a p-value corresponding to the likelihood ratio. The p-value is the probability of observing a likelihood ratio at least as large as the one observed, given that the defense hypothesis is true.
}
\usage{
LRp( sampleData, victimData, suspectData, db, hp, hd, prD, prC )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleData}{
Data frame or matrix with sample profile. Each column represent an allele, each row represent a marker. Only autosomal markers. Marker names that correspond with markers in allele frequency database must be given as row names.
}
 \item{victimData}{
Data frame or matrix with victim profile. Each column represent an allele, each row represent a marker. Only autosomal markers. Markers must be in the same order as for sampleData.
}
 \item{suspectData}{
Data frame or matrix with suspect profile. Each column represent an allele, each row represent a marker. Only autosomal markers. Markers must be in the same order as for victimData and sampleData.
}
\item{db}{
Data frame with allele frequencies. Data for the various markers are stacked. First column contains marker names, each name repeated as many times as there are alleles for the marker. Second column contains the allele names and third column contains the frequencies. 
}
\item{hp}{
Prosecution hypothesis. A character vector of all contributors under $H_p$, where S denotes suspect, V victim and U unknown. E.g. if the hypothesis is that the sample is a mixture of the suspect, the victim and one unknown, this is specified with the vector c('S','V,'U').
}
\item{hd}{
Defense hypothesis. A character vector of all contributors under $H_d$, specified like hp. E.g. if the hypothesis is that the sample is a mixture of the suspect, and two unknowns, this is specified with the vector c('S','U,'U').
}
\item{prD}{
Probability of drop-out. A number between 0 and 1.
}
\item{prC}{
Probability of drop-in. A number between 0 and 1.
}
}
\details{
The function is a wrapper for \link{pvalue.machine}. Likelihood ratios are computed with the \link{LR} function in \link{forensim}. Use \link{pvalue.machine} for a more generic function that is independent of LR model.
}
\value{
\item{LR}{Likelihood ratio}
\item{pvalue}{P-value corresponding to the likelihood ratio}
}
\references{
Dorum et al. Exact computation of the distribution of likelihood ratios with forensic
  applications. FSI: Genetics, 9, 2014, doi: http://dx.doi.org/10.1016/j.fsigen.2013.11.008
}
\author{
Guro Dorum <guro.dorum@nmbu.no>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pvalue.machine},\link{LRpvalue}}
}

\examples{
data(R,S,V)
data(db2)
LRp(sampleData=R,victimData=V,suspectData=S,db=db2,hp=c('V','S'),hd=c('V','U'),prD=0.47,prC=0.05 )
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
