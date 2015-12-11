\name{LRpvalue}
\alias{LRpvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the p-value corresponding to a likelihood ratio.
}
\description{
Reads mixture data from files, computes the likelihood ratio for the given hypotheses and finally calculates a p-value corresponding to the likelihood ratio. The p-value is the probability of observing a likelihood ratio at least as large as the one observed, given that the defense hypothesis is true.
}
\usage{
LRpvalue(samplefile, victimfile, suspectfile, freqfile, hp, hd, prD, prC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samplefile}{
CSV file with sample profile. The file can only contain data for autosomal markers, and apart from that the format is the same as required in LRmix. See the LRmix manual for details. The file name must contain the complete path if the file is not in the current working directory.
}
\item{victimfile}{
CSV file with victim profile. Same format as in LRmix. Only autosomal markers.
}
\item{suspectfile}{
CSV file with suspect profile. Same format as in LRmix. Only autosomal markers.
}
\item{freqfile}{
CSV file with allele frequencies. Same format as in LRmix.
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
The function is a wrapper for \code{\link{LRp}} which again is a wrapper for \code{\link{pvalue.machine}}. Likelihood ratios are computed with the \code{\link{LR}} function in \code{\link{forensim}}. For liberties regarding reading data from files, \code{\link{LRp}} can be used to compute p-values for already prepared data frames. For liberty also regarding the LR model used, \code{\link{pvalue.machine}} is the most generic function to compute a p-value. 
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
\code{\link{pvalue.machine},\link{LRp}}
}
\examples{
data(sample);data(suspect);data(victim);data(freqs)
samplefile <- tempfile(); write.table(sample, samplefile, sep=",", row.names=FALSE)
victimfile <- tempfile(); write.table(victim, victimfile, sep=",", row.names=FALSE)
suspectfile <- tempfile(); write.table(suspect, suspectfile, sep=",", row.names=FALSE)
freqfile <- tempfile(); write.table(freqs, freqfile, sep=",", row.names=FALSE)

LRpvalue(samplefile, victimfile, suspectfile, freqfile, hp=c("V","S"), hd=c("V","U"), 
prD=0.47, prC=0.05)
unlink(c(samplefile, victimfile, suspectfile, freqfile))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
