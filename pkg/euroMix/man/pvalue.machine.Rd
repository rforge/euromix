\name{pvalue.machine}
\alias{pvalue.machine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the p-value for LR.suspect
}
\description{
It is difficult to obtain accurate p-values based on simulation. This function provides
an exact alternative.
}
\usage{
pvalue.machine(LR.suspect, LR.table, P.table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LR.suspect}{
Numeric. Observed likeliood ratio (1x1 positive value) 
}
  \item{LR.table}{
Pre-computed likelihood ratios for every genotype of every marker (MxG matrix).
Each row corresponds to a marker. G is the maximum number of genotypes for any marker. 
Markers with fewer than G genotypes must have 0 in redundant columns
}
  \item{P.table}{
The population probabilities for every genotype of every marker (MxG matrix).
Must corresponds to the genotypes in LR.table. See description of LR.table
}
}

\value{
The p-value, where a value close to 0 indicates that the suspect is a contributor.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Dorum, Bleka, Snipen <guro.dorum@nmbu.no>
}

\seealso{
The function is obsolete. 

See
\code{\link[DNAprofiles]{dists.product}} and \code{\link[DNAprofiles]{dists.product.pair}} for efficient computation of likelihood ratio distributions.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
#Simple example, 2 markers, 3 genotypes. LR's and genotype probabilities precalculated
#The LR's for all possible genotypes for both markers. Each row corresponds to a marker
LR.table <- matrix(c(6,5,5,4,3,2),2,3)
#The population probabilities corresponding to the genotypes in LR.table
P.table <- rbind(c(0.2, 0.4, 0.4), c(0.1,0.6,0.3))
#LR observed for suspect
LR.suspect <- 20
pvalue <- pvalue.machine(LR.suspect, LR.table, P.table)
cat("p-value = ", pvalue, "\n")

% ##########################################################################
% #Example with 9 markers. Compute LR's and gentoype probabilities first
% require(forensim)
% data(db2)
% data(ex9m)
% M <- nrow(R)
% markers <- rownames(R)

% #LR for Hp: victim + suspect vs. Hd: victim + unknown
% lr <- A <- numeric()
% for(i in 1:M){
  % r <- as.numeric(R[i,!is.na(R[i,])])
  % v <- as.numeric(V[i,])
  % s <- as.numeric(S[i,])
  % afreq <- db2[db2[,1]==markers[i],3]
  % A[i] <- choose(length(afreq)+1,2)
  % names(afreq) <- db2[db2[,1]==markers[i],2]
  % lr[i] <- LR(r, Tp=c(v,s), Td=v, Vp=NULL, Vd=s, 
            % xp=0, xd=1, theta=0, prDHet=rep(0.47,4), prDHom=rep(0.47,4)^2, prC=0.05, freq=afreq)$LR
% }
% LR.suspect <- prod(lr)

% #LR when replacing suspect with all possible genotypes (i.e. all random men)
% LR.table <- P.table <- matrix(0,M,max(A))
% for(i in 1:M){
  % afreq <- db2[db2$Marker==markers[i],]$Frequency
  % names(afreq) <- db2[db2$Marker==markers[i],]$Allele
  % #All possible genotypes for the given marker
  % gtAll <- as.matrix(expand.grid(names(afreq),names(afreq)))
  % gtAll <- gtAll[gtAll[,1]<=gtAll[,2],]
  % #likelihood ratios and corresponding genotype probabilities
  % lrAll <- pAll <- numeric()
  % for(j in 1:nrow(gtAll)) {
    % rm <- gtAll[j,]
    % r <- as.numeric(R[i,!is.na(R[i,])])
    % v <- as.numeric(V[i,])
    % s <- as.numeric(S[i,])
    % lrAll[j] <- LR(r, Tp=c(v,rm), Td=v, Vp=NULL, Vd=rm, 
                % xp=0, xd=1, theta=0, prDHet=rep(0.47,4), prDHom=rep(0.47,4)^2, prC=0.05, freq=afreq)$LR
    % pAll[j] <- ifelse(rm[1]==rm[2], prod(afreq[as.character(rm)]), 2*prod(afreq[as.character(rm)]) )
  % }
  % #Sort LR's in decreasing order
  % idxSort <- order(lrAll, decreasing=TRUE)
  % LR.table[i,1:length(lrAll)] <- lrAll[idxSort]
  % P.table[i,1:length(lrAll)] <- pAll[idxSort]
% }

% #The p-value
% pvalue <- pvalue.machine(LR.suspect, LR.table, P.table)
% cat("p-value = ", pvalue, "\n")

}