
pvalue.machineC <- function( LR.suspect, LR.table, P.table ){
#
# Computes the p-value for LR.suspect
#
# Input:
# LR.suspect  - The likelihood ratio for the suspect genotype profile (1x1 positive value)
# LR.table    - The pre-computed likelihood ratios for every genotype of every marker (MxG matrix).
#               Each row corresponds to a marker. G is the maximum number of genotypes for any marker. 
#               Markers with fewer than G genotypes must have 0 in redundant columns.
# P.table     - The population probabilities for every genotype of every marker (MxG matrix)
#               Must corresponds to the genotypes in LR.table. See description of LR.table.
  
  #sys <- Sys.info()[1]
  #if(sys=="Windows") file <- system.file("files/prodrecfunction.dll", package = "euroMix")
  #if(sys=="Darwin") file <- system.file("files/prodrecfunction.so", package = "euroMix")
  #else stop("Unknown system")
  #dyn.load(file)
  
  # The number of markers
  M <- nrow( LR.table )
  
  # Preparing input
  prepped <- prepare.input( LR.table, P.table )
  LR.list <- prepped$LR.list
  Plist  <- as.numeric(unlist(prepped$P.list))
  inflationfactor <- as.numeric(prepped$Inflation.factor)
  
  # The number of genotypes for each marker
  G <- as.integer(sapply( LR.list, length ))

  # Initial values
  pvalue <- as.numeric(1)
  LRsuspect  = as.numeric(LR.suspect)
  LRlist <- as.numeric(unlist(LR.list))
 
  p.value <- .C("prodrecfunction",pvalue, LRsuspect, LRlist,Plist,M,G,inflationfactor)[[1]]
  #p.value <- .C("prodrecfunction",pvalue, LRsuspect, LRlist,Plist,M,G,inflationfactor,PACKAGE="euroMix")[[1]]
  return( p.value )
}