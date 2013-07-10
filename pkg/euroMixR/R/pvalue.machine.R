


pvalue.machine <- function( LR.suspect, LR.table, P.table ){
#
# Computes the p-value for LR.suspect
#
# Input:
# LR.suspect  - The likelihood ratio for the suspect genotype profile (1x1 positive value)
# LR.table    - The pre-computed likelihood ratios for every genotype of every marker (MxG matrix)
# P.table     - The population probabilities for every genotype of every marker (MxG matrix)
  
  # The number of markers
  M <- nrow( LR.table )
  
  # Preparing input
  prepped <- prepare.input( LR.table, P.table )
  LR.list <- prepped$LR.list
  P.list  <- prepped$P.list
  inflation.factor <- prepped$Inflation.factor
  
  # The number of genotypes for each marker
  G <- sapply( LR.list, length )

  # Initial values
  marker <- 1
  score <- 1.0
  pvalue <- 1.0
  
  p.value <- recfun( marker, score, pvalue, LR.suspect, LR.list, P.list, M, G, inflation.factor )
  return( p.value )
}


recfun <- function( marker, score, pvalue, LR.suspect, LR.list, P.list, M, G, inflation.factor ){
#
# The recursive function that computes the p-value.
#
  
  if( marker == M ){ # The terminating case, we are at the last marker
    pvals <- rep( 0, G[marker] ) # There are G[m] genotypes of this marker, all must have 0 p-value until we know they give score above LR.suspect
    for( genotype in 1:G[marker] ){ # Looping over all genotypes, compute LR-score, and p-value if score is above LR.suspect
      if( score * LR.list[[marker]][genotype] >= LR.suspect ){
        pvals[genotype] <- pvalue * P.list[[marker]][genotype]
      } else {
        break; # Breaking out of the loop, no need to compute more scores since genotypes are sorted in descending order
      }
    }
    return( sum( pvals ) ) # This is the p-value for the path leading here
  } else { # We are not yet at the final marker
    pvals <- rep( 0, G[marker] ) # There are G[m] genotypes of this marker, all must have 0 p-value until we know they give score above LR.suspect
    for( genotype in 1:G[marker] ){
      score.new <- score * LR.list[[marker]][genotype]
      pvalue.new <- pvalue * P.list[[marker]][genotype]
      if( score.new * inflation.factor[marker+1] >= LR.suspect ){ # Give up further recursions if we cannot inflate current score above LR.suspect
        pvals[genotype] <- recfun( marker+1, score.new, pvalue.new, LR.suspect, LR.list, P.list, M, G, inflation.factor ) # We dig one marker deeper into the tree...
      } else {
        break; # Breaking out of loop and terminating recursion because we cannot inflate current score
      }
    }
    return( sum( pvals ) ) # This is the p-value for the path leading here
  }
}



prepare.input <- function( LR.table, P.table ) {
#
# Converting the tables of likelihood ratios and probabilities to lists,
# keeping only the unique LR-scores for each marker, and sorting everything
# in proper order.
#
# Input:
# LR.table  - Matrix of likelihood ratio scores for all genotypes of all markers. 
#             Each row correspond to a marker, and it may contain 0's (MxG matrix)
# P.table   - Matrix of population probabilities for all genotypes of all markers.
#             Each row correspond to a marker, and it may contain 0's (MxG matrix)
#
# The function returns a list containing two lists, LR.list and P.list, containing the unique 
# non-zero values in LR.table and P.table. For each marker the elements are sorted in descending
# order by LR-score.
#
  M <- nrow( LR.table ) # The number of markers
  LR.list <- P.list <- vector( "list", M )
  LR.max <- apply( LR.table, 1, max )
  
  # Smart sorting of markers
  ix <- order( LR.max, decreasing=T )
  LR.table <- LR.table[ix,]
  P.table <- P.table[ix,]
  LR.max <- LR.max[ix]
  
  for( marker in 1:M ){
    # Finding unique likelihood ratios and summing probabilities accordingly
    lr <- unique( LR.table[marker,] )
    pr <- tapply( P.table[marker,], match( LR.table[marker,], lr ), sum )
    # Discarding 0's
    ix <- which( lr > 0 )
    lr <- lr[ix]
    pr <- pr[ix]
    # Sorting in descending order by LR
    ix <- order( lr, decreasing=T )
    LR.list[[marker]] <- lr[ix]
    P.list[[marker]] <- pr[ix]
  }
  Inflation.factor <- cumprod( LR.max[M:1] )[M:1]
  return( list( LR.list=LR.list, P.list=P.list, Inflation.factor=Inflation.factor ) )
}
