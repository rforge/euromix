\name{convertToFamilias}
\alias{convertToFamilias}
\title{
Convert genotype data to Familias format
}
\description{
Genotype data are transformed from lines two columns.
If there is only one column for each marker and two lines for each individual,
the data is transformed so that there is one line for each individual.
The pecularities of the input format of Familias is handled.
}
\usage{
convertToFamilias(infile, outfile = paste("out", infile, sep = ""))
}

\arguments{
  \item{infile}{
File name.
}
  \item{outfile}{
File name for output file.
}
}
\details{
The first column is the name of the individual, the second indicates sex (X X for females,
X Y for males) while the
remaining columns are marker names (no blanks are allowed in names;
usual rules for variable names apply). There are two lines for each individual.
The smallest example (below) is for a female called 32293 with genotypes 15/16 for Marker1:
\tabular{rrr}{
         Name  \tab Sex   \tab Marker \cr
         32293 \tab X     \tab 15 \cr
         32293 \tab X     \tab 16 \cr
       }
Typically, Familias is started by loading a file containing the database (markers and allele
frequencies). This file needs only be prepared once for each individual. The 
\code{Case Related DNA Data} can then be read from the file produced by \code{convertToFamilias}.
Note that allele must have precisely the same name (8 and 8.0 are different alleles for instance)
}
\value{
A file with default file name as input file name preceded by "out".
This file can be read by Familias in the \code{Case Related DNA Data} window.
}

\author{
Thore.Egeland@nmbu.no
}

\examples{
\dontrun{convertToFamilias("denise3.txt")} 
}
