\name{Background}
\alias{Background}
\title{
  Background
}
\description{
    Estimates the dependence of background noise on the number of amplicons
    sequenced for each gene.
}
\usage{
Background(geneNames,
           samplesNormalizedReadCounts,
           referenceNormalizedReadCounts,
           bootList,
           replicates = 1000)
}
\arguments{
  \item{geneNames}{
    Names of the genes
}
  \item{samplesNormalizedReadCounts}{
  A matrix with the normalized read counts of the samples of interest
}
  \item{referenceNormalizedReadCounts}{
  A matrix with the normalized reference read counts
}
  \item{bootList}{
  A list as returned by BootList
}
  \item{replicates}{
    an integer number of how many replicates should be performed
}
}

\value{
 Returns a list of data frames. One data frame for each sample of interest.
 The data frames report the background noise for each number of gene
 associated amplicons
 being sequenced.
}

\author{
    Thomas Wolf,
    Cristiano Oliveira 
}

\examples{

data(sampleReadCounts)
data(referenceReadCounts)
## Gene names should be same size as row columns
geneNames <- row.names(referenceReadCounts)

ampliconNames <- NULL

normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts,
                                                 referenceReadCounts,
                                                 ampliconNames = ampliconNames)

# After normalization data sets need to be splitted again to perform bootstrap
samplesNormalizedReadCounts = normalizedReadCounts["samples"][[1]]
referenceNormalizedReadCounts = normalizedReadCounts["reference"][[1]]

#Values above 10000 should be used
replicates <- 10

# Perform the bootstrap based analysis
bootList <- BootList(geneNames,
                     samplesNormalizedReadCounts,
                     referenceNormalizedReadCounts,
                     replicates = replicates)

background <- Background(geneNames,
                        samplesNormalizedReadCounts,
                        referenceNormalizedReadCounts,
                        bootList,
                        replicates = replicates)
}

