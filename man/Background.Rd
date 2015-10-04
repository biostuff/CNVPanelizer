\name{Background}
\alias{Background}
\title{
  Background
}
\description{
    Makes use of a subsampling approach to estimate the background
    noise when sequencing a gene with a specific number of amplicons.
    The 95 percent confidence interval is returned for each unique number
    of amplicons in the experiment.
}
\usage{
Background(geneNames,
           samplesNormalizedReadCounts,
           referenceNormalizedReadCounts,
           bootList,
           replicates = 1000,
           significanceLevel = 0.05,
           robust = FALSE)
}
\arguments{
  \item{geneNames}{
    A vector of gene names, with one entry for each sequenced
    amplicon.
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
  \item{significanceLevel}{
    The significance level for the calculated confidence interval
}
  \item{robust}{
    If set to true the confidence interval is calculated replacing mean
    with median and sd with mad.
}
}

\value{
 Returns a list of data frames. One data frame for each sample of interest.
 The data frames report the 95 percent confidence interval of the background
 noise for each number of amplicons and sample combination.

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
                        replicates = replicates,
                        significanceLevel = 0.1)
}

