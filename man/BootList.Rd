\name{BootList}
\alias{BootList}
\title{
  BootList
}
\description{
    Performs a hybrid bootstrapping subsampling
    procedure similar to random forest. It bootstraps the reference
    samples and subsamples the amplicons associated with each gene.
    Returns a distribution of sample/reference ratios for
    each gene and sample of interest combination.
}
\usage{
    BootList(geneNames, sampleMatrix, refmat, replicates)
}
\arguments{
  \item{geneNames}{
     A vector of gene names, with one entry for each sequenced
     amplicon.
}
  \item{sampleMatrix}{
     A vector or matrix of the read counts from the sample of interest.
     In the case of a matrix columns represent samples and rows amplicons.
}
  \item{refmat}{
     A matrix of the read counts obtianed from the reference samples.
     Columns represent reference samples and rows amplicons.
}
  \item{replicates}{
     How many bootstrap replicates should be performed.
}

}

\value{

 Returns a list of numeric matrices:
 For each matrix a row represent a gene while each column represents
 a bootstrapping/subsampling iteration.

}

\author{
  Thomas Wolf, Cristiano Oliveira
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

# Should be used values above 10000
replicates <- 10

# Perform the bootstrap based analysis
bootList <- BootList(geneNames,
         samplesNormalizedReadCounts,
         referenceNormalizedReadCounts,
         replicates = replicates)
}
