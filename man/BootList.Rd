\name{BootList}
\alias{BootList}
\title{
  BootList
}
\description{
    Performs a hybrid bootstrapping subsampling procedure similar to random forest.
    A distribution of sample/reference ratios for each gene and sample of interest combination.

}
\usage{
BootList(geneNames, sampleMatrix, refmat, reps, refWeights)
}
\arguments{
  \item{geneNames}{
     A vector of the gene names used in the panel.
}
  \item{sampleMatrix}{
     A vector or matrix of the read counts from the sample of interest.
     In the case of a matrix columns represent samples and rows amplicons.
}
  \item{refmat}{
     A matrix of the read counts obtianed from the reference samples.
     Columns represent reference samples and rows amplicons.
}
  \item{reps}{
     How many bootstrap replicates should be performed.
}
  \item{refWeights}{
     A numeric vector that allows for the weighting of reference samples.
     Reference samples with a higher weight will be sampled more often.
}
}

\value{

 Returns a list of numeric matrices:
 For each matrix a row represent a gene while each column represents
 bootstrapping sumsampling iteration.

}

\author{
  Thomas Wolf, Cristiano Oliveira
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
%% TODO create some

\dontrun{
  
   bootList = BootList(genesPositionsIndex,
                    samplesNormalizedReadCounts,
                    referenceNormalizedReadCounts,
                    reps = replicates,
                    refWeights = NULL)

}


}
