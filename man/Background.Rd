\name{Background}
\alias{Background}
\title{
  Background
}
\description{
  Estimates the dependence of background noise on the number of amplicons sequenced for each gene.
}
\usage{
Background(geneNames,
           samplesNormalizedReadCounts,
           referenceNormalizedReadCounts,
           sigList,
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
  \item{sigList}{
  a list as returned by CheckSignificance
}
  \item{replicates}{
    an integer number of how many replicates should be performed
}
}

\value{
 Returns a list of data frames. One data frame for each sample of interest.
 Thr data frames report the background noise for each number of gene associated amplicons
 being sequenced.
}

\author{
  Cristiano Oliveira, Thomas Wolf
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

\dontrun{

Background(geneNames,
           samplesNormalizedReadCounts,
           referenceNormalizedReadCounts,
           sigList,
           replicates = 1000)

}
}

