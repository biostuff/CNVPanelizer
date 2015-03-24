\name{CombinedNormalizedCounts}
\alias{CombinedNormalizedCounts}
\title{
  CombinedNormalizedCounts
}
\description{
  This function makes use of cn.mops::normalizeChromosomes to normalize the read counts of all samples and references
  to the same median read count.
}
\usage{
CombinedNormalizedCounts(sampleCounts, referenceCounts, ampliconNames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleCounts}{
      Vector or matrix  with sample read counts (rows: amplicons, columns: samples)
  }
    \item{referenceCounts}{
      Matrix with reference read counts (rows: amplicons, columns: samples)
    }
    \item{ampliconNames}{
      A vector with the amplicon names used for the reference and sample matrices
  }
}

\value{
A list object with
\item{comp1 }{The  normalized sample matrix/vector}
\item{comp2 }{the median normalized read count matrix of the reference samples}
}

\author{
  Thomas Wolf, Cristiano Oliveira
}


\examples{

%\dontshow{
%}
data(sampleReadCounts)
data(referenceReadCounts)

normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts,
                                                 referenceReadCounts)
}
