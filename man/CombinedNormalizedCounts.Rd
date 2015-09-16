\name{CombinedNormalizedCounts}
\alias{CombinedNormalizedCounts}
\title{
  CombinedNormalizedCounts
}
\description{
    This function makes use of NOISeq::tmm to normalize
    the read counts of all samples and references to the same median read count
}
\usage{
    CombinedNormalizedCounts(sampleCounts,
                            referenceCounts,
                            ampliconNames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{sampleCounts}{
         Matrix or vector with sample read counts (rows: amplicons, columns: samples)
    }
    \item{referenceCounts}{
        Matrix with reference read counts (rows: amplicons, columns: samples)
    }
    \item{ampliconNames}{
        A vector with  amplicon defining names for the reference and sample matrices
    }
}

\value{
A list object with  two matrices
\item{samples}{The samples matrix normalized}
\item{reference}{The reference matrix normalized}
}

\author{
  Cristiano Oliveira, Thomas Wolf
}

\examples{
data(sampleReadCounts)
data(referenceReadCounts)

normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts,
                                                 referenceReadCounts)
}
