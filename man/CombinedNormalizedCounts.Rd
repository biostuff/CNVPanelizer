\name{CombinedNormalizedCounts}
\alias{CombinedNormalizedCounts}
\title{
  CombinedNormalizedCounts
}
\description{
    This function makes use of cn.mops::normalizeChromosomes to normalize
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
         Matrix  with sample read counts (rows: amplicons, columns: samples)
    }
    \item{referenceCounts}{
        Matrix with reference read counts (rows: amplicons, columns: samples)
    }
    \item{ampliconNames}{
        A vector with the amplicon names for the ref. and sample matrices
    }
}

\value{
A list object with matrix
\item{samples}{The samples matrix normalized}
\item{reference}{The reference matrix normalized}
}

\author{
  Thomas Wolf, Cristiano Oliveira
}

\examples{
data(sampleReadCounts)
data(referenceReadCounts)

normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts,
                                                 referenceReadCounts)
}
