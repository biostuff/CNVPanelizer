\name{NormalizeCounts}
\alias{NormalizeCounts}
\title{
  NormalizeCounts
}
\description{
    This function normalize counts use of Total sum scaling or NOISeq::tmm to normalize
    the read counts
}
\usage{
    NormalizeCounts(allCounts,
                    method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{allCounts}{
         Matrix or vector with sample read counts (rows: amplicons, columns: samples)
    }
    \item{method}{
        either "tmm" (trimmed mean of m values) or "tss"(total sum scaling)
    }
}

\value{
A matrice
\item{samples}{The samples matrix normalized}
}

\author{
  Cristiano Oliveira, Thomas Wolf
}

\examples{
data(sampleReadCounts)

normalizedReadCounts <- NormalizeCounts(sampleReadCounts)
}
