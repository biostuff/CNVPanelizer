\name{SelectReferenceSetByPercentil}
\alias{SelectReferenceSetByPercentil}

\title{
    SelectReferenceSetByPercentil
}
\description{
    Select a reference set using percentiles
}
\usage{
    SelectReferenceSetByPercentil(allSamplesReadCounts,
                                  normalizationMethod = "tmm",
                                  lowerBoundPercentage = 1,
                                  upperBoundPercentage = 99)
}

\arguments{
    \item{allSamplesReadCounts}{
        All samples read counts matrix
    }
    \item{normalizationMethod}{
        tmm (trimmed mean of m values) or tss (total sum scaling)
    }
    \item{lowerBoundPercentage}{
        Lower bound percentage
    }
    \item{upperBoundPercentage}{
        Upper bound percentage
    }
}

\value{
    Returns a list of sample identifiers to be used as reference
}

\author{
    Cristiano Oliveira
}

\examples{
    \dontshow{
      TRUE
    }
    \donttest{
    SelectReferenceSetByPercentil(allSamplesReadCounts,
                                  normalizationMethod = "tmm",
                                  lowerBoundPercentage = 1,
                                  upperBoundPercentage = 99)
    }
}

