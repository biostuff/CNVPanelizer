\name{SelectReferenceSetByInterquartileRange}
\alias{SelectReferenceSetByInterquartileRange}

\title{
    SelectReferenceSetByInterquartileRange
}
\description{
    Select a reference set using a factor of the Interquartile Range
}
\usage{
    SelectReferenceSetByInterquartileRange(allSamplesReadCounts,
                                  	   normalizationMethod = "tmm",
                                           iqrFactor = 1)
}

\arguments{
    \item{allSamplesReadCounts}{
        All samples read counts matrix
    }
    \item{normalizationMethod}{
        tmm (trimmed mean of m values) or tss (total sum scaling)
    }
    \item{iqrFactor}{
        Interquantile range factor
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
                                  iqrFactor = 1)
    }
}

