\name{SelectReferenceSetByKmeans}
\alias{SelectReferenceSetByKmeans}

\title{
    SelectReferenceSetByKmeans
}
\description{
    Select a reference set using Kmeans
}
\usage{
    SelectReferenceSetByKmeans(allSamplesReadCounts,
			       normalizationMethod = "tmm",
			       referenceNumberOfElements)
}

\arguments{
    \item{allSamplesReadCounts}{
        All samples read counts matrix
    }
    \item{normalizationMethod}{
        tmm (trimmed mean of m values) or tss (total sum scaling)
    }
    \item{referenceNumberOfElements}{
        Number of elements to select for the reference set
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
    SelectReferenceSetByKmeans(allSamplesReadCounts, 
                               normalizationMethod = "tmm", 
                               referenceNumberOfElements)
    }
}

