\name{SelectReferenceSetFromReadCounts}
\alias{SelectReferenceSetFromReadCounts}

\title{
    SelectReferenceSetFromReadCounts
}
\description{
    Select a reference set from read counts
}
\usage{
SelectReferenceSetFromReadCounts(allSamplesReadCounts,
                                 normalizationMethod = "tmm",
                                 referenceMaximumNumberOfElements = 30,
                                 referenceSelectionMethod = "kmeans",
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
    \item{referenceMaximumNumberOfElements}{
        Maximum number of elements to consider as reference (only to be used in case interquantile reference selection method)
    }
    \item{referenceSelectionMethod}{
        Reference selection method ("kmeans", ...)
    }
    \item{lowerBoundPercentage}{
    	Lower bound percentage (only to be used in case interquantile reference selection method)
    }
    \item{upperBoundPercentage}{
        Upper bound percentage (only to be used in case interquantile reference selection method)
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
SelectReferenceSetFromReadCounts(allSamplesReadCounts,
                                 normalizationMethod = "tmm",
                                 referenceMaximumNumberOfElements = 30,
                                 referenceSelectionMethod = "kmeans")
    }
}

