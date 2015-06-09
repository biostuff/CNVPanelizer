\name{WriteReadCountsToXLSX}
\alias{WriteReadCountsToXLSX}

\title{
    WriteReadCountsToXLSX
}

\description{
    Writes both the read counts from reference and test samples to an xlsx file.
}

\usage{
    WriteReadCountsToXLSX(sampleReadCounts,
                        referenceReadCounts,
                        filepath = "readCounts.xls")
}

\arguments{
    \item{sampleReadCounts}{
        sampleReadCounts
    }
    \item{referenceReadCounts}{
        referenceReadCounts
    }
    \item{filepath}{
        filepath
    }
}

\value{
    Not returning any value
}

\author{
    Thomas Wolf, Cristiano Oliveira
}

\examples{
    \dontrun{
        WriteReadCountsToXLSX(sampleReadCounts,
                            referenceReadCounts,
                            filepath = "readCounts.xls")
    }
}
