\name{WriteReadCountsToXLSX}
\alias{WriteReadCountsToXLSX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  WriteReadCountsToXLSX
}
\description{
  Writes both the read counts from reference and test samples to an xlsx file.
}
\usage{
WriteReadCountsToXLSX(sampleReadCounts, referenceReadCounts, filepath = "readCounts.xls")
}
%- maybe also 'usage' for other objects documented here.
%%TODO create some
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


\author{
  Thomas Wolf, Cristiano Oliveira
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
  WriteReadCountsToXLSX(sampleReadCounts, referenceReadCounts, filepath = "readCounts.xls")


}

}
