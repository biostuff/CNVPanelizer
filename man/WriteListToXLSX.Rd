\name{WriteListToXLSX}
\alias{WriteListToXLSX}

\title{
    WriteListToXLSX
}
\description{
    Writes list of data frames to an xlsx file
}
\usage{
    WriteListToXLSX(listOfDataFrames, filepath = "list.xlsx")
}

\arguments{
    \item{listOfDataFrames}{
        list of dataframes
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
    \dontshow{
      TRUE
    }
    \donttest{
        WriteListToXLSX(listOfDataFrames = exampleList, filepath = "list.xlsx")
    }
}
