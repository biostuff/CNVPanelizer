\name{WriteListToXLSX}
\alias{WriteListToXLSX}

\title{
    WriteListToXLSX
}
\description{
    Writes list of data frames to an xlsx file
}
\usage{
    WriteListToXLSX(listOfDataFrames,
                    multipleFiles = FALSE,
                    outputFolder = file.path(getwd(), "xlsx"),
                    filepath = "list.xlsx")
}

\arguments{
    \item{listOfDataFrames}{
        list of dataframes
    }
    \item{multipleFiles}{
        If should be generated on single file with all results or multiple files
    }
    \item{outputFolder}{
        Output folder
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
