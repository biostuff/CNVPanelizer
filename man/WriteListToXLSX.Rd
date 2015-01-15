\name{WriteListToXLSX}
\alias{WriteListToXLSX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  WriteListToXLSX
}
\description{
  Writes list of data frames to an xlsx file
}
\usage{
WriteListToXLSX(listOfDataFrames, filepath = "list.xlsx")
}
%- maybe also 'usage' for other objects documented here.
%%TODO create some
\arguments{
  \item{listOfDataFrames}{
    list of dataframes
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

  WriteListToXLSX(listOfDataFrames = exampleList, filepath = "list.xlsx")


}
}
