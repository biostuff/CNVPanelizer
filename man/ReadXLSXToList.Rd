\name{ReadXLSXToList}
\alias{ReadXLSXToList}
\title{
  ReadXLSXToList
}
\description{
    Reads a list of read count matrices from a xlsx as generated
    by \code{WriteReadCountsToXLSX}
}
\usage{
ReadXLSXToList(filepath, rowNames = TRUE, colNames = TRUE)
}

\arguments{
  \item{filepath}{
    filepath
  }
  \item{rowNames}{
    if row names should be included
  }
  \item{colNames}{
    if col names should be included
  }
}

\value{
    A list of read count matrices
}

\author{
    Thomas Wolf, Cristiano Oliveira
}

\examples{
    % To avoid biocCheck complaints..
    \dontshow{
        TRUE
    }
    \donttest{
        ReadXLSXToList(filepath)
    }
}
