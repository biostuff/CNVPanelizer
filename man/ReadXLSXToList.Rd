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
ReadXLSXToList(filepath)
}

\arguments{
  \item{filepath}{
    filepath
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
