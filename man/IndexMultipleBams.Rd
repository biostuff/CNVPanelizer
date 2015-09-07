\name{IndexMultipleBams}
\alias{IndexMultipleBams}
\title{
  IndexMultipleBams
}
\description{
    Index a list of bam files if there is no index exists for the file entries in the list.
}
\usage{
IndexMultipleBams(bams, index_type = ".bam.bai")
}
\arguments{
  \item{bams}{
    A character vector of bam files to be indexed
  }
  \item{index_type}{
  The index file type extension
  }
}

\value{
Not returning any value
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
        files = c("file1.bam","file2.bam","file3.bam")
        IndexMultipleBams(bams = files)
    }
}
