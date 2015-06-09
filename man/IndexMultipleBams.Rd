\name{IndexMultipleBams}
\alias{IndexMultipleBams}
\title{
  IndexMultipleBams
}
\description{
    Index a list of bam files if there is no index for the entries in the list.
}
\usage{
IndexMultipleBams(bams, multicore = FALSE, ncores = 2, index_type = ".bam.bai")
}
\arguments{
  \item{bams}{
    A character vector of bam files to be indexed
  }
  \item{multicore}{
  Use multicore
  }
  \item{ncores}{
  The number of cores, if multicore is set to TRUE
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
  \dontrun{
    files = c("file1.bam","file2.bam","file3.bam")
    IndexMultipleBams(bams = files)
  }
}
