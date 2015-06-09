\name{ReadCountsFromBam}
\alias{ReadCountsFromBam}
\title{
  ReadCountsFromBam
}
\description{
  Returns a matrix with the read counts from a set of bam files.
}
\usage{
ReadCountsFromBam(bamFilenames,
                sampleNames,
                gr,
                ampliconNames,
                removeDup = FALSE)
}
\arguments{
  \item{bamFilenames}{
    Vector of bamfile filepaths
  }
  \item{sampleNames}{
    Vector of sample names to be used as colums names instead of bam filepaths
  }
  \item{gr}{
    Genomic Range object
  }
  \item{ampliconNames}{
    List of amplicon names
  }
  \item{removeDup}{
  Boolean value in to remove duplicates
  }
}

\value{
    A matrix with read counts where the rows represents the Amplicons
    and the columns represents the samples.
}

\examples{

\dontshow{
TRUE
}

\dontrun{
    ReadCountsFromBam(bamFilenames, sampleNames, gr, ampliconNames, removeDup)
}
}
