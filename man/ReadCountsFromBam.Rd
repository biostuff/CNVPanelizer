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
                minimumMappingQuality = 20,
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
    Genomic Range object as created by \code{BedToGenomicRanges}
  }
  \item{ampliconNames}{
    List of amplicon defining names
  }
  \item{minimumMappingQuality}{
     Specifies the minimum mapping quality
  }
  \item{removeDup}{
  Boolean value to remove duplicates. For reads with the same start site, end site
  and orientation only one is kept. For IonTorrent data this can be used to
  as an additional quality control. For Illumina data too many reads are being removed.
  }
}

\value{
    A matrix with read counts where the rows represents the Amplicons
    and the columns represents the samples.
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
        ReadCountsFromBam(bamFilenames,
                            sampleNames,
                            gr,
                            ampliconNames,
                            minimumMappingQuality = 20,
                            removeDup)
    }
}
