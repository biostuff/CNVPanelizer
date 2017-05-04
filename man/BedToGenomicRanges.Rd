\name{BedToGenomicRanges}
\alias{BedToGenomicRanges}
\title{
  BedToGenomicRanges
}
\description{
It generates a GenomicRanges object from a bed file.
Needs to be passed the correct number of the gene name column.
If the strings contain more information then just the gene name,
a splitting character (split) has to be defined. I.e GeneName1;Amplicon2
}
\usage{
BedToGenomicRanges(panelBedFilepath,
                   ampliconColumn,
                   split,
                   doReduce,
                   rangeExtend,
                   dropChromossomes,
                   skip)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{panelBedFilepath}{
      Filepath of the bed file.
}
  \item{ampliconColumn}{
    Number of the column that identifies the gene name in the bed file
    passed through \code{panelBedFilepath}.
}
  \item{split}{
     The character used as separator in the \code{ampliconColumn}.
     It is ";" by default.
}
  \item{doReduce}{
     Should overlapping ranges be merged.
}
  \item{rangeExtend}{
     Should the defined ranges be extended left and right
     by the given value.
     Affects the merging of overlapping regions and
     also read counting.
}

  \item{dropChromossomes}{
     Exclude vector of chromossomes (by the identifier used at the BED file)
}

  \item{skip}{
     How many lines should be skipped from the top of the bed file.
     The function assumes a bed file with column names.
     Thus default is skip = 1.
}

}

\value{
A \code{GenomicRanges} object containing information about the amplicons 
described in the bed file.
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
    bedFilepath <- file.path("someFile.bed")
    ampliconColumn <- 4
    genomicRangesFromBed <- BedToGenomicRanges(bedFilepath, ampliconColumn)}
}

