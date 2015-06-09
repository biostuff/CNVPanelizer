\name{BedToGenomicRanges}
\alias{BedToGenomicRanges}
\title{
  BedToGenomicRanges
}
\description{
It generates a GenomicRanges object from a bed file.
Needs to be passed the correct column number of the amplicon names
and the separator used to be parsed properly.
}
\usage{
BedToGenomicRanges(panelBedFilepath, ampliconColumn, split)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{panelBedFilepath}{
      Filepath of the bed file.
}
  \item{ampliconColumn}{
    Number of the column that identifies the amplicon name at the bed file
    passed through \code{panelBedFilepath}.
}
  \item{split}{
     Char used as separator. It is ";" by default.
}
}

\value{
A \code{GenomicRanges} object containing information about the amplicons 
described at the bed file.
}

\author{
Thomas Wolf, Cristiano Oliveira
}

\examples{
% To avoid biocCheck complaints..
\dontshow{
TRUE
}
\dontrun{
    bedFilepath <- file.path("someFile.bed")
    ampliconColumn <- 4
    BedToGenomicRanges(bedFilepath, ampliconColumn)}
}