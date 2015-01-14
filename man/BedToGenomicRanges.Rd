\name{BedToGenomicRanges}
\alias{BedToGenomicRanges}
\title{
  BedToGenomicRanges
}
\description{
It generates a GenomicRanges object from a bed file. Needs to be passed the correct column number of the amplicon names and the separator used to be parsed properly.
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
    Number of the column that identifies the amplicon name at the bed file passed through \code{panelBedFilepath}.
}
  \item{split}{
     Char used as separator. It is ";" by default.
}
}
\details{
TODO create some
}
\value{
A \code{GenomicRanges} object containing information about the amplicons described at the bed file.
}
\references{
TODO create some
}
\author{
Thomas Wolf, Cristiano Oliveira
}
\note{
TODO create some
}

\seealso{
 %%  \code{\link{GenomicRanges}}
 TODO create some
}
\examples{

bedFilepath <- file.path("someFile.bed")
ampliconColumn <- 4
# TODO how to make this code to not be executed work..?!?!!?
\dontrun{
BedToGenomicRanges(bedFilepath, ampliconColumn)
}
}
