\name{ReadCountsFromBam}
\alias{ReadCountsFromBam}
\title{
  ReadCountsFromBam
}
\description{
  Returns a matrix with the read counts from a list of bam files.
}
\usage{
ReadCountsFromBam(bamFilenames, sampleNames, gr, ampliconNames, removeDup = FALSE)
}
\arguments{
  \item{bamFilenames}{
    List of bamfiles filepaths
}
  \item{sampleNames}{
    List of sample names to be used instead of bam filepaths
}
  \item{gr}{
    Genomic Range
}
  \item{ampliconNames}{
    List of amplicon names
}
  \item{removeDup}{
  Boolean value in to remove duplicates
}
}
\details{
TODO create some
}
\value{
A matrix with read counts where the rows represents the Amplicons and the columns represents the samples.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 TODO create some
 objects to See Also as \code{\link{help}},
}
\examples{
# ReadCountsFromBam(bamFilenames, sampleNames, gr, ampliconNames, removeDup = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
