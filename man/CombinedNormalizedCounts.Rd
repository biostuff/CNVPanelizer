\name{CombinedNormalizedCounts}
\alias{CombinedNormalizedCounts}
\title{
  CombinedNormalizedCounts
}
\description{
  It combines the rows of the matrix of reference read counts and the matrix of sample reads and normalized them together.
}
\usage{
CombinedNormalizedCounts(sampleCounts, referenceCounts, gr, amplicons)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleCounts}{
      Matrix with sample read counts
  }
    \item{referenceCounts}{
      Matrix with reference read counts
  }
    \item{gr}{
      A genomic ranges object
  }
    \item{amplicons}{
      A vector with the amplicon names used at the reference and sample matrices
  }
}
\details{
TODO create some 
  If necessary, more details than the description above
}
\value{
A normalized matrix composed of data from the reference and sample matrices
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
%% TODO create some
# normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts,
#                                                     referenceReadCounts,
#                                                     genomicRangesFromBed,
#                                                     amplicons = ampliconNames)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
