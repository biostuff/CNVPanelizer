\name{IndexMultipleBams}
\alias{IndexMultipleBams}
\title{
  IndexMultipleBams
}
\description{
# TODO \cr
  Index a list of bam files if there is no index yet for the entries in the list
}
\usage{
IndexMultipleBams(bams, multicore = FALSE, ncores = 2, index_type = ".bam.bai")
}
\arguments{
  \item{bams}{
    List of bam files to be indexed
  }
  \item{multicore}{
  Use multicore
  }
  \item{ncores}{
  The Number of cores, case you multicore is set to TRUE
  }
  \item{index_type}{
  The index file type extension
  }
}
\details{
TODO create some
  If necessary, more details than the description above
}
\value{
TODO create some
  Describe the value returned
  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
 ...
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
#IndexMultipleBams(bams,
#                   multicore = FALSE,
#                   ncores = 2,
#                   index_type = ".bam.bai")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
