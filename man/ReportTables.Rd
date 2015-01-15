\name{ReportTables}
\alias{ReportTables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ReportTables
}
\description{
  This function generates the final report of the copy number detection procedure.
  One data frame is generated for each sample of interest.
}
\usage{
ReportTables(bootList,
             geneNames,
             backgroundNoise,
             referenceNormalizedReadCounts,
             samplesNormalizedReadCounts,
             tableNames)
}
%% TODO create some argument list
\arguments{
  \item{bootList}{
     A list as returned by the  \code{BootList} function
}
  \item{geneNames}{
     Describe \code{geneNames} here
}
  \item{backgroundNoise}{
     A list of background noise as  returned by  the \code{Background} function.
}
  \item{referenceNormalizedReadCounts}{
     Describe \code{referenceNormalizedReadCounts} here
}
  \item{samplesNormalizedReadCounts}{
     Describe \code{samplesNormalizedReadCounts} here
}
  \item{tableNames}{
     Describe \code{tableNames} here
}
}

\value{
   Returns a list of tables, one for each sample of interest.
   Each of these tables contains numerical information
   of the aberration status of each gene.
}

\author{
  Thomas Wolf, Cristiano Oliveira
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

\dontrun{
 ReportTables(bootList,
             geneNames,
             backgroundNoise,
             referenceNormalizedReadCounts,
             samplesNormalizedReadCounts,
             tableNames)
}

}
