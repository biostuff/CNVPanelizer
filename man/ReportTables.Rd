\name{ReportTables}
\alias{ReportTables}
\title{
  ReportTables
}
\description{
    This function generates the final report of the CNV detection procedure.
    One data frame is generated for each sample of interest.
}
\usage{

ReportTables(geneNames,
             samplesNormalizedReadCounts,
             referenceNormalizedReadCounts,
             bootList,
             backgroundNoise)
}
\arguments{
    \item{geneNames}{
        Describe \code{geneNames} here
    }
    \item{samplesNormalizedReadCounts}{
        Describe \code{samplesNormalizedReadCounts} here
    }
    \item{referenceNormalizedReadCounts}{
        Describe \code{referenceNormalizedReadCounts} here
    }
    \item{bootList}{
        A list as returned by the  \code{BootList} function
    }
    \item{backgroundNoise}{
        A list of background noise as returned by the
        \code{Background} function
    }
}

\value{
   Returns a list of tables, one for each sample of interest.
   Each of these tables contains numerical information
   of the aberration status of each gene. For a detailed description
   see the Vignette.
}

\author{
  Thomas Wolf, Cristiano Oliveira
}

\examples{

data(sampleReadCounts)
data(referenceReadCounts)
## Gene names should be same size as row columns
geneNames <- row.names(referenceReadCounts)

ampliconNames <- NULL

normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts,
                                                 referenceReadCounts,
                                                 ampliconNames = ampliconNames)

# After normalization data sets need to be splitted again to perform bootstrap
samplesNormalizedReadCounts = normalizedReadCounts["samples"][[1]]
referenceNormalizedReadCounts = normalizedReadCounts["reference"][[1]]

# Should be used values above 10000
replicates <- 10

# Perform the bootstrap based analysis
bootList <- BootList(geneNames,
                     samplesNormalizedReadCounts,
                     referenceNormalizedReadCounts,
                     replicates = replicates)

backgroundNoise = Background(geneNames,
                             samplesNormalizedReadCounts,
                             referenceNormalizedReadCounts,
                             bootList,
                             replicates = replicates)

reportTables <- ReportTables(geneNames,
             samplesNormalizedReadCounts,
             referenceNormalizedReadCounts,
             bootList,
             backgroundNoise)
}
