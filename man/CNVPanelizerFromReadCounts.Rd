\name{CNVPanelizerFromReadCounts}
\alias{CNVPanelizerFromReadCounts}

\title{
    CNVPanelizerFromReadCounts
}
\description{
    Performs the workflow analysis with CNVPanelizer from the read counts and splitting the batch of samples analyzed
}
\usage{
            CNVPanelizerFromReadCounts(sampleReadCounts,
                                       referenceReadCounts,
                                       genomicRangesFromBed,
                                       numberOfBootstrapReplicates = 10000,
                                       normalizationMethod = "tmm",
                                       robust = TRUE,
                                       backgroundSignificanceLevel = 0.05,
                                       outputDir = file.path(getwd(), "CNVPanelizer"))
}

\arguments{
    \item{sampleReadCounts}{
        samples read counts matrix
    }
    \item{referenceReadCounts}{
        reference read counts matrix
    }
    \item{genomicRangesFromBed}{
        genomic ranges from bed
    }
    \item{numberOfBootstrapReplicates}{
        number of bootstrap replicates
    }
    \item{normalizationMethod}{
      Normalization method ("tmm" or "tss") 
    }
    \item{robust}{
       if TRUE, the median is used instead of mean
    }
    \item{backgroundSignificanceLevel}{
       The background Significance Level
    }
    \item{outputDir}{
       Output directory
    }
}

\value{
    Returns a list with the results of each samples analyzed
}

\author{
    Cristiano Oliveira
}

\examples{
    \dontshow{
      TRUE
    }
    \donttest{
            CNVPanelizerFromReadCounts(sampleReadCounts,
                                       referenceReadCounts,
                                       genomicRangesFromBed,
                                       numberOfBootstrapReplicates = 10000,
                                       normalizationMethod = "tmm",
                                       robust = TRUE,
                                       backgroundSignificanceLevel = 0.05,
                                       outputDir = file.path(getwd(), "CNVPanelizer"))
    }
}

