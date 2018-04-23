\name{CNVPanelizerFromReadCountsHELPER}
\alias{CNVPanelizerFromReadCountsHELPER}

\title{
    CNVPanelizerFromReadCountsHELPER
}
\description{
    Helper to performs the workflow analysis with CNVPanelizer from the read counts and splitting the batch of samples analyzed
}
\usage{
            CNVPanelizerFromReadCountsHELPER(sampleReadCounts,
                                             referenceReadCounts,
                                             genomicRangesFromBed,
                                             numberOfBootstrapReplicates = 10000,
                                             normalizationMethod = "tmm",
                                             robust = TRUE,
                                             backgroundSignificanceLevel = 0.05,
                                             outputDir = file.path(getwd(), "CNVPanelizer"),
                                             splitSize = 5)
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
    \item{splitSize}{
       Split size of the batches analyzed
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
            CNVPanelizerFromReadCountsHELPER(sampleReadCounts,
                                             referenceReadCounts,
                                             genomicRangesFromBed,
                                             numberOfBootstrapReplicates = 10000,
                                             normalizationMethod = "tmm",
                                             robust = TRUE,
                                             backgroundSignificanceLevel = 0.05,
                                             outputDir = file.path(getwd(), "CNVPanelizer"),
                                             splitSize = 5)
    }
}

