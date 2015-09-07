\name{PlotBootstrapDistributions}
\alias{PlotBootstrapDistributions}
\title{
  PlotBootstrapDistributions
}

\description{
    Plots the generated bootstrap distribution as violin plots.
    Genes showing significant values
    are marked in a different color.
}

\usage{
PlotBootstrapDistributions(bootList,
                           reportTables,
                           outputFolder = getwd(),
                           sampleNames = NULL,
                           save = FALSE,
                           scale = 7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bootList}{
      List of bootstrapped read counts for each sample data
  }
  \item{reportTables}{
      List of report tables for each sample data
  }
  \item{outputFolder}{
    Path to the folder where the data plots will be created
  }
  \item{sampleNames}{
    List with sample names
  }
  \item{save}{
    Boolean to save the plots to the output folder
  }
  \item{scale}{
    Numeric scale factor
  }
}

\value{
A list with ggplot2 objects.
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

backgroundNoise <- Background(geneNames,
           samplesNormalizedReadCounts,
           referenceNormalizedReadCounts,
           bootList,
           replicates = replicates)

reportTables <- ReportTables(geneNames,
             samplesNormalizedReadCounts,
             referenceNormalizedReadCounts,
             bootList,
             backgroundNoise)

PlotBootstrapDistributions(bootList, reportTables, save = FALSE)
}
