\name{StatusHeatmap}
\alias{StatusHeatmap}

\title{
    StatusHeatmap
}
\description{
    Generates a status heapmap for all samples analyzed
}
\usage{
            StatusHeatmap(dfData,
                          statusColors = c("Deletion" = "blue",
                                           "Normal" = "green",
                                           "Amplification" = "red"),
                          header = "Status Heatmap",
                          filepath = "CNVPanelizerHeatMap.png")
}

\arguments{
    \item{dfData}{
        data frame with the "Amplification", "Deletion" and "Normal" status
    }
    \item{statusColors}{
        A named vector with the colors associated with each level
    }
    \item{header}{
        Header text at the plot
    }
    \item{filepath}{
        Filepath where the generated heatmap is saved
    }
}

\value{
    Returns the filepath of the saved Heatmap
}

\author{
    Cristiano Oliveira
}

\examples{
    \dontshow{
      TRUE
    }
    \donttest{
            StatusHeatmap(dfData,
                          statusColors = c("Deletion" = "blue",
                                           "Normal" = "green",
                                           "Amplification" = "red"),
                          header = "Status Heatmap",
                          filepath = "CNVPanelizerHeatMap.png")
    }
}

