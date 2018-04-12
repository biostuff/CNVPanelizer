\name{RunCNVPanelizerShiny}
\alias{RunCNVPanelizerShiny}

\title{
    RunCNVPanelizerShiny
}
\description{
    Run CNVPanelizer as a shiny app
}
\usage{
    RunCNVPanelizerShiny(port = 8100)
}

\arguments{
    \item{port}{
        Port where the app will be listening
    }
}

\value{
    Not returning any value
}

\author{
    Thomas Wolf, Cristiano Oliveira
}

\examples{
    \dontshow{
      TRUE
    }
    \donttest{
       RunCNVPanelizerShiny(port=8080)
    }
}

