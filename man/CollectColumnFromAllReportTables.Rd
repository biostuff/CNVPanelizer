\name{CollectColumnFromAllReportTables}
\alias{CollectColumnFromAllReportTables}

\title{
    CollectColumnFromAllReportTables
}
\description{
    Collect a single column from all report tables at the list
}
\usage{
	CollectColumnFromAllReportTables(reportTables, columnName)
}

\arguments{
    \item{reportTables}{
        A list of report tables
    }
    \item{columnName}{
	The column name
    }
}

\value{
    Returns a data frame with where the columns were collected from the entire list of report tables
}

\author{
    Cristiano Oliveira
}

\examples{
    \dontshow{
      TRUE
    }
    \donttest{
	CollectColumnFromAllReportTables(reportTables, columnName)
    }
}

