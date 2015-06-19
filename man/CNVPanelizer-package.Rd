\name{CNVPanelizer-package}
\alias{CNVPanelizer-package}
\alias{CNVPanelizer}
\docType{package}
\title{
Reliable CNV detection in targeted sequencing applications
}
\description{
This package implements uses a collection of non-matched normal tissue samples
as reference set to detect CNV aberrations on a test data set.
}
\details{
\tabular{lc}{
Our approach uses a non-parametric bootstrap subsampling of the available 
reference samples, to estimate the distribution of read counts from targeted
sequencing. As inspired by random forest, this is combined with a procedure
that subsamples the amplicons associated with each of the targeted genes. 
The obtained information allows us to reliably classify the copy number
aberrations on the gene level.
}
For a complete list of functions, use library(help = "CNVPanelizer").

\tabular{ll}{
Package: \tab CNVPanelizer\cr
Type: \tab Package\cr
License: \tab GPL-3\cr
}
}
\author{
Thomas Wolf <t.wolf@dkfz.de>\cr
Cristiano Oliveira <cristiano.oliveira@med.uni-heidelberg.de>\cr
}
\keyword{ R, Random Forest, CNV, Bootstraping, Panel Sequencing, Ion Torrent }
