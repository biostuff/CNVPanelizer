\name{CNVPanelizer-package}
\alias{CNVPanelizer-package}
\alias{CNVPanelizer}
\docType{package}
\title{
Reliable CNV detection in targeted sequencing applications
}
\description{
This package implements an algorithm that uses a collection of non-matched
normal tissue samples as a reference set to detect CNV aberrations in data
generated from amplicon based targeted sequencing.
}
\details{
\tabular{lc}{
Our approach uses a non-parametric bootstrap subsampling of the available 
reference samples, to estimate the distribution of read counts from targeted
sequencing. As inspired by random forest, this is combined at each iteration
with a procedure that subsamples the amplicons associated with each of the targeted genes.
To estimate the background noise of sequencing genes with a low number of amplicons
a second subsampling step is performed. Both steps are combined to make a decision
on the CNV status. Thus classifying the copy number aberrations on the gene level.
}
For a complete list of functions, use library(help = "CNVPanelizer").

\tabular{ll}{
Package: \tab CNVPanelizer\cr
Type: \tab Package\cr
License: \tab GPL-3\cr
}
}
\author{
Thomas Wolf <thomas_wolf71f@gmx.de>\cr
Cristiano Oliveira <cristiano.oliveira@med.uni-heidelberg.de>\cr
}
\keyword{ R, Random Forest, CNV, Bootstraping, Panel Sequencing, Ion Torrent }
