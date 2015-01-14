\name{CNVPanelizer-package}
\alias{CNVPanelizer-package}
\alias{CNVPanelizer}
\docType{package}
\title{
Reliable CNV detection in targeted sequencing applications
}
\description{
This package implements uses a collection of non-matched normal tissue samples as reference set to detect CNV aberrations on a test data set.
}
\details{
\tabular{lc}{
Our approach uses a non-parametric bootstrap subsampling of the available reference samples, to estimate the distribution of read counts from targeted sequencing. As inspired by random forest, this is combined with a procedure that subsamples the amplicons associated with each of the targeted genes. The obtained information allows us to reliably classify the copy number aberrations on the gene level.
}
For a complete list of functions, use library(help = "CNVPanelizer").

\tabular{ll}{
Package: \tab CNVPanelizer\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2014-11-11\cr
License: \tab GPL\cr
}
}
\author{
Thomas Wolf <t.wolf@dkfz.de>\cr
Cristiano Oliveira <cristiano.oliveira@med.uni-heidelberg.de>\cr

}
\references{
TODO create some

@Book{abramowitz+stegun,
 author    = "Milton {Abramowitz} and Irene A. {Stegun}",
 title     = "Handbook of Mathematical Functions with
              Formulas, Graphs, and Mathematical Tables",
 publisher = "Dover",
 year      =  1964,
 address   = "New York",
 edition   = "ninth Dover printing, tenth GPO printing"
}




---
title: "CNVPanelizer"
output: html_document
bibliography: bibliography.bib
---


 @article{lin1973,
    author  = {Shen Lin and Brian W. Kernighan},
    title   = {An Effective Heuristic Algorithm for the Travelling-Salesman Problem},
    journal = {Operations Research},
    volume  = {21},
    year    = {1973},
    pages   = {498--516},
 }

}
\keyword{ R, Random Forest, CNV, Bootstraping, Panel Sequencing, Ion Torrent }
\seealso{
TODO create some
~~ Optional links to other man pages, e.g. ~~
}
