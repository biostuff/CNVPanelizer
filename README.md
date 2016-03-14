# CNVPanelizer

This is an R bioconductor package to use targeted sequencing data to reliably
detect CNVs from clinical samples. To assess how reliable a change in reads
counts in a specific region correlates with the presence of CNVs we implemented
an algorithm which uses a subsampling strategy similar to Random Forest to
predict the presence of reliable CNVs. We also introduce a novel method to
correct for the background noise introduced by sequencing genes with a low
number of amplicons. We describe the implementation of these models in the
package <b>CNVPanelizer</b> and illustrate its usage to reliably detect CNVs
on several simulation and real data examples including several code snippets
when dealing with clinical data.

## Installation

You can install the stable version on
[Bioconductor](http://www.bioconductor.org/packages/release/bioc/html/CNVPanelizer.html)

```r
source("http://bioconductor.org/biocLite.R")
biocLite("CNVPanelizer")

```

If you want the latest version, install it directly from GitHub:

```r
library(devtools)
install_github("biostuff/CNVPanelizer")

```
## Motivation

Targeted sequencing, over the last few years, has become a mainstay in the
clinical use of next generation sequencing technologies. For the detection of
somatic and germline SNPs this has been proven to be a highly robust
methodology. One area of genomic analysis which is usually not covered by
targeted sequencing, is the detection of copy number variations (CNVs). While
a large number of available algorithms and software address the problem of CNV
detection in whole genome or whole exome sequencing, there are no such
established tools for targeted sequencing.

To assess how reliable a change in reads counts in a specific region correlates
with the presence of CNVs, we implemented an algorithm which uses a subsampling
strategy similar to Random Forest to predict the presence of reliable CNVs. We
also introduce a novel method to correct for the background noise introduced by
sequencing genes with a low number of amplicons.

## Usage

```r
library(CNVPanelizer)
?CNVPanelizer
```

## Contributors

cristiano.oliveira@med.uni-heidelberg.de

thomas_wolf71@gmx.de

## License

This package is free and open source software, licensed under GPL-3.
