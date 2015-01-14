# CNVPanelizer


This is an R bioconductor package to use targeted sequencing data to reliably detect CNVs from clinical samples. To assess how reliable a change in reads counts in a specific region correlates with the presence of CNVs we implemented an algorithm which uses a subsampling strategy similar to Random Forest to predict the presence of reliable CNVs. We also introduce a novel method to correct for the background noise introduced by sequencing genes with a low number of amplicons. We describe the implementation of these models in the package <b>CNVPanelizer</b> and illustrate its usage to reliably detect CNVs on several simulation and real data examples including several code snippets when dealing with clinical data.


## Installation

You can install the stable version on
[CRAN](http://cran.rstudio.com/package=knitr):

```r
install.packages('CNVPanelizer', dependencies = TRUE)
```





## Usage

```r
library(CNVPanelizer)
?CNVPanelizer
```


## License

This package is free and open source software, licensed under GPL.