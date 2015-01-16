\name{GenerateSynthetic}
\alias{GenerateSynthetic}
\title{
  GenerateSynthetic
}
\description{
    For testing and example purposes it is often good to have synthetically generated data.
    This function generates a synthetic read count matrix for such testing purposes.
    These examplem data sets are also being used in the function examples. The read count values
    are being drawn from poisson distribution.
}
\usage{
GenerateSynthetic(numberOfSamples,numberOfGenes,status = NULL,cnvProb = 1/10,label = "Sample",lambdaNor = 1000,lamdaAmp = 1.5 * lambdaNor,lambdaDel = 1/1.5 * lambdaNor,rdist = rpois,varyCoverage =TRUE,seed = 123)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numberOfSamples}{
       How many samples should be generated for the data set.
  }
  
  \item{numberOfGenes}{
       How many genes should the data set have.
       1. Gene has one amplicon
       2. Gene has two amplicons
       and so on ....
  }
  
  \item{status}{
        A vector of double values
  }
  
  \item{cnvProb}{
       The probability that a gene shows a copy number aberration.
  }
  
   \item{label}{
       The base string for sample naming.
  }
  
  \item{lambdaNor}{
       The standard read count if no abberation is present
  }
  
  
  \item{lambdaAmp}{
       The lambda read count if an amplification is present
  }
  
  \item{lambdaDel}{
       A lambda read count if a deletion is present.
  }
  
  \item{rdist}{
       The function used to generate a read count distribution.
       Standard is rdist, but other one parameter distribution functions
       can be passed to the function.
  }
  
  \item{varyCoverage}{
       Should all samples have the same average read count or should this vary.
  }
  
   \item{seed}{
       To make results reproducible, which seed should be used.
  }
}
\value{

A synthetically generated read count matrix, with genes(amplicons) in the rows
and the samples in the columns

}

\author{
  Thomas Wolf, Cristiano Oliveira
}



\examples{

#generating a reference set
#distribution setting
lambdaNor = 1000


#number of genes
numberOfGenes = 10

#number of reference samples
numberOfReferenceSamples = 100

referenceSamples = GenerateSynthetic(numberOfReferenceSamples, numberOfGenes, cnvProb = 1/100, label = "Reference")


#generating a sample set

#number of tumor samples
numberOfTestSamples = 2


testSamples =GenerateSynthetic(numberOfTestSamples,numberOfGenes, cnvProb = 1/5)

#generating a customized sample set

#number of tumor samples
numberOfTestSamples = 2

status =list(c(1,2,0.5,0.75,1),c(1,2,0.5,0.75,1))

testSamples =GenerateSynthetic(numberOfTestSamples,numberOfGenes, cnvProb = 1/5)



}

