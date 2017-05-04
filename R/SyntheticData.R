#distribution setting
# nocov start
lambdaNor = 1000
lambdaAmp = 1.5 * lambdaNor 
lambdaDel = 1/1.5 * lambdaNor

#number of genes
numberOfGenes = 5

#number of tumor samples
numberOfTestSamples = 2

#number of reference samples
numberOfReferenceSamples = 100

GenerateSynthetic <- function(numberOfSamples,
                              numberOfGenes,
                              status = NULL,
                              cnvProb = 1/10,
                              label = "Sample",
                              lambdaNor = 1000,
                              lamdaAmp = 1.5 * lambdaNor,
                              lambdaDel = 1/1.5 * lambdaNor,
                              rdist = rpois,
                              seed = 123) {

    set.seed(seed)
  
    # To avoid complaints from "R CMD check" ..
    i <- NULL
    j <- NULL
    samples = foreach(i = 1:numberOfSamples,.combine=cbind) %:% 
        foreach(j = 1:numberOfGenes,.combine =  c) %do% {
            if(is.null(status)) {
                lambda = sample(c(lambdaDel, lambdaNor, lambdaAmp),
                                1,
                                prob = c(cnvProb,1,cnvProb))
            } else {
                lambda =  status[[i]][j] * lambdaNor
            }
            rdist(j,lambda)
        }

    colnames(samples) = paste0(label,"_",1:numberOfSamples)
    #name of genes
    geneNames = paste0("Gene_",rep(1:numberOfGenes,1:numberOfGenes))
    rownames(samples) = geneNames
    return(samples)
}

status = list(c(1,2,0.5,0.75,1), c(1,2,0.5,0.75,1))
testSamples = GenerateSynthetic(numberOfTestSamples,
                                numberOfGenes,
                                cnvProb = 1/10,
                                status = status)

referenceSamples = GenerateSynthetic(numberOfReferenceSamples,
                                     numberOfGenes,
                                     cnvProb = 1/100,
                                     label = "Reference")
# nocov end

