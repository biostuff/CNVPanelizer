#library(stringr)
#library(snow)
#library(data.table)
#library(WriteXLS)
library(plyr)
library(cn.mops)
library(foreach)
library(matrixStats)
library(Rsamtools)
library(exomeCopy)
library(xlsx)
library(ggplot2)


# TODO generate functions to generate fake data 
# lambda <- 1
# rpois(1000, 1000*lambda)
# hist(rpois(1000, 1000*lambda))

# rJava causing errors..
# Sys.setenv(JAVA_HOME='C:\\programs\\Java\\jre1.8.0_20')

# Have added the directory where jvm.dll was..
# C:\programs\Java\jre1.8.0_20\bin\server

# Find the lib folder (where the packages are installed)
# .libPaths()

# set same proxy as the browser
# Method1: Invoking R using --internet2

# THIS DOES NOT WORK!?!
#setInternet2(TRUE)

# THIS WORKS!!!!!!
#Setinternet2=TRUE

# If you want internet2 to be used everytime you use R you could add the following line to the Rprofile.site file which is located in R.x.x\etc\Rprofile.site
# utils::setInternet2(TRUE)



##################################################################################
# functions to read the read counts
##################################################################################


# 
# SampleNameFromBam <-function(bamFileNames) {
#     sampleNames = foreach(i = seq_along(bamFileNames),.combine = c) %do% {
#     header <- scanBamHeader(bamFileNames[i])
#     #print(header)
# missing the dollar symbol before text because of groovy..
#     #sm <- header[[1]] text[2][[1]][5]
#     sm <- header[[1]][2] text["@RG"][[1]][10]
#     strsplit(sm,split = ":")[[1]][2]
#   }
#   return(sampleNames)
# }
# 
# bamFileNames = "Z:\\Projekt Results\\Copy Numer Analysis\\ReferenceBam\\LCPv1\\E-2014-01020_1PG-94_IonXpress_019_rawlib.bam"
# SampleNameFromBam(bamFileNames)
#
# BedToGenomicRanges <- function(bedFolder, panel, ampliconColumn, split = ";") {
BedToGenomicRanges <- function(panelBedFilepath, ampliconColumn, split = ";") {
    #load the bed file
#     segments <- read.table(paste0(bedFolder, panel, ".bed"), sep = "\t", as.is = TRUE, skip = 1)
    segments <- read.table(panelBedFilepath, sep = "\t", as.is = TRUE, skip = 1)
    gr <- GRanges(segments[, 1], IRanges(segments[, 2], segments[, 3]))
    #get the amplicon names form the segments
    amplicons = segments[, ampliconColumn]
    #get the genes from segments
    splitted = strsplit(amplicons, split = split)
    # TODO Required because of package checking complaints
    i <- NULL
    genes = foreach(i = seq_along(splitted)) %do% {
      splitted[[i]][1]
    }
    
    
    genes = unlist(genes)
    # missing the dollar symbol before text because of groovy..
#    elementMetadata(gr) geneNames = genes
#    elementMetadata(gr) ampliconNames = paste(seqnames(gr), start(gr), end(gr), sep = "_", amplicons)

    elementMetadata(gr)["geneNames"] = genes
    elementMetadata(gr)["ampliconNames"] = paste(seqnames(gr), start(gr), end(gr), sep = "_", amplicons)

    return(gr)
}

#load the files you want to analyze
ReadCountsFromBam <- function(bamFilenames, sampleNames, gr, ampliconNames = elementMetadata(gr)$ampliconNames, removeDup = FALSE) {
#ReadCountsFromBam <- function(bamFilenames, sampleNames, gr, ampliconNames = elementMetadata(gr)["ampliconNames"], removeDup = FALSE) {
    # TODO Because of package check complaints
    i <- NULL
    curbam = foreach(i = seq_along(bamFilenames), .combine = cbind) %do% {
        countBamInGRanges(bamFilenames[i], gr, remove.dup = removeDup, min.mapq = 20, get.width = TRUE)
    }
    colnames(curbam) = sampleNames
    rownames(curbam) = ampliconNames
    return(curbam)
}




# ReadCountsFromBam <- function(bamFilenames, sampleNames, gr, ampliconNames, removeDup = FALSE) {
#   # TODO Because of package check complaints
#   i <- NULL
#   curbam = foreach(i = seq_along(bamFilenames), .combine = cbind) %do% {
#       countBamInGRanges(bamFilenames[i], gr, remove.dup = removeDup, min.mapq = 20, get.width = TRUE)
#   }
#   colnames(curbam) <- sampleNames
#   rownames(curbam) <- ampliconNames
#   return(curbam)
# }

GenerateSynthetic <- function(numberOfSamples,
                              numberOfGenes,
                              amplCount,
                              status = NULL,
                              cnvProb = 1/10,
                              label = "Sample",
                              lambdaNor = 1000,
                              lambdaAmp = 1.5 * lambdaNor,
                              lambdaDel = 1/1.5 * lambdaNor,
                              rdist = rpois,
                              varyCoverage =TRUE,
                              seed = 123) {
  set.seed(seed)
  if(varyCoverage) {
    coverageVariability = sample(1:10,numberOfSamples,replace = TRUE)
  } else{
    coverageVariability = rep(1,numberOfSamples)
  }

  #how many amplicons for each sample
  #minAmpl = 1
  #maxAmpl = 10
  #seqAmpl = minAmpl:maxAmpl
  #amplCount = sample(seqAmpl,numberOfGenes,replace = TRUE,prob = 1/seqAmpl)

    # TO avoid complaints from "R CMD check" ..
    i <- NULL
    j <- NULL
  
    samples = foreach(i = 1:numberOfSamples, .combine=cbind) %:%
    foreach(j = 1:numberOfGenes, .combine = c) %do% {
      if(is.null(status)) {
        lambda = sample(c(coverageVariability[i] * lambdaDel, coverageVariability[i] * lambdaNor, coverageVariability[i] * lambdaAmp), 1, prob = c(cnvProb, 1, cnvProb))
      } else {
        lambda =  status[[i]][j] * lambdaNor * coverageVariability[i]
      }
    rdist(amplCount[j],lambda)
    }

  colnames(samples) = paste0(label, "_", 1:numberOfSamples)
  geneNames = paste0("Gene_", rep(1:numberOfGenes, amplCount))
  rownames(samples) = geneNames
  return(samples)
}

#get the combined counts
#CombinedNormalizedCounts <- function(sampleCounts, referenceCounts, gr, amplicons = elementMetadata(gr) ampliconNames) {
# TODO change amplicons to ampliconNames to have same parameter with the same name as in the other functions..
CombinedNormalizedCounts <- function(sampleCounts, referenceCounts, amplicons = NULL) {
  #combine call samples and reference
  allCounts = cbind(sampleCounts, referenceCounts)
  classes = rep(c("samples", "reference"), c(ncol(sampleCounts), ncol(referenceCounts)))
  #normalize all samples
  bamDataRangesNorm = as.matrix((normalizeGenome(allCounts, normType = "median")))
  if(!is.null(amplicons)) {
    #add the amplicon names as rownames
    rownames(bamDataRangesNorm) = amplicons
  }
  #get the count information for one sample
  #add small value to prevent a zero read count
  normalizedSamples = bamDataRangesNorm[, which(classes == "samples")] + 0.00001
  normalizedReference = bamDataRangesNorm[, which(classes == "reference")] +  0.00001
  return(list(samples = normalizedSamples, reference = normalizedReference))
}





#get the combined counts
#CombinedNormalizedCounts <- function(sampleCounts, referenceCounts, gr, amplicons = elementMetadata(gr) ampliconNames) {
# TODO change amplicons to ampliconNames to have same parameter with the same name as in the other functions..
#CombinedNormalizedCounts <- function(sampleCounts, referenceCounts, gr, amplicons = elementMetadata(gr)["ampliconNames"]) {
    #combine call samples and reference
#    allCounts = cbind(sampleCounts, referenceCounts)
#    classes = rep(c("samples", "reference"), c(ncol(sampleCounts), ncol(referenceCounts)))
    #normalize all samples
#    bamDataRangesNorm = as.matrix((normalizeGenome(allCounts, normType = "median")))
#    if(!is.null(amplicons)) {
#        #add the amplicon names as rownames
#        rownames(bamDataRangesNorm) = amplicons
#    }
    #get the count information for one sample
    #add one tp revent a zero read count
#    normalizedSamples = bamDataRangesNorm[, which(classes == "samples")] + 0.00001
#    normalizedReference = bamDataRangesNorm[, which(classes == "reference")] +  0.00001
#    return(list(samples = normalizedSamples, reference = normalizedReference))
#}

# #save the normalized and unnormalized tables as an excel file
# saveReadCountTable <- function(countTable,file) {
#   
#   
#   write.csv(countTable,file = file)
#   
#   
# }




##################################################################################
# helper functions
##################################################################################

#get the positions for each gene as a list
IndexGenesPositions = function(genes) {
    positions = seq_along(genes)
    genesPos = split(positions, genes)
    return(genesPos)
}

#how many numbers in a vector are above a given cutoff
PercentAboveValue <- function(vector, value) {
    sum(vector > value)/length(vector)
}

#how many numbers in a vector are below a given cutoff
PercentBelowValue <- function(vector,value) {
    sum(vector < value)/length(vector)
}

#index the bam files if there is no index yet
IndexMultipleBams <- function(bams, multicore = FALSE, ncores = 2, index_type = ".bam.bai") {
    #check if the index already exists and need to be indexed
    potentialBaiFilenames <- gsub(".bam", bams, replacement = index_type, ignore.case = TRUE)
    bamsToBeIndexed <- bams[!sapply(potentialBaiFilenames, file.exists)]
    if(length(bamsToBeIndexed) > 0) {
        #index the bams
        if(multicore) {
            mclapply(bamsToBeIndexed, indexBam, mc.cores = ncores)
        } else {
            lapply(bamsToBeIndexed, indexBam)
        }
    }
}

# # extension of filepath determines if the format will be excel 2003 (.xls) or excel 2007 (.xlsx)
# WriteReadCountsToXLS <- function(sampleReadCounts, referenceReadCounts, filepath = "readCounts.xls") {
#   # writeCounts <- function(sampleReadCounts,referenceReadCounts,file = "") {
#   
#   writeSampleReadCounts <- as.data.frame(sampleReadCounts)
#   #rownames(writeSampleReadCounts) <- paste0(seq_along(ampliconNames), "_",ampliconNames)
#   
#   writeReferenceReadCounts <- as.data.frame(referenceReadCounts)
#   
#   writeList <- list(SampleReadCounts = writeSampleReadCounts, referenceReadCounts = writeReferenceReadCounts)
# 
# #  WriteXLS("writeList", paste0(outputFolder, file), AdjWidth = TRUE, BoldHeaderRow = TRUE, row.names = TRUE, col.names = TRUE)
# 
# #  WriteXLS("writeList", filepath, AdjWidth = TRUE, BoldHeaderRow = TRUE, row.names = TRUE, col.names = TRUE)
# WriteXLS("writeReferenceReadCounts", filepath, AdjWidth = TRUE, col.names = TRUE)
# #WriteXLS("writeList", filepath, AdjWidth = TRUE, BoldHeaderRow = TRUE, row.names = TRUE, col.names = TRUE)
# }

WriteReadCountsToXLSX <- function(sampleReadCounts, referenceReadCounts, filepath = "readCounts.xls") {
    write.xlsx(sampleReadCounts,
               filepath,
               sheetName = "sampleReadCounts")
    write.xlsx(referenceReadCounts,
               filepath,
               sheetName = "referenceReadCounts",
               append = TRUE)
}

WriteListToXLSX <- function(listOfDataFrames, filepath = "list.xlsx") {
  sizeOfList <- length(listOfDataFrames)
  dataFrameNames <- names(listOfDataFrames)
  index <- 1
  write.xlsx(listOfDataFrames[[index]],
             filepath,
             sheetName = dataFrameNames[index])
  for(i in (sizeOfList - 1)) {
    index <- i + 1
    write.xlsx(listOfDataFrames[[index]],
               filepath,               
               sheetName = dataFrameNames[index],
               append = TRUE)
  }
}

ReadXLSXToList <- function(filepath) {
  inputList <- list(reference = unfactorize(read.xlsx(filepath, header = FALSE, sheetName = "reference")),
                    sample = unfactorize(read.xlsx(filepath, header = FALSE, sheetName = "sample")))
#  names(inputList) <- c("reference", "sample")
#  inputList[] <- lapply(inputList, as.character)
  return(inputList)
}

# ReferenceBamFiles <- function(filepath) {
#   unfactorize(ReadXLSXToList(filepath)$reference)[,1]
# }
#  
# SampleBamFiles <- function(filepath) {
#   unfactorize(ReadXLSXToList(filepath)$sample)[,1]
# }
#   
unfactorize <- function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.character(df[[i]])
  return(df)
}


# 
# 
# 
##################################################################################
# functions for bootstrapping
##################################################################################
#a function that randomly samples positions from a vector using the bootstrap principle
BootstrapPositions <-  function(pos) {
    return(sample(pos, length(pos), replace = TRUE))
}
 
#a function that randomly samples positions from a vector using the subsampling principle
SubsamplingPositions <- function(pos, mtry = round(sqrt(length(pos)))) {
    return(sample(pos, mtry, replace = FALSE))
}

#calculates the median from a vector given specific positions
PositionMedian <- function(positon, vector) {
    return(median(vector[positon]))
}

#calculate the median for each gene
GenePositionMedian <- function(genesPos, vector) {
    medians = sapply(genesPos, PositionMedian, vector = vector)
    names(medians) = names(genesPos)
    return(medians)
}

RatioMatrix <- function(sampleMat, refMedian) {
    division <- function(a, b) {
        return(a/b)
    }
    ratioMatrix = apply(sampleMat, 2, division, b = refMedian)
    return(ratioMatrix)
}

#calculates the median from a vector given specific positions for a matrix
GeneMedianRatioMatrix <- function(genesPos, ratioMatrix) {
    # Package check complaints
    i <- NULL
    ratioList =  foreach(i = 1:ncol(ratioMatrix), .combine = cbind) %do% {
        GenePositionMedian(genesPos,ratioMatrix[, i])
    }
    return(ratioList)
}

ReferenceWeights <- function(refmat, varianceFunction = IQR) {
    variance = apply(refmat, 2, varianceFunction)
    weights = 1/variance
    return(weights)
}

# #split into a function generating the bootstrap distribution #and a function averaging over amplicons.
# define the function to generate bootstrap distributions and return a list with the genewise bootstrapping
BootList <- function(geneNames, sampleMatrix, refmat, replicates) {
  
  # get the genes positions in the matrix as a list from a gene name vector
  genesPos <- IndexGenesPositions(geneNames)
  
  # a vector and matrix are not the same and for a vector iterating over the column makes no sense so we have
  # to check if a mtrix or a vector was passed. ncol only works for matrix not for vector
  if (class(sampleMatrix) == "matrix") {
    iterator <- 1:ncol(sampleMatrix)
  } else {
    iterator <- 1
  }
  
  i <- NULL
  j <- NULL
  bootListSamples <- foreach(i = iterator) %:% foreach(j = rep(1, replicates), .combine = rbind) %dopar% {
    # a vector and matrix are not the same and for a vector iterating over the column makes no sense so we have
    # to chekc if a mtrix or a vector was passed.
    if (class(sampleMatrix) == "matrix") {
      testSample <- sampleMatrix[, i]
    } else {
      testSample <- sampleMatrix
    }
    
    
    # for each gene subsample the amplicon positions independently
    
    # sample the samples using bootstrapping
    sampleBootPos <- sample(1:ncol(refmat), ncol(refmat), replace = TRUE)
    
    
    geneBootPos <- c(lapply(genesPos, SubsamplingPositions), recursive = TRUE)
    
    
    # given the obtained sampling using the bootstraps calulated above
    bootRatio <- testSample[geneBootPos]/rowMedians(refmat[geneBootPos, sampleBootPos])
    
    
    # after the bootstrapping the gene positions in the vector changes so recalculate them
    adjustedLength <- function(a) {
      round(sqrt(length(a)))
    }
    
    splitClass <- rep(names(genesPos), sapply(genesPos, adjustedLength))
    
    newGenesPos <- split(seq_along(splitClass), splitClass)
    
    
    sapply(newGenesPos, PositionMedian, vector = bootRatio)
  }
  
  names(bootListSamples) <- basename(colnames(sampleMatrix))
  return(bootListSamples)
}

#calculate significance 
CheckSignificance <- function(bootList) {
    # TODO package check complains
    i <- NULL
    j <- NULL
    sigTables =  foreach(i = seq_along(bootList)) %:% 
        foreach(j = 1:ncol(bootList[[i]]),.combine = rbind) %do% {
            quantiles = quantile(bootList[[i]][,j],c(0.05,0.95))
            significant = (quantiles[1] > 1 & quantiles[2]  > 1) | (quantiles[1] < 1 & quantiles[2] < 1)
            c(quantiles,significant)
        }
    
    return(sigTables)
}


SignificantGenes <- function(sigList, genesPositionsIndex) {
    # TODO package check complains
    i <- NULL
    sigGenes = foreach(i = seq_along(sigList)) %do% {
        names(genesPositionsIndex)[sigList[[i]][, 3] == 1]
    }
    return(sigGenes)
}


# ##################################################################################
# # functions for background noise estimation
# ##################################################################################
NonSignificantGeneIndex <- function(sigList, genesPositionsIndex) {
    # TODO package check complains
    i <- NULL
    genePosNonSig = foreach(i = seq_along(sigList)) %do% {
        sigGenes = names(genesPositionsIndex)[sigList[[i]][, 3] == 1]
        RemSigGenes(genesPositionsIndex, sigGenes)
    }
    return(genePosNonSig)
}

AmplProbMultipeSamples <- function(genePosNonSig) {
    # TODO package check complains
    i <- NULL
    amplWeights = foreach(i = seq_along(genePosNonSig)) %do% {
        AmplProb(genePosNonSig[[i]]) 
    }
    return(amplWeights)
}

Background <- function(geneNames, samplesNormalizedReadCounts, referenceNormalizedReadCounts, bootList, replicates = 1000) {

  #which genes showed significant changes
  sigList <- CheckSignificance(bootList)
  # gene index
  genesPositionsIndex <- IndexGenesPositions(geneNames)
  # calculate the reference median
  #refMedian <- apply(referenceNormalizedReadCounts, 1, median)
  refMedian <- rowMedians(referenceNormalizedReadCounts)
  # calculate the ratio matrix for each sample
  ratioMatrix <- RatioMat(samplesNormalizedReadCounts, refMedian)
  # remove the significant genes from the noise estimation
  genesPosNonSig <- NonSignificantGeneIndex(sigList, genesPositionsIndex)
  # calculate the weight for each amplicon of significant genes
  amplWeights <- AmplProbMultipeSamples(genesPosNonSig)
  uniqueAmpliconNumbers <- NumberOfUniqueAmplicons(genesPositionsIndex)
  
  i <- NULL
  backgroundObject <- foreach(i = seq_along(genesPosNonSig)) %dopar% {
    IterateAmplNum(uniqueAmpliconNumbers, ratioMatrix[unlist(genesPosNonSig[[i]]), i], replicates = replicates, 
                   probs = amplWeights[[i]])
  }
  return(backgroundObject)
}

round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
#  numeric_columns <- sapply(x, mode) == 'numeric'
  numeric_columns <- sapply(x, is.numeric)
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

ReportTables <- function(bootList,
                         geneNames,
                         backgroundNoise,
                         referenceNormalizedReadCounts,
                         samplesNormalizedReadCounts) {
  
  # get the background noise in a format that can be used for a report table
  backgroundReport <- BackgroundReport(backgroundNoise, geneNames)
  
  
  # gene index
  genesPositionsIndex <- IndexGenesPositions(geneNames)
  
  # calculate the reference median
  #refMedian <- apply(referenceNormalizedReadCounts, 1, median)
  refMedian <- rowMedians(referenceNormalizedReadCounts)
  
  # calculate the ratio matrix for each sample
  ratioMatrix <- RatioMat(samplesNormalizedReadCounts, refMedian)
  
  # calculate the genewise ratio matrix from the ratio_mat
  ratioMatGene <- GeneMedianRatioMatrix(genesPositionsIndex, ratioMatrix)
  
  sigList <- CheckSignificance(bootList)
  # TODO because package generaton complains..
  i <- NULL
  reportTables <- foreach(i = seq_along(backgroundReport)) %do% {
    isSig <- (sigList[[i]][, 3] < 1) | (sigList[[i]][, 1] > 1)
    backgroundUp <- 1 * (backgroundReport[[i]][, 2] - 1)
    backgroundDown <- 1 * (1 - backgroundReport[[i]][, 1])
    aboveNoise <- (ratioMatGene[, i] > 1 & (ratioMatGene[, i] - 1) > backgroundUp) | (ratioMatGene[, i] < 
                                                                                        1 & (1 - ratioMatGene[, i]) > backgroundDown)
    up <- apply(bootList[[i]], 2, PercentAboveValue, value = 1.5)
    down <- apply(bootList[[i]], 2, PercentBelowValue, value = 0.5)


#     data.frame(medianRatio = ratioMatGene[, i],
#                up = up,
#                down = down,
#                fivePercentQuantile = sigList[[i]][, 1],
#                fiftyPercentQuantile = sigList[[i]][, 2],
#                ninetyFivePercentQuantile = sigList[[i]][, 3],
#                fivePercentNoise = backgroundReport[[i]][, 1],
#                ninetyfivePercentNoise = backgroundReport[[i]][, 2],
#                significant = isSig,
#                aboveNoise = aboveNoise,
#                amplNum = table(geneNames),
#                passed = isSig + aboveNoise)

    dfTemp <- data.frame(medianRatio = ratioMatGene[, i],
               up = up,
               down = down,
               fivePercentQuantile = sigList[[i]][, 1],
               fiftyPercentQuantile = sigList[[i]][, 2],
               ninetyFivePercentQuantile = sigList[[i]][, 3],
               fivePercentNoise = backgroundReport[[i]][, 1],
               ninetyfivePercentNoise = backgroundReport[[i]][, 2],
               significant = isSig,
               aboveNoise = aboveNoise,
               amplNum = as.vector(table(geneNames)),
               passed = isSig + aboveNoise)

    significativeNumbers <- 2

    dfTemp <- round_df(dfTemp, significativeNumbers)

#     names(dfTemp) <- c("MedRat",
#                       "Up",
#                       "Down",
#                       "5%Qte",
#                       "50%Qte",
#                       "95%Qte",
#                       "5%Nse",
#                       "95%Nse",
#                       "Signif.",
#                       "AbvNse",
#                       "ANum",
#                       "Psd")


names(dfTemp) <- c("MedianRatio",
                   "Up",
                   "Down",
                   "5% Quantile",
                   "50% Quantile",
                   "95% Quantile",
                   "5% Noise",
                   "95% Noise",
                   "Signif.",
                   "AboveNoise",
                   "Amplicons",
                   "Passed")

    dfTemp
  }
  names(reportTables) <- seq_along(reportTables)
  return(reportTables)
}

BackgroundReport <- function(background, geneNames) {
    # TODO because package generaton complains..
    i <- NULL
    j <- NULL
    backgroundReport = foreach(i = seq_along(background)) %:% 
        foreach(j = table(geneNames), .combine = rbind) %do% {         
          background[[i]][[as.character(j)]]
        }
}

#remove genes that were considered significant by the algorithm
RemSigGenes <- function(genesPos, sigGenes) {
    #if some of the genes were reported with a significantly different read count the should not be used for
    #background testing
    #which genes are in the genes pos object
    geneNames= names(genesPos)
    #define a function for not in
    `%ni%` = Negate(`%in%`) 
    #which of these genes did not show a significant read count  change
    nonSigGenes = which(geneNames %ni% sigGenes)
    #keep only the genes that did not show a significant read coount change
    nonSigGenesPos = genesPos[nonSigGenes]
    return(nonSigGenesPos)
}

MultiRemSigGenes <- function(sigList,genesPos) {
    i <- NULL
    #remove the significant genes from the noise estimation
    genePosNonSig = foreach(i = seq_along(sigList)) %do% {
        sigGenes = names(genesPos)[sigList[[i]][,3] == 1]
        RemSigGenes(genesPos, sigGenes)
    }
    return(genePosNonSig)
}

# #extract the postions from a ratios vector that belong to specified genes
# getGeneAmpl <- function(ratios,genes_pos) {
#   
#   posVector = unlist(genesPos)
#   posRatios = ratios[posVector]
#   
#   return(posRatios)
#   
# }

AmplProb <- function(genesPos) {
    #how many amplicons where used for each ofhte genes
    geneCounts = elementLengths(genesPos)
    #adjust the probablity for depending on the number of amplicons for each gene
    genePerc =  1/geneCounts
    #this information has to be available for each stable position
    ampliconProb = rep(genePerc,geneCounts)
    return(ampliconProb)
}

# #calculate the weight for each amplicon
# multiGetAmplProb <- function(genesPos) {
# 
#     amplWeights = foreach(i = seq_along(genePosNonSig)) %do%  {
#         AmplProb(genePosNonSig[[i]]) 
#     }
#     
#     return(amplWeights)
# }


RatioMat <- function(sampleMat, refMedian) {
  division <- function(a,b) {
    return(a/b)
  }
  ratioMatrix = apply(sampleMat, 2, division, b =refMedian)
  return(ratioMatrix)
}


SampleRatio <- function(ratios, numAmpl, amplWeights = NULL) {
    randomPos = sample(seq_along(ratios),numAmpl,prob = amplWeights)
    randomlySelectedRatios = ratios[randomPos]
    randomMedian = median(randomlySelectedRatios)
    return(randomMedian)
}

SampleNoiseGenes <- function(numAmpl = 2, ratios,replicates = 100, probs = NULL) {
  #now repeat the sampling for the selected number of amplicons replicates times
  sampleNoiseDistribution = replicate(replicates, SampleRatio(ratios, numAmpl, amplWeights = probs))
  #get the 5 and 95 percent quantiles of this noise distribution    
  distributionQuantiles = quantile(sampleNoiseDistribution,c(0.05,0.95))
  return(distributionQuantiles)
}


#how many unique amplicon numbers are there
NumberOfUniqueAmplicons <- function(genesPos) {
  # calcUniqueAmpliconNumbers <-function(genesPos) {
    ampliconNumbers = elementLengths(genesPos)
    uniqueAmpliconNumbers = sort(unique(ampliconNumbers))
    return(uniqueAmpliconNumbers)
}

#now calculate the background for each of the unique amplicon numbers
IterateAmplNum <- function(uniqueAmpliconNumbers, ratios, replicates = 100, probs = NULL) {
    # TODO Needed because of package check complaints..
    i <- NULL
    noiseResults =  foreach(i = seq_along(uniqueAmpliconNumbers)) %do% { 
        sampledNoise = SampleNoiseGenes(uniqueAmpliconNumbers[i], ratios = ratios,replicates = replicates, probs = probs)
        names(sampledNoise) = c("5%","95%")  
        sampledNoise
    }
    names(noiseResults) = as.character(uniqueAmpliconNumbers)    
    return(noiseResults)
}


# #check function not sure what it does
# #multiIterateAmplNum <- function(unique_ampliconNumbers, genePosNonSig, ratioMatrix, replicates = 100, probs = NULL) {
#   
#   
# #   backgroundObject = foreach(i = seq_along(genePosNonSig)) %do% {
#     
#     
# #     IterateAmplNum(uniqueAmpliconNumbers, ratioMatrix[unlist(genePosNonSig[[i]]), i], replicates = replicates, probs = probs[[i]])
#     
#     
# #   }
#  
# #   return(backgroundObject)
#    
# #}
# 
# 
# #######################################################################################
# ## functions for report generation
# #######################################################################################
# SelectNoiseInfo <- function(amplNum,noiseResults) {
#     noise =  foreach(i = seq_along,.combine = rbind) %do% {
#         noiseResults[[as.character(amplNum[i])]] 
#     }     
#     return(noise)
# }
# 
# #define a function to generate the results tables
# getResultTables <- function(bootListSamples,genesPos,samplemat,refmat,reps,genes){
#   
#   
#   resultTables =  foreach(i = seq_along(bootListSamples)) %dopar% {
#     
#     bootFive = apply(bootListSamples[[i]],2,quantile,0.05)
#     bootNinetyFive = apply(bootListSamples[[i]],2,quantile,0.95)
#     bootMedian = apply(bootListSamples[[i]],2,median)
#     
#     calcMedian = geneMedianRatio(genesPos,ratio)
#     
#     
#     if(class(samplemat) == "matrix") {
#       
#       testsample = samplemat[,i]
#       
#     }else{
#       
#       testsample = samplemat
#       
#     }     
#     
#     medianReference = rowMedians(refmat)
#     ratio =   testsample/medianReference
#     
#     #diff_sig_result = diff_sig_test(genes_pos,ratio,reps,genes)
#     noiseInfo =  SelectNoiseInfo(amplNum,noiseResults)
#     
#     
#     
#     isSigBoot = (bootNinetyFive < 1 & bootMedian  < 1) | (bootFive > 1 & bootMedian > 1)
#     isSigSampling = (diffSigResult[,"5%"] > bootMedian & bootMedian  < 1) | (diffSigResult[,"95%"] < bootMedian & bootMedian > 1)
#     
#     isSig = isSigBoot & isSigSampling & (bootMedian < 0.75 | bootMedian > 1.25)
#     
#     cbind(bootFive,bootNinetyFive,bootMedian,diffSigResult,isSig,table(genes)) 
#     
#   }
#   
# }
# 
# 
# ####################################################################################################################################
# ## Visualize the results
# NOT BEING USED AT THIS POINT .. MAYBE IN FUTURE..
# ####################################################################################################################################
# bedToGeneOrder <- function(gr) {
#   
#   geneNames <- gr$geneNames
#   chromosomes <- as.character(seqnames(gr)) 
#   starts <-  start(gr)
#   
#   geneStart <- aggregate(starts,list(geneNames),min)[,2]
#   geneChromosome <- aggregate(chromosomes,list(geneNames),unique)[,2]
#   
#   chromosomes <- gsub("chr", "",geneChromosome,ignore.case = TRUE)
#   chromosomes <- gsub("X", "22",chromosomes,ignore.case = TRUE)
#   chromosomes <- as.numeric(gsub("Y", "23",chromosomes,ignore.case = TRUE))
#   divider =  10 ^max(nchar(geneStart)) 
#   geneOrderer <- chromosomes + geneStart/divider
#   
#   return(rank(geneOrderer))  
# }

PlotBootstrapDistributions <- function(bootList, reportTables, outputFolder = getwd(), sampleNames = NULL, save = FALSE) {

  i <- NULL
  plotList <- foreach(i = seq_along(bootList)) %do% {
    #for (i in seq_along(bootList)) {
    selSample <- i
    test <- as.factor(reportTables[[selSample]][, "Passed"])
    test <- revalue(test, c(`0` = "noChange", `1` = "nonReliableChange", `2` = "ReliableChange"))

    
    ratios <- NULL
    testsPassed <- NULL
    
    df <- data.frame(class = as.factor(colnames(bootList[[selSample]])), ratios = (as.vector(t(bootList[[selSample]]))), 
                     testsPassed = test)
    
    
    #if a genomic ranges object has been 
    #if(!is.null(gr)) {
    # newGeneOrder <- bedToGeneOrder(gr)
    # df$class = with(df,factor(class,levels(class)[newGeneOrder]))
    #}
    
    
    ylim1 <- boxplot.stats((df$ratios))$stats[c(1, 5)]
    ylim1 <- c(max(ylim1),max(ylim1))
    ylim1 <- log(ylim1)
    ylim1 <- ylim1 * c(-7,7)
    
    if(is.null(sampleNames)) {
      
      filename <- names(bootList[selSample])
      
    }else{
      
      filename <- sampleNames[selSample]
      
    } 
    
    
    
    test <- ggplot(df, aes(x = class, y = log(ratios), fill = testsPassed)) + geom_boxplot() + ggtitle(filename) + 
      theme(plot.title = element_text(lineheight = 0.8, face = "bold"), text = element_text(size = 15), 
            axis.text.x = element_text(angle = 90)) + scale_fill_manual(name = "CNV Reliability", values = c("#56B4E9", "#CC79A7" ,"#D55E00")
                                                                        , labels = c("0" = "Foo", "1" = "Bar")) + coord_cartesian(ylim = ylim1) + theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold")) + 
      scale_x_discrete("Gene Names") + geom_hline(yintercept=log(1.5), color="#009E73") + geom_hline(yintercept=log(0.5), color="#009E73") + geom_hline(yintercept=0, color="#009E73") 
    
    if(save == TRUE) {
      
      dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
      
      
      filepath <- paste0(outputFolder, "/", filename, "_plot.pdf")
      
      
      ggsave(filename = filepath, plot = test, height = 7.42 * 1, width = 8.11 * 2, limitsize = FALSE)
      
    }
    
    return(test)
    
  }
  
  return(plotList)
  
}

#
# PlotBootstrapDistributions <- function(bootList, reportTables, outputFolder) {
#   print("BEGIN of PlotBootstrapDistributions")
#   # TODO Needed because of package check complaints..
#   i <- NULL
#   for(i in seq_along(bootList)) {
#     print("inside cycle")
#     selSample = i
#     test = as.factor(reportTables[[selSample]][, "passed"])
#     test = revalue(test,c("0" = "noChange","1" =  "significantChange","2" = "aboveNoise"))
#     
#     # Package check complaints...
#     #        ratios <- NULL
#     #        testsPassed <- NULL
#     
#     df = data.frame(class = as.factor(colnames(bootList[[selSample]])), ratios = (as.vector(t(bootList[[selSample]]))), testsPassed = test)
#     
#     print("STEP 0")
#     
#     filename <- names(bootList[selSample])
#     
#     print("STEP 0.5")
#     #        print(df)
#     #        test =  ggplot(df)
#     # test = qplot(Sepal.Length, Petal.Length, data = iris, color = Species)
#     print("STEP 0.6")
#     test =  ggplot(df, aes(x=class, y=log(ratios), fill=testsPassed)) +
#       geom_boxplot() + 
#       ggtitle(filename) +
#       theme(plot.title = element_text(lineheight=.8, face="bold"), text = element_text(size=15), axis.text.x = element_text(angle=90, vjust=1))
#
#     #     test + ggtitle(filename)
#     
#     print("STEP 1")
#     print(outputFolder)
#     print(filename)
#     print("STEP 1.1")
#     dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
#     
#     filepath = paste0(outputFolder, "/",filename, "_plot.pdf")
#     print(paste0("filepath : ", filepath))
#     
#     print("STEP 2")
#     
#     #     ggsave(filename = paste0(i,"_plot.pdf"),plot = test,height = 7.42 * 1,width = 8.11 * 2,limitsize = FALSE)
#     ggsave(filename=filepath, plot=test, height=7.42 * 1, width=8.11 * 2, limitsize=FALSE)
#     
#     print("STEP 3")
#     
#   }
#   print("END of PlotBootstrapDistributions")
# }
#
