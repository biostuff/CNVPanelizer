###############################################################################
# functions to read the read counts
###############################################################################

BedToGenomicRanges <- function(panelBedFilepath,
                               ampliconColumn,
                               split = "_",
                               doReduce = TRUE,
                               rangeExtend = 0,
                               dropChromossomes = NA,
                               skip = 1) {
  #load the bed file
  segments <- read.table(panelBedFilepath,
                         sep = "\t",
                         as.is = TRUE,
                         skip = skip)
  gr <- GRanges(segments[, 1], IRanges(segments[, 2], segments[, 3]))

  #make sure that all regions are more than one read apart
  GenomicRanges::start(gr) <- GenomicRanges::start(gr) - rangeExtend
  GenomicRanges::end(gr) <- GenomicRanges::end(gr) + rangeExtend

  if(doReduce) {
    reducedGr <- GenomicRanges::reduce(gr, ignore.strand = TRUE)

    #get the genes for each reduced range
    hitsAnnotation <- GenomicRanges::findOverlaps(reducedGr, gr,
                                                  select = "first",
                                                  ignore.strand=TRUE)
    #get the amplicon names from the segments
    amplicons = segments[hitsAnnotation, ampliconColumn]
    #replace the original gr with the reduced gr
    gr <- reducedGr
  } else{
    amplicons = segments[ , ampliconColumn]
  }

  #get the genes from segments
  splitted = strsplit(amplicons, split = split)
  # Required because of package checking complaints
  i <- NULL
  genes = foreach(i = seq_along(splitted)) %do% {
    splitted[[i]][1]
  }

  genes = unlist(genes)
  # missing the dollar symbol before text because of groovy script
  elementMetadata(gr)["geneNames"] = genes
  elementMetadata(gr)["ampliconNames"] = paste(seqnames(gr),
                                               start(gr),
                                               end(gr),
                                               sep = "_",
                                               amplicons)

  if (!is.na(dropChromossomes)) {
    gr <- dropSeqlevels(gr, dropChromossomes)
  }
  return(gr)
}

#load the files you want to analyze
ReadCountsFromBam <- function(bamFilenames,
                              sampleNames,
                              gr,
                              ampliconNames,
                              minimumMappingQuality = 20,
                              removeDup = FALSE) {

  if (missing(sampleNames)) {
    sampleNames = bamFilenames
  }

  if (missing(ampliconNames)) {
    ampliconNames = elementMetadata(gr)$ampliconNames
  }

  bamIndexFilepaths <- BamIndexFilepaths(bamFilenames)

  # Because of package check complaints
  i <- NULL
  curbam = foreach(i = seq_along(bamFilenames), .combine = cbind) %do% {

    if (!file.exists(bamIndexFilepaths[i])) {
      message(paste("Bai file not found. Creating bai file for", bamFilenames[i], "..."))
      IndexMultipleBams(bamFilenames)
    }

    message(paste("Reading counts for bam file: ", bamFilenames[i]))
    countBamInGRanges(bamFilenames[i],
                      gr,
                      remove.dup = removeDup,
                      min.mapq = minimumMappingQuality,
                      get.width = TRUE)
  }

  listOfBamsHasOnlyOneElement <- length(bamFilenames) == 1
  if (listOfBamsHasOnlyOneElement) {
    curbam <- matrix(curbam)
  }

  colnames(curbam) = sampleNames
  rownames(curbam) = ampliconNames
  return(curbam)
}

#get the combined counts
CombinedNormalizedCounts <- function (sampleCounts,
                                      referenceCounts,
                                      method = "tmm",
                                      ampliconNames = NULL) {
  sampleCounts <- as.matrix(sampleCounts)
  referenceCounts <- as.matrix(referenceCounts)
  allCounts <- cbind(sampleCounts, referenceCounts)
  classes <- rep(c("samples", "reference"), c(ncol(sampleCounts),
                                              ncol(referenceCounts)))
  bamDataRangesNorm <- NormalizeCounts(allCounts, method = method)
  #  bamDataRangesNorm <- NormalizeCounts(allCounts)
  if (!is.null(ampliconNames)) {
    rownames(bamDataRangesNorm) = ampliconNames
  }
  normalizedSamples <- bamDataRangesNorm[, which(classes == "samples"), drop = FALSE]
  normalizedReference <- bamDataRangesNorm[, which(classes == "reference"), drop = FALSE]
  return(list(samples = normalizedSamples, reference = normalizedReference))
}

###############################################################################
# helper functions
###############################################################################

VerifiyIfOutputDirectoryExistsOrIsNotEmpty <- function(outputDir) {
  outputDirectoryExists <- file.exists(outputDir)
  outputDirectoryIsNotEmpty <- length(dir(outputDir, all.files=TRUE)) != 0
  if (outputDirectoryExists || outputDirectoryIsNotEmpty) {
    stop("Output directory already exists or is not empty. Please delete or change output directory")
  }
}

tss <- function(allCounts) {
  return(apply(allCounts, 2, function(x) {x / sum(x)}))
}

NormalizeCounts <- function(allCounts, method = "tmm") {
  # TODO check if normalization method is supported
  allowedNormalizationMethods <- c("tmm", # trimmed mean of M-values
                                   "tss") # total sum scaling
  tinyValueToAvoidZeroReadCounts <- 1e-05
  bamDataRangesNorm <- NULL
  if (method == "tmm") {
    bamDataRangesNorm <- tmm(allCounts, refColumn = NULL, k = tinyValueToAvoidZeroReadCounts)
  } else if (method == "tss") {
    #    bamDataRangesNorm <- apply(allCounts, 2, function(x) {x / sum(x)})
    bamDataRangesNorm <- tss(allCounts)
  } else {
    message(paste("method", method, "not included in the allowed Normalization methods", allowedNormalizationMethods))
  }
  return(bamDataRangesNorm)
}

#get the positions for each gene as a list
IndexGenesPositions = function(genes) {
  positions = seq_along(genes)
  genesPos = split(positions, genes)
  return(genesPos)
}

# #how many numbers in a vector are above a given cutoff
# PercentAboveValue <- function(vector, value) {
#     sum(vector > value)/length(vector)
# }
#
# #how many numbers in a vector are below a given cutoff
# PercentBelowValue <- function(vector, value) {
#     sum(vector < value)/length(vector)
# }

BamIndexFilepaths <- function(bamFilepaths, indexType = ".bam.bai", ignoreCase = FALSE) {
  return (gsub(".bam$",
               bamFilepaths,
               replacement = indexType,
               ignore.case = ignoreCase))
}

# check this params with Thomas
#index the bam files if there is no index yet
IndexMultipleBams <- function(bams,
                              index_type = ".bam.bai") {
  #check if the index already exists and need to be indexed
  #     potentialBaiFilenames <- gsub(".bam$",
  #                                 bams,
  #                                 replacement = index_type,
  #                                 ignore.case = TRUE)

  potentialBaiFilenames <- BamIndexFilepaths(bams, indexType = index_type)

  print(paste("potentialBaiFilenames", potentialBaiFilenames))

  bamsToBeIndexed <- bams[!sapply(potentialBaiFilenames, file.exists)]
  if(length(bamsToBeIndexed) > 0) {
    #index the bams
    #if(multicore) {
    # check with Thomas if we can replace by this..
    #  bplapply(bamsToBeIndexed, indexBam)
    #mclapply(bamsToBeIndexed, indexBam, mc.cores = ncores)
    #} else {
    lapply(bamsToBeIndexed, indexBam)
    #}
  }
}

WriteListToXLSX <- function(listOfDataFrames, multipleFiles = FALSE, outputFolder = file.path(getwd(), "xlsx"), filepath = "list.xlsx") {
  if (multipleFiles) {
    dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
    for(name in names(listOfDataFrames)){
      message(paste0("Saving file to '", file.path(outputFolder, paste0(name, ".xlsx"))))
      write.xlsx(listOfDataFrames[[name]], file = file.path(outputFolder, paste0(name, ".xlsx")), rowNames = TRUE)
    }
  } else {
    message(paste0("Saving file to '", filepath, "'"))
    write.xlsx(listOfDataFrames, file = filepath, rowNames = TRUE)
  }
}

ReadXLSXToList <- function(filepath, rowNames = TRUE, colNames = TRUE) {
  listOfDataFrames <- list()
  for(name in getSheetNames(filepath)) {
    listOfDataFrames[[name]] <- read.xlsx(filepath,
                                          rowNames = rowNames,
                                          colNames = colNames,
                                          sheet = name)
  }
  return(listOfDataFrames)
}

###############################################################################
# functions for bootstrapping
###############################################################################

adjustedLength <- function(a) {
  round(sqrt(length(a)))
}

ratiosMean <- function(ratios) {
  return(exp(mean(log(ratios))))
}

#a function that randomly samples positions from a vector using subsampling
SubsamplingPositions <- function(pos, mtry = adjustedLength(pos)) {
  return(sample(pos, mtry, replace = FALSE))
}

#calculates the mean from a vector given specific positions
PositionMean <- function(position, vector) {
  return(ratiosMean(vector[position]))
}

#calculate the mean for each gene
GenePositionMean <- function(genesPos, vector) {
  means = sapply(genesPos, PositionMean, vector = vector)
  names(means) = names(genesPos)
  return(means)
}

#calculates the mean from a vector given specific positions for a matrix
GeneMeanRatioMatrix <- function(genesPos, ratioMatrix) {
  # Package check complaints
  i <- NULL
  ratioList =  foreach(i = 1:ncol(ratioMatrix), .combine = cbind) %do% {
    GenePositionMean(genesPos,ratioMatrix[, i])
  }
  ratioList <- as.matrix(ratioList)
  colnames(ratioList) <- colnames(ratioMatrix)
  return(ratioList)
}

# ReferenceWeights <- function(refmat, varianceFunction = sd) {
#     variance = apply(refmat, 2, varianceFunction)
#     weights = 1/variance
#     return(weights)
# }

# split into a function generating the bootstrap distribution
# and a function averaging over amplicons.
# define the function to generate bootstrap distributions and
# return a list with the genewise bootstrapping
BootList <- function(geneNames, sampleMatrix, refmat, replicates) {
  enforceDeterministicResult()
  # get the genes positions in the matrix as a list from a gene name vector
  genesPos <- IndexGenesPositions(geneNames)

  if (class(sampleMatrix) != "matrix") {
    stop(paste("Parameter 'sampleMatrix' has to be of class matrix and is of class", class(sampleMatrix)))
  }

  if (is.null(colnames(sampleMatrix)) | length(colnames(sampleMatrix))!=ncol(sampleMatrix)) {
    stop("All columns of 'sampleMatrix' have to be named")
  }

  # a vector and matrix are not the same and for a vector iterating over
  # the column makes no sense so we have
  # to check if a matrix or a vector was passed. ncol only works for matrix
  # not for vector
  if (class(sampleMatrix) == "matrix") {
    iterator <- 1:ncol(sampleMatrix)
    # } else {
    #     iterator <- 1
  }

  i <- NULL
  j <- NULL
  bootListSamples <- foreach(i = iterator) %:%
    foreach(j = rep(1, replicates), .combine = rbind) %do% {
      # a vector and matrix are not the same and for a vector
      # iterating over the column makes no sense so we have to check
      # if a mtrix or a vector was passed.
      if (class(sampleMatrix) == "matrix") {
        testSample <- sampleMatrix[, i]
        # } else {
        #     testSample <- sampleMatrix
      }
      # for each gene subsample the amplicon positions independently
      # sample the samples using bootstrapping
      sampleBootPos <- sample(1:ncol(refmat),
                              ncol(refmat),
                              replace = TRUE)
      geneBootPos <- c(lapply(genesPos,
                              SubsamplingPositions),
                       recursive = TRUE)
      # given the obtained sampling using the bootstraps calculated
      refMatPos <- refmat[geneBootPos, sampleBootPos, drop=FALSE]

      bootRatio <- testSample[geneBootPos]/rowMeans(refMatPos)

      # after the bootstrapping the gene positions in the vector
      # changes so recalculate them
      splitClass <- rep(names(genesPos), sapply(genesPos,
                                                adjustedLength))
      newGenesPos <- split(seq_along(splitClass), splitClass)
      sapply(newGenesPos, PositionMean, vector = bootRatio)
    }

  names(bootListSamples) <- basename(colnames(sampleMatrix))

  # In case a single step of bootstrap is done (for testing purposes),
  # the elements of each sample will be converted to a matrix
  for (i in seq_along(bootListSamples)) {
    bootListSamples[[i]] <- matrix(bootListSamples[[i]], ncol = length(names(genesPos)))
    colnames(bootListSamples[[i]]) <- names(genesPos)
  }
  return(bootListSamples)
}



ConsensusCheck <- function(genomicRangesFromBed, sampleMatrix, referenceMatrix) {
  #  sampleMatrix <- samplesNormalizedReadCounts
  #  referenceMatrix <- referenceNormalizedReadCounts

  #  View(sampleMatrix[, c(1,2)])
  #  View(apply(referenceMatrix, 1, median))

  #  ConsensusCheck(genomicRangesFromBed, results_1PGM_sequencingLibraries1@sampleReadCounts, results_1PGM_sequencingLibraries1@referenceReadColunts),

  #  sampleMatrix <- results_1PGM_sequencingLibraries1@sampleReadCounts
  #  referenceMatrix <- results_1PGM_sequencingLibraries1@referenceReadColunts


  # This should not reach here.. but just in case..
  colnames(sampleMatrix) <- basename(colnames(sampleMatrix))
  colnames(referenceMatrix) <- basename(colnames(referenceMatrix))

  geneNames <- genomicRangesFromBed$geneNames

  allSampleRatios <- data.frame()
  for (smpl in colnames(sampleMatrix)) {
    #    smpl <- colnames(sampleMatrix)[which(grepl("28208", colnames(sampleMatrix)))]
    ampliconsRatio <- sampleMatrix[, smpl]/apply(referenceMatrix, 1, median)
    allSampleRatios <- rbind(allSampleRatios, ampliconsRatio)
  }

  allSampleRatios <- t(allSampleRatios)
  rownames(allSampleRatios) <- rownames(sampleMatrix)
  colnames(allSampleRatios) <- colnames(sampleMatrix)

  genesPositions <- IndexGenesPositions(geneNames)

  geneConsensus <- data.frame(matrix(nrow=length(genesPositions), ncol=ncol(sampleMatrix)))
  geneConsensus <- data.frame()
  for (smpl in colnames(allSampleRatios)) {
    for (gene in names(genesPositions)) {
      geneAmpliconsRatios <- allSampleRatios[,smpl][genesPositions[[gene]]]
      ampliconConsensus <- all(geneAmpliconsRatios>1) | all(geneAmpliconsRatios<1)
      geneConsensus[gene, smpl] <- ampliconConsensus
    }
  }
  # TODO case the reference has amplicons with 0 read counts, it returns NA instead of boolean
  geneConsensus[is.na(geneConsensus)] <- FALSE
  return(geneConsensus)
}

HaveMininumNumberOfAmplicons <- function(genomicRangesFromBed, numberOfAmplicons) {
  geneNames <- genomicRangesFromBed$geneNames
  genesPositions <- IndexGenesPositions(geneNames)
  sapply(genesPositions, function(x) {length(x)>=numberOfAmplicons})
}

#calculate significance
CheckSignificance <- function(bootList, significanceLevel = 0.05) {
  # Bonferroni correction
  significanceLevel <- significanceLevel/ncol(bootList[[1]])
  margin <- significanceLevel/2

  # package check complains
  i <- NULL
  j <- NULL
  sigTables =  foreach(i = seq_along(bootList)) %:%
    foreach(j = 1:ncol(bootList[[i]]), .combine = rbind) %do% {
      bootRatioDistribution <- bootList[[i]][,j]
      #      meanNoise <- exp(mean(log(bootRatioDistribution)))
      meanNoise <- ratiosMean(bootRatioDistribution)
      lowerBound <- quantile(bootRatioDistribution, margin, type = 1)
      upperBound <- quantile(bootRatioDistribution, 1 - margin, type = 1)
      bounds <- c(lowerBound, meanNoise, upperBound)
      roundedBounds <- round(bounds, digits = 2)
      maybeAmplification <- roundedBounds[1] > 1
      maybeDeletion <- roundedBounds[3] < 1
      roundedSignificant <- maybeAmplification | maybeDeletion
      copyNumberPutativeStatus <- "Normal"
      if(maybeAmplification) {
        copyNumberPutativeStatus <- "Amplification"
      } else {
        if (maybeDeletion) {
          copyNumberPutativeStatus <- "Deletion"
        }
      }
      sigTable <- c(roundedBounds, roundedSignificant, copyNumberPutativeStatus)
      names(sigTable)  <- c("lowerBound", "mean", "upperBound", "isSig", "putativeStatus")
      sigTable
    }
  names(sigTables) <- names(bootList)

  # name the genes in the list. All samples should share the same gene names
  geneNames <- colnames(bootList[[1]])
  for (i in seq_along(sigTables)) {
    rownames(sigTables[[i]]) <- geneNames
  }
  return(sigTables)
}

#SignificantGenes <- function(sigList, genesPositionsIndex) {
#  # package check complains
#  i <- NULL
#  sigGenes = foreach(i = seq_along(sigList)) %do% {
#    names(genesPositionsIndex)[sigList[[i]][, "isSig"] == TRUE]
#  }
#  return(sigGenes)
#}

# #############################################################################
# # functions for background noise estimation
# #############################################################################
NonSignificantGeneIndex <- function(sigList, genesPositionsIndex) {
  # package check complains
  i <- NULL
  genePosNonSig = foreach(i = seq_along(sigList)) %do% {
    sigGenes = names(genesPositionsIndex)[sigList[[i]][, "isSig"] == TRUE]
    selectedIndex <- RemSigGenes(genesPositionsIndex, sigGenes)
    if (length(selectedIndex) == 0) {
      genesPositionsIndex
    } else {
      selectedIndex
    }
  }
  names(genePosNonSig) <- names(sigList)
  return(genePosNonSig)
}

AmplProbMultipeSamples <- function(genePosNonSig) {
  # package check complains
  i <- NULL
  amplWeights = foreach(i = seq_along(genePosNonSig)) %do% {
    AmplProb(genePosNonSig[[i]])
  }
  names(amplWeights) <- names(genePosNonSig)
  return(amplWeights)
}

enforceDeterministicResult <- function() {
  set.seed(1)
}

Background <- function(geneNames,
                       samplesNormalizedReadCounts,
                       referenceNormalizedReadCounts,
                       bootList,
                       replicates = 1000,
                       significanceLevel = 0.05,
                       robust = FALSE) {
  enforceDeterministicResult()
  #which genes showed significant changes
  sigList <- CheckSignificance(bootList)
  # gene index
  genesPositionsIndex <- IndexGenesPositions(geneNames)
  # calculate the reference mean
  #refMean <- apply(referenceNormalizedReadCounts, 1, mean)
  refMean <- rowMeans(referenceNormalizedReadCounts)
  # calculate the ratio matrix for each sample
  ratioMatrix <- RatioMatrix(samplesNormalizedReadCounts, refMean)
  # remove the significant genes from the noise estimation
  genesPosNonSig <- NonSignificantGeneIndex(sigList, genesPositionsIndex)
  # calculate the weight for each amplicon of non significant changed genes
  amplWeights <- AmplProbMultipeSamples(genesPosNonSig)
  uniqueAmpliconNumbers <- NumberOfUniqueAmplicons(genesPositionsIndex)

  sampleIndex <- NULL
  backgroundObject <- foreach(sampleIndex = seq_along(genesPosNonSig)) %do% {
    message(paste0("Calculating Background for ", colnames(samplesNormalizedReadCounts)[sampleIndex]))
    nonSigAmpliconRatios <- ratioMatrix[unlist(genesPosNonSig[[sampleIndex]]), sampleIndex]
    IterateAmplNum(uniqueAmpliconNumbers,
                   nonSigAmpliconRatios,
                   replicates = replicates,
                   probs = amplWeights[[sampleIndex]],
                   significanceLevel = significanceLevel,
                   robust = robust)
  }
  return(backgroundObject)
}

roundDf <- function(x, digits) {
  # round all numeric variables
  # x: data frame
  # digits: number of digits to round
  #  numeric_columns <- sapply(x, mode) == 'numeric'
  numeric_columns <- sapply(x, is.numeric)
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

# Helper function to add column to reportTables
ReliableStatus <- function(putative, numberOfThresholdsPassed) {
  status <- as.vector(putative)
  #  print(status)
  if (numberOfThresholdsPassed < 2) {
    status <- "Normal"
  }
  return(status)
}

ReportTables <- function(geneNames,
                         samplesNormalizedReadCounts,
                         referenceNormalizedReadCounts,
                         bootList,
                         backgroundNoise) {
  # get the background noise in a format that can be used
  # for a report table
  backgroundReport <- BackgroundReport(backgroundNoise, geneNames)
  # gene index
  genesPositionsIndex <- IndexGenesPositions(geneNames)
  # calculate the reference mean
  #refMean <- apply(referenceNormalizedReadCounts, 1, mean)
  refMean <- rowMeans(referenceNormalizedReadCounts)
  # calculate the ratio matrix for each sample
  ratioMatrix <- RatioMatrix(samplesNormalizedReadCounts, refMean)
  # calculate the genewise ratio matrix from the ratio_mat
  ratioMatGene <- GeneMeanRatioMatrix(genesPositionsIndex, ratioMatrix)
  sigList <- CheckSignificance(bootList)
  # because package generaton complains..
  i <- NULL
  reportTables <- foreach(i = seq_along(backgroundReport)) %do% {
    #isSig <- (sigList[[i]][, "upperBound"] < 1) | (sigList[[i]][, "lowerBound"] > 1)
    putativeStatus <- sigList[[i]][, "putativeStatus"]
    isSig <- putativeStatus != "Normal"

    backgroundUp <- backgroundReport[[i]][, "UpperNoise"]
    backgroundDown <- backgroundReport[[i]][, "LowerNoise"]
    rMatGene <- ratioMatGene[, i]
    # TODO is this correct?! why compare rMatGene > 1 ?
    aboveNoise <- (rMatGene > 1 & sigList[[i]][, "lowerBound"] > backgroundUp) |
      (rMatGene < 1 & sigList[[i]][, "upperBound"] < backgroundDown)

    dfTemp <- data.frame(meanRatio = ratioMatGene[, i],
                         lowerBoundBootstrapRatio = sigList[[i]][, "lowerBound"],
                         meanBootstrapRatio = sigList[[i]][, "mean"],
                         upperBoundBootstrapRatio = sigList[[i]][, "upperBound"],
                         #                          lowerNoise = backgroundReport[[i]][, 1],
                         #                          meanNoise = backgroundReport[[i]][, 2],
                         #                          upperNoise = backgroundReport[[i]][, 3],
                         lowerNoise = backgroundReport[[i]][, "LowerNoise"],
                         meanNoise = backgroundReport[[i]][, "MeanNoise"],
                         upperNoise = backgroundReport[[i]][, "UpperNoise"],
                         significant = isSig,
                         aboveNoise = aboveNoise,
                         amplNum = as.vector(table(geneNames)),
                         putativeStatus = putativeStatus,
                         passed = isSig + aboveNoise)

    significativeNumbers <- 2
    dfTemp <- roundDf(dfTemp, significativeNumbers)
    names(dfTemp) <- c("MeanRatio",
                       "LowerBoundBoot",
                       "MeanBoot",
                       "UpperBoundBoot",
                       "LowerNoise",
                       "MeanNoise",
                       "UpperNoise",
                       "Signif.",
                       "AboveNoise",
                       "Amplicons",
                       "PutativeStatus",
                       "Passed")
    dfTemp
  }

  # TODO Check if this column added is correct
  for(i in 1:length(reportTables)) {
    tmpReportTable <- reportTables[[i]]
    #  tmpReportTable <- reportTables[[1]]
    status <- mapply(ReliableStatus, tmpReportTable$PutativeStatus, tmpReportTable$Passed)
    #  print(i)
    reportTables[[i]] <- cbind(reportTables[[i]], "ReliableStatus" = status)
  }

  #    names(reportTables) <- colnames(samplesNormalizedReadCounts)
  names(reportTables) <- names(bootList)
  return(reportTables)
}

BackgroundReport <- function(background, geneNames) {
  # because package generation complains..
  i <- NULL
  j <- NULL
  backgroundReport = foreach(i = seq_along(background)) %:%
    foreach(j = table(geneNames), .combine = rbind) %do% {
      background[[i]][[as.character(j)]]
    }
}

#remove genes that were considered significant by the algorithm
RemSigGenes <- function(genesPos, sigGenes) {
  #if some of the genes were reported with a significantly different read count
  #should not be used for
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

AmplProb <- function(genesPos) {
  #how many amplicons where used for each ofhte genes
  geneCounts = countAmplicons(genesPos)
  #adjust the probablity for depending on the number of amplicons for each gene
  genePerc =  1/geneCounts
  #this information has to be available for each stable position
  ampliconProb = rep(genePerc, geneCounts)
  return(ampliconProb)
}

RatioMatrix <- function(sampleMat, refMean) {
  apply(sampleMat, 2, `/`, refMean)
}

SampleRatio <- function(ratios, numAmpl, amplWeights = NULL) {
  replace <- FALSE
  if(numAmpl >= length(ratios)) {
    replace <- TRUE
  }
  randomPos = sample(seq_along(ratios),
                     numAmpl,
                     prob = amplWeights,
                     replace = replace)
  randomlySelectedRatios = ratios[randomPos]
  #  randomMean = exp(mean(log(randomlySelectedRatios)))
  randomMean <- ratiosMean(randomlySelectedRatios)
  return(randomMean)
}

DescriptiveStatistics <- function(distribution, robust = FALSE) {
  #	print(paste("Robust is ", robust))
  if (robust) {
    centralTendency <- median(distribution)
    variability <- mad(distribution, constant = 1) # median absolute deviation
  } else {
    centralTendency <- mean(distribution)
    variability <- sd(distribution)
  }
  return(list(centralTendency = centralTendency,
              variability = variability))
}

SampleNoiseGenes <- function(numAmpl = 2,
                             ratios,
                             replicates = 100,
                             probs = NULL,
                             significanceLevel = 0.05,
                             robust = FALSE) {
  margin <- significanceLevel / 2
  # now repeat the sampling for the selected number of
  # amplicons replicates
  sampleNoiseDistribution = replicate(replicates,
                                      SampleRatio(ratios,
                                                  numAmpl,
                                                  amplWeights = probs))
  #get the upper, mean and lower values of the noise distribution
  logSampleNoiseDistribution <- log(sampleNoiseDistribution)

  ds <- DescriptiveStatistics(logSampleNoiseDistribution, robust)
  logMeanNoise <- ds$centralTendency
  logSdNoise <- ds$variability

  # if (robust) {
  #     logMeanNoise <- median(logSampleNoiseDistribution)
  #     logSdNoise <- mad(logSampleNoiseDistribution, constant = 1)
  # } else {
  #     logMeanNoise <- mean(logSampleNoiseDistribution)
  #     logSdNoise <- sd(logSampleNoiseDistribution)
  # }

  lowerBound <- exp(logMeanNoise + qnorm(margin) * logSdNoise)
  meanNoise <- exp(logMeanNoise)
  upperBound <- exp(logMeanNoise + qnorm(1 - margin) * logSdNoise)

  sampledNoise = c(lowerBound, meanNoise, upperBound)
  names(sampledNoise) = paste0(c("Lower", "Mean", "Upper"), "Noise")
  return(sampledNoise)
}

countAmplicons <- function(x) {
  return(sapply(x, NROW))
}

# Returns a list with the unique number of amplicons for all genes
NumberOfUniqueAmplicons <- function(genesPos) {
  # calcUniqueAmpliconNumbers <-function(genesPos) {
  ampliconNumbers = countAmplicons(genesPos)
  uniqueAmpliconNumbers = sort(unique(ampliconNumbers))
  return(uniqueAmpliconNumbers)
}

#now calculate the background for each of the unique amplicon numbers
# ratios: the
IterateAmplNum <- function(uniqueAmpliconNumbers,
                           ratios,
                           replicates = 100,
                           probs = NULL,
                           significanceLevel = 0.05,
                           robust = FALSE) {
  # Needed because of package check complaints..
  i <- NULL
  noiseResults = foreach(i = seq_along(uniqueAmpliconNumbers)) %do% {
    sampledNoise = SampleNoiseGenes(uniqueAmpliconNumbers[i],
                                    ratios = ratios,
                                    replicates = replicates,
                                    probs = probs,
                                    significanceLevel = significanceLevel,
                                    robust = robust)
    sampledNoise
  }
  names(noiseResults) = as.character(uniqueAmpliconNumbers)
  return(noiseResults)
}

SelectReferenceSetByPercentil <- function(allSamplesReadCounts,
                                          normalizationMethod = "tmm",
                                          lowerBoundPercentage = 1,
                                          upperBoundPercentage = 99) {

  allSamplesReadCountsNormalized <-NormalizeCounts(allSamplesReadCounts,
                                                   method = normalizationMethod)

  selectedSamplesFilepath <- NULL
  allSamples <- round(allSamplesReadCountsNormalized)
  samplesColnames <- colnames(allSamples)

  selectedRange <- list()
  for (i in 1:nrow(allSamples)) {
    #    selectedRange[[i]] <- quantile(allSamples[i,])
    selectedRange[[i]] <- quantile(allSamples[i,], c(lowerBoundPercentage, upperBoundPercentage) / 100)
  }

  selectedSamplesIndex <- NULL
  excludedSamplesIndex <- NULL
  for (j in 1:ncol(allSamples)) {
    addSample <- TRUE
    for(i in 1:nrow(allSamples)) {
      if ((allSamples[i, j] < selectedRange[[i]][1]) | (allSamples[i, j] > selectedRange[[i]][2])) {
        addSample <- FALSE
        #        print(paste("Add sample?", addSample, i, selectedRange[[i]][1], allSamples[i, j], selectedRange[[i]][2]))
        break
      }
    }
    if (addSample) {
      selectedSamplesIndex <- cbind(selectedSamplesIndex, j)
    } else {
      excludedSamplesIndex <- cbind(excludedSamplesIndex, j)
    }
  }

  selectedSamples <- allSamples[, selectedSamplesIndex]
  excludedSamples <- allSamples[, -excludedSamplesIndex]
  selectedSamplesFilepath <- colnames(selectedSamples)
  return(selectedSamplesFilepath)
}

SelectReferenceSetByInterquartileRange <- function(allSamplesReadCounts, normalizationMethod = "tmm", iqrFactor = 1) {
  allSamplesReadCounts <- NormalizeCounts(allSamplesReadCounts, method = normalizationMethod)
  fivenumReads <- t(apply(allSamplesReadCounts, 1, fivenum))
  colnames(fivenumReads) <- c("SampleMinimum", "LowerQuartile", "Median", "UpperQuartile", "SampleMaximum")
  iqr <- fivenumReads[, "UpperQuartile"] - fivenumReads[, "LowerQuartile"]
  dispersion <- iqrFactor * iqr
  sampleIndexes <- NULL
  for (i in 1:ncol(allSamplesReadCounts)) {
    if (all((allSamplesReadCounts[, i] >= (fivenumReads[,"LowerQuartile"] - dispersion ))  & ((dispersion + fivenumReads[,"UpperQuartile"])  >= allSamplesReadCounts[, i]))) {
      sampleIndexes <- c(sampleIndexes, i)
    }
  }
  return(colnames(allSamplesReadCounts)[sampleIndexes])
}

SelectReferenceSetByKmeans <- function(allSamplesReadCounts, normalizationMethod = "tmm", referenceNumberOfElements) {
  # referenceSamples <- function(referenceNumberOfElements, allSamples) {

  allSamples <- NormalizeCounts(allSamplesReadCounts, method = normalizationMethod)

  #  allSamples <- allSamplesReadCounts
  numberOfSamples <- ncol(allSamples)
  numberOfOutliers <- numberOfSamples - referenceNumberOfElements
  if (numberOfOutliers < 0) {
    message(paste("Number of samples", numberOfSamples, "is smaller than the size of the reference set requested", referenceNumberOfElements, ", so all samples will be selected"))
    return(colnames(allSamples))
  } else if (numberOfOutliers == 0) {
    message(paste("Number of samples", numberOfSamples, "is equal to the size of the reference set requested", referenceNumberOfElements))
    return(colnames(allSamples))
  } else {
    allSamples <- t(allSamples)
    kmeans.result <- kmeans(allSamples, centers=1)
    # "centers" is a data frame with one center
    centers <- kmeans.result$centers[kmeans.result$cluster, ]
    distances <- sqrt(rowSums((allSamples - centers)^2))
    outliers <- order(distances, decreasing=TRUE)[1:numberOfOutliers]
    # these rows are top outliers
    message(cat("outliers:", colnames(allSamplesReadCounts)[outliers]))
    allSamples <- t(allSamples)
    #    return(colnames(allSamples[, -outliers, drop = FALSE]))
    return(colnames(allSamples)[-outliers])
  }
}

SelectReferenceSet <- function(samplesFilepaths,
                               genomicRanges,
                               removePcrDuplicates = FALSE,
                               normalizationMethod = "tmm",
                               referenceMaximumNumberOfElements = 10,
                               referenceSelectionMethod = "percentil",
                               lowerBoundPercentage = 1,
                               upperBoundPercentage = 99) {

  #  samplesFilepaths <- allSamplesFilteredFilepaths[1:3]
  #  genomicRanges <- genomicRangesFromBed

  allSamplesReadCounts <- ReadCountsFromBam(bamFilenames = samplesFilepaths,
                                            sampleNames = samplesFilepaths,
                                            gr = genomicRanges,
                                            #                                            ampliconNames = ampliconNames,
                                            removeDup = removePcrDuplicates)

  selectedSamplesFilepath <- SelectReferenceSetFromReadCounts(allSamplesReadCounts,
                                                              #                                                              genomicRanges,
                                                              #                                                              removePcrDuplicates = removePcrDuplicates,
                                                              normalizationMethod = normalizationMethod,
                                                              referenceMaximumNumberOfElements = referenceMaximumNumberOfElements,
                                                              referenceSelectionMethod = referenceSelectionMethod,
                                                              lowerBoundPercentage = lowerBoundPercentage,
                                                              upperBoundPercentage = lowerBoundPercentage)
  return(selectedSamplesFilepath)
}

SelectReferenceSetFromReadCounts <- function(allSamplesReadCounts,
                                             #                               genomicRanges,
                                             #                               removePcrDuplicates = FALSE,
                                             normalizationMethod = "tmm",
                                             referenceMaximumNumberOfElements = 30,
                                             #                               referenceSelectionMethod = "percentil",
                                             referenceSelectionMethod = "kmeans",
                                             lowerBoundPercentage = 1,
                                             upperBoundPercentage = 99) {


  selectedSamplesFilepath <- NULL
  if (referenceSelectionMethod == "percentil") {
    message("Selecting reference by percentil")
    selectedSamplesFilepath <- SelectReferenceSetByPercentil(allSamplesReadCounts,
                                                             normalizationMethod = normalizationMethod,
                                                             #                                                              referenceSelectionMethod = referenceSelectionMethod,
                                                             lowerBoundPercentage = lowerBoundPercentage,
                                                             upperBoundPercentage = upperBoundPercentage)
  } else if (referenceSelectionMethod == "kmeans") {
    message("Selecting reference by kmeans")
    selectedSamplesFilepath <- SelectReferenceSetByKmeans(allSamplesReadCounts,
                                                          normalizationMethod = normalizationMethod,
                                                          referenceMaximumNumberOfElements)
  } else {
    message(paste("method", referenceSelectionMethod, "not supported"))
  }
  return(selectedSamplesFilepath)
}

CollectColumnFromAllReportTables <- function(reportTables, columnName) {
  compilation <- NULL
  mySampleNames <- names(reportTables)
  for(name in mySampleNames) {
    tmpReportTable <- reportTables[[name]]
    compilation <- cbind(compilation, as.vector(tmpReportTable[, columnName]))
  }
  colnames(compilation) <- mySampleNames
  rownames(compilation) <- rownames(reportTables[[1]])
  return(compilation)
}

revalueDF <- function(myDataFrame, mappings) {
  myDataFrameColnames <- colnames(myDataFrame)
  for(columnName in myDataFrameColnames) {
    myDataFrame[, columnName] <- revalue(myDataFrame[, columnName], mappings)
  }
  mode(myDataFrame) <- "numeric"
  return(myDataFrame)
}

## TODO check this function
# ReliableAberrationStatusHeatMap <- function(meanBootResults,
#                                             realiableResults,
#                                             filepath = NULL) {
#   meanBootResults <- melanomaMeanBootResults
#   reliableResults <- melanomaCNVPanelizerResults
#   test_that("testing correct dimensions", {
#     expect_equal(colnames(meanBootResults), colnames(reliableResults))
#     expect_equal(rownames(meanBootResults), rownames(reliableResults))
#   })
#
#   kNormalDefaultValue <- 1.0
#   kAmplificationMaximalValue <- 4.0
#   kAmplificationMinimumValue <- 2.0
#
#   for(i in rownames(reliableResults)) {
#     for(j in colnames(reliableResults)) {
#       if (reliableResults[i,j] == "Normal") {
#         #        print(paste("Normal", i, reliableResults[i,j], "setting ", meanBootResults[i,j],"as", kNormalDefaultValue))
#         meanBootResults[i,j] <- kNormalDefaultValue
#       }
#
#       if ((reliableResults[i,j] == "Amplification") & (meanBootResults[i,j] > kAmplificationMaximalValue)) {
#         print(paste("Too high ", i, reliableResults[i,j], "setting ", meanBootResults[i,j],"as ", kAmplificationMaximalValue))
#         meanBootResults[i,j] <- kAmplificationMaximalValue
#       }
#
#       if ((reliableResults[i,j] == "Amplification") & (meanBootResults[i,j] < kAmplificationMinimumValue)) {
#         print(paste("Amplification too low", i, reliableResults[i,j], "setting ", meanBootResults[i,j],"as", kNormalDefaultValue))
#         meanBootResults[i,j] <- kNormalDefaultValue
#       }
#     }
#   }
#
#   df.team_data <- melt(meanBootResults)
#   colnames(df.team_data) <- c("Gene", "Sample", "Bootratio")
#   g <- ggplot(data = df.team_data,
#               aes(x = Gene, y = Sample)) +
#     #    geom_tile(aes(fill = scale(Bootratio))) +
#     geom_tile(aes(fill = log(Bootratio))) +
#     #      + opts(theme(axis.text.x = element_text(angle = 180, hjust = 1)))
#     #       scale_fill_gradient2(low="darkblue", high="darkgreen", guide="colorbar") +
#     #      scale_fill_gradient2(low="green4", high="red4", guide="colorbar") +
#     scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     theme(axis.text.y = element_text(size=5))
#   if (!missing(filepath)) {
#     ggsave(filepath)
#   }
# }



# why these colors by default:
StatusHeatmap <- function(dfData,
                          statusColors = c("Deletion" = "blue",
                                           "Normal" = "green",
                                           #                                         " " = "green",
                                           "Amplification" = "red"),
                          header = "Status Heatmap",
                          filepath = "CNVPanelizerHeatMap.png") {

  mappings <- seq(length(statusColors))
  #  mappings <- c(-1, 0, 1)
  names(mappings) <- names(statusColors)
  matrixData <- revalueDF(dfData, mappings)

  valuesUsed <- table(matrixData)
  indexOfColorsUsed <- as.numeric(names(valuesUsed))
  colorsRequired <- unname(statusColors[indexOfColorsUsed])

  # TODO find out why do I need this horrible hack when only one status is available
  if (length(colorsRequired) <= 1) {
    colorsRequired <- c("green", "yellow")
  }

  #  colorsRequired <- as.vector(statusColors[as.numeric(names(table(matrixData)))])

  png(filename = file.path(filepath),    # create PNG for the heat map
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size

  lmat <- rbind( c(5,3,4), c(2,1,4) )
  lhei <- c(1.5, 4)
  lhei <- c(1, 10)
  lwid <- c(1.5, 4, 0.75)

  myCexRow <- 0.2

  heatmap.2(matrixData,
            main = header,
            # plot layout
            lmat = lmat,
            lhei = lhei,
            lwid = lwid,
            #          scale = "row",
            #          ColSideColors=heatColors,
            key = FALSE, # whether a color-key should be shown
            #          key = TRUE, # whether a color-key should be shown
            #          key.par=list(match=c(0),
            #             amplification_deletion=c(1, 2, 3, 4),
            #             deletion_amplification=c(5, 6)),
            dendrogram = "none",
            tracecol=NA,
            #          density.info=c("histogram","density","none"),
            #            col = unname(statusColors),
            #            col = c("green"),
            col = colorsRequired,    # if the values are not all revalued the colors have to be passed selectively

            #          col = heatColors,
            #          col = heatColors[matchPoints],
            colsep=1:ncol(matrixData),
            rowsep=1:nrow(matrixData),
            margins = c(15, 5),
            Rowv=FALSE,
            Colv=FALSE,
            #            cexRow=myCexRow
  )

  #   legend("topleft", inset=.02, title= "UOUOU",
  #          "this is legend",
  #          #         fill=topo.colors(3),
  # #         fill=unname(statusColors),
  #          horiz=TRUE, cex=0.8)
  #bottomleft
  #  legend("bottom", inset=.01, title="CNV Status",





  #
  #
  #          legend("left", inset=.01, title="CNV Status",
  # #                legend("left", inset=0, title="CNV Status",
  #          names(statusColors),
  #          #         fill=topo.colors(3),
  #          fill=unname(statusColors),
  #          cex=0.8)
  # #         horiz=TRUE, cex=0.8)
  #




  #   legend("bottomleft", legend=c("Line 1", "Line 2"),
  #          col=c("red", "blue"), lty=1:2, cex=0.8)

  #   legend(1, 95, legend=c("Line 1", "Line 2"),
  #          col=c("red", "blue"), lty=1:2, cex=0.8)

  # in case of mistake this instruction should be called until return null instead of png..
  dev.off()
}

## TODO check
# StatusStability <- function(geneNames, sampleNormalizedReadCounts, tmpReferenceNormalizedReadCounts) {
#   centralTendency <-  matrix(rowMeans(tmpReferenceNormalizedReadCounts), ncol = 1)
#   rownames(centralTendency) <- rownames(tmpReferenceNormalizedReadCounts)
#   referenceAndSampleDifferences <- centralTendency - sampleNormalizedReadCounts
#   rownames(referenceAndSampleDifferences) <- rownames(tmpReferenceNormalizedReadCounts)
#   colnames(referenceAndSampleDifferences) <- "difference"
#   uniqueGeneNames <- unique(geneNames)
#   geneStabilityStatus <- c()
#   for (i in seq(uniqueGeneNames)) {
#     geneAmpliconIndexes <- grep(uniqueGeneNames[i], rownames(myReference))
#     geneAmpliconStatus <- referenceAndSampleDifferences[geneAmpliconIndexes, ]
#     allGeneAmpliconsHaveSameStatus <- length(unique(geneAmpliconStatus >= 0))==1
#     geneStabilityStatus <- c(geneStabilityStatus, allGeneAmpliconsHaveSameStatus)
#   }
#   names(geneStabilityStatus) <- uniqueGeneNames
#   return(geneStabilityStatus)
# }

## TODO CHECK
# StatusStabilityTable <- function(geneNames, tmpSamplesNormalizedReadCounts, tmpReferenceNormalizedReadCounts) {
#   statusStabilityTable <- NULL
#   for (i in seq(ncol(tmpSamplesNormalizedReadCounts))) {
#     statusStabilityTable <- cbind(statusStabilityTable, StatusStability(geneNames, matrix(tmpSamplesNormalizedReadCounts[, i], ncol = 1), tmpReferenceNormalizedReadCounts))
#   }
#   colnames(statusStabilityTable) <- colnames(tmpSamplesNormalizedReadCounts)
#   #   statusStabilityTable <- t(statusStabilityTable)
#   #   colnames(statusStabilityTable) <- colnames(tmpSamplesNormalizedReadCounts)
#   return(statusStabilityTable)
# }

PlotBootstrapDistributions  <- function(bootList,
                                        reportTables,
                                        outputFolder = getwd(),
                                        sampleNames = NULL,
                                        save = FALSE,
                                        #                                        scale = 7) {
                                        scale = 10) {
  selSample <- NULL
  plotList <- foreach(selSample = seq_along(bootList)) %do% {
    test <- as.factor(reportTables[[selSample]][, "Passed"])

    levelLabels <- c("NoChange", "NonReliableChange", "ReliableChange")
    names(levelLabels) <- c(0, 1, 2)
    test <- suppressMessages(revalue(test, levelLabels))
    namedColors <- c("#56B4E9", "#CC79A7" ,"#D55E00")
    names(namedColors) <- levelLabels

    ratios <- NULL
    testsPassed <- NULL
    df <- data.frame(class = as.factor(colnames(bootList[[selSample]])),
                     ratios = (as.vector(t(bootList[[selSample]]))),
                     testsPassed = test)

    #if a genomic ranges object has been
    #if(!is.null(gr)) {
    # newGeneOrder <- bedToGeneOrder(gr)
    # df$class = with(df,factor(class,levels(class)[newGeneOrder]))
    #}

    ylim1 <- boxplot.stats((df$ratios))$stats[c(1, 5)]
    ylim1 <- c(max(ylim1),max(ylim1))
    ylim1 <- log2(ylim1)
    ylim1 <- ylim1 * c(-scale,scale)

    if(is.null(sampleNames)) {
      filename <- names(bootList[selSample])
      filepath <- paste0(outputFolder, "/", filename, "_plot.pdf")
    } else {
      filename <- sampleNames[selSample]
      filepath <- paste0(outputFolder, "/", filename, ".pdf")
    }
    bootPlot <- ggplot(df, aes(x = class, y = log2(ratios), fill = testsPassed)) +
      geom_violin() + ggtitle(filename) +
      theme(plot.title = element_text(lineheight = 0.8, face = "bold"),
            text = element_text(size = 15),
            axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(name = "CNV Reliability",
                        values = namedColors,
                        labels = c("0" = "Foo", "1" = "Bar")) +
      coord_cartesian(ylim = ylim1) +
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=20,face="bold")) +
      scale_x_discrete("Gene Names") +
      #       geom_hline(yintercept=log2(1.5), color="#009E73") +
      #       geom_hline(yintercept=log2(0.5), color="#009E73") +
      geom_hline(yintercept=0, color="#009E73")

    if(save == TRUE) {
      message(paste0("Saving plot to '", filepath, "'"))
      #            dir.create(outputFolder, recursive = TRUE, showWarnings = TRUE)
      dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
      ggsave(filename = filepath,
             plot = bootPlot,
             height = 7.42 * 1,
             width = 8.11 * 2,
             limitsize = FALSE)
    }
    #        dev.off()
    return(bootPlot)
  }
  names(plotList) <- names(reportTables)
  return(plotList)
}

CNVPanelizerFromReadCountsHELPER <- function(sampleReadCounts,
                                             referenceReadCounts,
                                             genomicRangesFromBed,
                                             numberOfBootstrapReplicates = 10000,
                                             normalizationMethod = "tmm",
                                             robust = TRUE,
                                             backgroundSignificanceLevel = 0.05,
                                             outputDir = file.path(getwd(), "CNVPanelizer"),
                                             splitSize = 5) {
  myCNVPanelizerResults <- list()

  splits <- split(1:ncol(sampleReadCounts), ceiling(seq_along(1:ncol(sampleReadCounts))/splitSize))

  #    for (i in 1:ncol(sampleReadCounts)) {

  for (i in 1:length(splits)) {
    myCNVPanelizerResults[[i]] <- CNVPanelizerFromReadCounts(sampleReadCounts[, splits[[i]]],
                                                             referenceReadCounts,
                                                             genomicRangesFromBed,
                                                             numberOfBootstrapReplicates = numberOfBootstrapReplicates,
                                                             normalizationMethod = normalizationMethod,
                                                             robust = robust,
                                                             backgroundSignificanceLevel = backgroundSignificanceLevel,
                                                             outputDir = outputDir)
  }

  n <- length(myCNVPanelizerResults)
  res <- NULL
  for (i in seq(n)) {
    res <- cbind(res, myCNVPanelizerResults[[i]])
  }

  myCNVPanelizerResultsCompiled <- myCNVPanelizerResults[[1]]

  if(length(myCNVPanelizerResults) > 1) {
    for (i in 2:length(myCNVPanelizerResults)) {
      myCNVPanelizerResultsCompiled$sampleReadCounts <- cbind(myCNVPanelizerResultsCompiled$sampleReadCounts, myCNVPanelizerResults[[i]]$sampleReadCounts)
      myCNVPanelizerResultsCompiled$bootList <- append(myCNVPanelizerResultsCompiled$bootList, myCNVPanelizerResults[[i]]$bootList)
      myCNVPanelizerResultsCompiled$backgroundNoise <- append(myCNVPanelizerResultsCompiled$backgroundNoise, myCNVPanelizerResults[[i]]$backgroundNoise)
      myCNVPanelizerResultsCompiled$reportTables <- append(myCNVPanelizerResultsCompiled$reportTables, myCNVPanelizerResults[[i]]$reportTables)
    }
  }
  return(myCNVPanelizerResultsCompiled)
}

######################################################################################
# still on testing..
######################################################################################
# CNVPanelizerFromReadCountsHelperParallel <- function(sampleReadCounts,
#                                                      referenceReadCounts,
#                                                      genomicRangesFromBed,
#                                                      numberOfBootstrapReplicates = 10000,
#                                                      normalizationMethod = "tmm",
#                                                      robust = TRUE,
#                                                      backgroundSignificanceLevel = 0.05,
#                                                      outputDir = file.path(getwd(), "CNVPanelizer"),
#                                                      splitSize = 5) {
#
#   myCNVPanelizerResults <- list()
#
#   splits <- split(1:ncol(sampleReadCounts), ceiling(seq_along(1:ncol(sampleReadCounts))/splitSize))
#
#   # Use the detectCores() function to find the number of cores in system
#   no_cores <- detectCores()
#   # Setup cluster
#   clust <- makeCluster(no_cores - 1)
#   clusterExport(clust,
#                 varlist = c(
#                   "sampleReadCounts",
#                   "referenceReadCounts",
#                   "genomicRangesFromBed",
#                   "numberOfBootstrapReplicates",
#                   "normalizationMethod",
#                   "robust",
#                   "backgroundSignificanceLevel",
#                   "outputDir",
#                   "CNVPanelizerFromReadCounts"),
#                 envir = environment())
#
#   myCNVPanelizerResults <- parLapply(clust, 1:length(splits), function(i) {
#     CNVPanelizerFromReadCounts(sampleReadCounts[, splits[[i]]],
#                                referenceReadCounts,
#                                genomicRangesFromBed,
#                                numberOfBootstrapReplicates = numberOfBootstrapReplicates,
#                                normalizationMethod = normalizationMethod,
#                                robust = robust,
#                                backgroundSignificanceLevel = backgroundSignificanceLevel,
#                                outputDir = outputDir)
#   })
#   stopCluster(clust)
#
#   n <- length(myCNVPanelizerResults)
#   res <- NULL
#   for (i in seq(n)) {
#     res <- cbind(res, myCNVPanelizerResults[[i]])
#   }
#
#   myCNVPanelizerResultsCompiled <- myCNVPanelizerResults[[1]]
#
#   if(length(myCNVPanelizerResults) > 1) {
#     for (i in 2:length(myCNVPanelizerResults)) {
#       myCNVPanelizerResultsCompiled$sampleReadCounts <- cbind(myCNVPanelizerResultsCompiled$sampleReadCounts, myCNVPanelizerResults[[i]]$sampleReadCounts)
#       myCNVPanelizerResultsCompiled$bootList <- append(myCNVPanelizerResultsCompiled$bootList, myCNVPanelizerResults[[i]]$bootList)
#       myCNVPanelizerResultsCompiled$backgroundNoise <- append(myCNVPanelizerResultsCompiled$backgroundNoise, myCNVPanelizerResults[[i]]$backgroundNoise)
#       myCNVPanelizerResultsCompiled$reportTables <- append(myCNVPanelizerResultsCompiled$reportTables, myCNVPanelizerResults[[i]]$reportTables)
#     }
#   }
#   return(myCNVPanelizerResultsCompiled)
# }

CNVPanelizerFromReadCounts <- function(sampleReadCounts,
                                       referenceReadCounts,
                                       genomicRangesFromBed,
                                       numberOfBootstrapReplicates = 10000,
                                       normalizationMethod = "tmm",
                                       robust = TRUE,
                                       backgroundSignificanceLevel = 0.05,
                                       outputDir = file.path(getwd(), "CNVPanelizer")) {

  #  VerifiyIfOutputDirectoryExistsOrIsNotEmpty(outputDir)


  metadataFromGenomicRanges <- elementMetadata(genomicRangesFromBed)
  geneNames = metadataFromGenomicRanges["geneNames"][, 1]

  #
  # #   metadataFromGenomicRanges <- elementMetadata(genomicRangesFromBed)
  # #   geneNames = metadataFromGenomicRanges["geneNames"][, 1]
  # #   ampliconNames = metadataFromGenomicRanges["ampliconNames"][, 1]
  #
  #   sampleReadCounts <- ReadCountsFromBam(bamFilenames = sampleBamFilepaths,
  #                                         sampleNames = sampleBamFilepaths,
  #                                         gr = genomicRangesFromBed,
  #                                         removeDup = FALSE)
  #
  #   sampleReadCounts <- ReadCountsFromBam(bamFilenames = sampleBamFilepaths,
  #                                         sampleNames = sampleBamFilepaths,
  #                                         gr = genomicRangesFromBed,
  #                                         removeDup = FALSE)

  normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts,
                                                   referenceReadCounts,
                                                   method = normalizationMethod
                                                   #                                                   ,                                                   ampliconNames = ampliconNames
  )

  # After normalization data sets need to be splitted again to perform bootstrap
  samplesNormalizedReadCounts = normalizedReadCounts["samples"][[1]]
  referenceNormalizedReadCounts = normalizedReadCounts["reference"][[1]]

  bootList <- BootList(geneNames,
                       samplesNormalizedReadCounts,
                       referenceNormalizedReadCounts,
                       replicates = numberOfBootstrapReplicates)

  # Estimate the background noise left after normalization
  backgroundNoise <- Background(geneNames,
                                samplesNormalizedReadCounts,
                                referenceNormalizedReadCounts,
                                bootList,
                                replicates = numberOfBootstrapReplicates,
                                significanceLevel = backgroundSignificanceLevel,
                                robust = robust)

  # Build report tables
  reportTables <- ReportTables(geneNames,
                               samplesNormalizedReadCounts,
                               referenceNormalizedReadCounts,
                               bootList,
                               backgroundNoise)

  plots = PlotBootstrapDistributions(bootList,
                                     reportTables,
                                     sampleNames = names(reportTables),
                                     outputFolder = file.path(outputDir, "plots"),
                                     save = TRUE)

  WriteListToXLSX(reportTables,
                  multipleFiles = TRUE,
                  outputFolder = file.path(outputDir, "xlsx"))

  # WriteListToXLSX(reportTables, outputFolder =  "", filepath = "list.xlsx") {
  #   write.xlsx(listOfDataFrames, , rowNames = TRUE)
  # }
  #
  #
  #
  #
  # WriteListToXLSX(tmp1@reportTables, filepath = file.path(outputDir, "reportTables.xlsx"))
  # WriteListToXLSX(tmp1@reportTables, filepath = file.path(getwd(), "CNVPanelizer", "reportTables.xlsx"))

  # CNVPanelizerResults <- setClass("CNVPanelizerResults",
  #                                 #                                  setClass("CNVPanelizerResults",
  #
  #                                 # use this function to extract the slot names of the slotNames(results_CNVPanelizer[[1]])
  #                                 slots = c(sampleReadCounts="matrix",
  #                                           referenceReadColunts="matrix",
  #                                           genomicRangesFromBed="GRanges",
  #                                           bootList="list",
  #                                           backgroundNoise="list",
  #                                           plots = "list",
  #                                           reportTables="list"))
  #
  # myCNVPanelizerResults <- CNVPanelizerResults(sampleReadCounts = sampleReadCounts,
  #                                              referenceReadColunts = referenceReadCounts,
  #                                              genomicRangesFromBed = genomicRangesFromBed,
  #                                              bootList = bootList,
  #                                              backgroundNoise = backgroundNoise,
  #                                              plots = plots,
  #                                              reportTables = reportTables)

  myCNVPanelizerResults <- list(sampleReadCounts = sampleReadCounts,
                               referenceReadCounts = referenceReadCounts,
                               genomicRangesFromBed = genomicRangesFromBed,
                               bootList = bootList,
                               backgroundNoise = backgroundNoise,
                               plots = plots,
                               reportTables = reportTables)



  return(myCNVPanelizerResults)
}

CNVPanelizer <- function(sampleBamFilepaths,
                         referenceBamFilepaths,
                         bedFilepath,
                         amplColumnNumber = 6,
                         minimumMappingQuality = 20,
                         numberOfBootstrapReplicates = 10000,
                         removePcrDuplicates = TRUE,
                         #                         analysisMode = "gene",   # analysisMode can be "gene" or "amplicon"
                         robust = TRUE,
                         backgroundSignificanceLevel = 0.05,
                         outputDir = file.path(getwd(), "CNVPanelizer")) {

  #  VerifiyIfOutputDirectoryExistsOrIsNotEmpty(outputDir)

  genomicRangesFromBed <- BedToGenomicRanges(bedFilepath,
                                             ampliconColumn = amplColumnNumber,
                                             split = "_")

  #  metadataFromGenomicRanges <- elementMetadata(genomicRangesFromBed)
  #  geneNames = metadataFromGenomicRanges["geneNames"][, 1]
  #  ampliconNames = metadataFromGenomicRanges["ampliconNames"][, 1]

  sampleReadCounts <- ReadCountsFromBam(bamFilenames = sampleBamFilepaths,
                                        sampleNames = sampleBamFilepaths,
                                        gr = genomicRangesFromBed,
                                        minimumMappingQuality = minimumMappingQuality,
                                        removeDup = removePcrDuplicates)

  referenceReadCounts <- ReadCountsFromBam(bamFilenames = referenceBamFilepaths,
                                           sampleNames = referenceBamFilepaths,
                                           gr = genomicRangesFromBed,
                                           minimumMappingQuality = minimumMappingQuality,
                                           removeDup = removePcrDuplicates)

  results <- CNVPanelizerFromReadCounts(sampleReadCounts = sampleReadCounts,
                                        referenceReadCounts = referenceReadCounts,
                                        genomicRangesFromBed = genomicRangesFromBed,
                                        #                                        bedFilepath = bedFilepath,
                                        #                                        amplColumnNumber = amplColumnNumber,
                                        #                                       minimumMappingQuality = minimumMappingQuality,
                                        numberOfBootstrapReplicates = numberOfBootstrapReplicates,
                                        #                                       removePcrDuplicates = removePcrDuplicates,
                                        #                                       analysisMode = analysisMode,
                                        robust = robust,
                                        backgroundSignificanceLevel = backgroundSignificanceLevel,
                                        outputDir = outputDir)
  return(results)
}


Strict <- function(cnvPanelizerResults,
                   amplificationMinimumThreshold = 2,
                   deletionMaximumThreshold = 0.6,
                   minimumNumberOfAmplicons = 2
                   #                   , fullConsensusAmongAmplicons = TRUE
) {
#  myCNVPanelizerTableResults <- CollectColumnFromAllReportTables(cnvPanelizerResults@reportTables, "ReliableStatus")
  myCNVPanelizerTableResults <- CollectColumnFromAllReportTables(cnvPanelizerResults$reportTables, "ReliableStatus")
  myCNVPanelizerTableResults <- t(myCNVPanelizerTableResults)
  myCNVPanelizerTableResults <- revalueDF(myCNVPanelizerTableResults, c("Normal"=""))

  interestingColumn <- "MeanBoot"
#  meanBootResults <- CollectColumnFromAllReportTables(cnvPanelizerResults@reportTables, interestingColumn)
  meanBootResults <- CollectColumnFromAllReportTables(cnvPanelizerResults$reportTables, interestingColumn)
  class(meanBootResults) <- "numeric"
  meanBootResults <- t(meanBootResults)

#  consensusResults <- ConsensusCheck(cnvPanelizerResults@genomicRangesFromBed, cnvPanelizerResults@sampleReadCounts, cnvPanelizerResults@referenceReadCounts)
  consensusResults <- ConsensusCheck(cnvPanelizerResults$genomicRangesFromBed, cnvPanelizerResults$sampleReadCounts, cnvPanelizerResults$referenceReadCounts)
  consensusResults <- t(consensusResults)

  test_that("colnames from CNVPanelizer results and Consensus match", {
    expect_equal(colnames(consensusResults), colnames(myCNVPanelizerTableResults))
  })

  test_that("rownames from CNVPanelizer results and Consensus match", {
    expect_equal(rownames(consensusResults), rownames(myCNVPanelizerTableResults))
  })

  cNVPanelizerResultsConsensus <- myCNVPanelizerTableResults
  for (rowname in rownames(myCNVPanelizerTableResults)) {
    for (colname in colnames(myCNVPanelizerTableResults)) {
      #        print(paste(rowname, colname))
      if (!consensusResults[rowname, colname]) {
        cNVPanelizerResultsConsensus[rowname, colname] <- ""
        print(paste("Setting ", myCNVPanelizerTableResults[rowname, colname], "to", cNVPanelizerResultsConsensus[rowname, colname]))
      }
    }
  }

  #  table(unlist(myCNVPanelizerTableResults))
  #  table(unlist(cNVPanelizerResultsConsensus))

  test_that("colnames from CNVPanelizer results and meanBootResults match", {
    expect_equal(colnames(meanBootResults), colnames(cNVPanelizerResultsConsensus))
  })

  test_that("rownames from CNVPanelizer results and meanBootResults match", {
    expect_equal(rownames(meanBootResults), rownames(cNVPanelizerResultsConsensus))
  })

  cNVPanelizerResultsConsensusAndThreshold <- cNVPanelizerResultsConsensus
  for (rowname in rownames(cNVPanelizerResultsConsensusAndThreshold)) {
    for (colname in colnames(cNVPanelizerResultsConsensusAndThreshold)) {
      #      print(cNVPanelizerResultsConsensus[rowname, colname])}}
      if ((cNVPanelizerResultsConsensusAndThreshold[rowname, colname] == "Amplification")
          & ((meanBootResults[rowname, colname] <= amplificationMinimumThreshold))
          #          & ((meanBootResults[rowname, colname] <= 3))
      ) {
        cNVPanelizerResultsConsensusAndThreshold[rowname, colname] <- ""
        print(paste("Setting ", cNVPanelizerResultsConsensus[rowname, colname], "to", cNVPanelizerResultsConsensusAndThreshold[rowname, colname], "because", meanBootResults[rowname, colname]))
      }

      if ((cNVPanelizerResultsConsensusAndThreshold[rowname, colname] == "Deletion")
          & ((meanBootResults[rowname, colname] >= deletionMaximumThreshold))
          #& ((meanBootResults[rowname, colname] >= 0.7))
      ) {
        cNVPanelizerResultsConsensusAndThreshold[rowname, colname] <- ""
        print(paste("Setting ", cNVPanelizerResultsConsensus[rowname, colname], "to", cNVPanelizerResultsConsensusAndThreshold[rowname, colname]))
      }
    }
  }

  #  table(unlist(cNVPanelizerResultsConsensus))
  #  table(unlist(cNVPanelizerResultsConsensusAndThreshold))

  cNVPanelizerResultsConsensusAndThresholdAndMinimumNumberOfAmplicons <- cNVPanelizerResultsConsensusAndThreshold
#  genesWithMinimumNumberOfAmplicons <- HaveMininumNumberOfAmplicons(cnvPanelizerResults@genomicRangesFromBed, minimumNumberOfAmplicons)
  genesWithMinimumNumberOfAmplicons <- HaveMininumNumberOfAmplicons(cnvPanelizerResults$genomicRangesFromBed, minimumNumberOfAmplicons)
  for (rowname in rownames(cNVPanelizerResultsConsensusAndThresholdAndMinimumNumberOfAmplicons)) {
    for (colname in names(genesWithMinimumNumberOfAmplicons)) {
      if ((cNVPanelizerResultsConsensusAndThreshold[rowname, colname] == "Deletion") &
          !genesWithMinimumNumberOfAmplicons[colname]) {
        cNVPanelizerResultsConsensusAndThresholdAndMinimumNumberOfAmplicons[rowname, colname] <- ""
        print(paste("Setting ", cNVPanelizerResultsConsensusAndThreshold[rowname, colname], "to", cNVPanelizerResultsConsensusAndThresholdAndMinimumNumberOfAmplicons[rowname, colname]))
      }
    }
  }

  #  table(unlist(cNVPanelizerResultsConsensusAndThreshold))
  #  table(unlist(cNVPanelizerResultsConsensusAndThresholdAndMinimumNumberOfAmplicons))

  #  write.xlsx(cNVPanelizerResults, file=file.path(outputDir, paste0("cNVPanelizerResults", ".xlsx")), row.names = TRUE)
  #  write.xlsx(cNVPanelizerResultsConsensus, file=file.path(outputDir, paste0("cNVPanelizerResultsConsensus", ".xlsx")), row.names = TRUE)
  #  write.xlsx(cNVPanelizerResultsConsensusAndThreshold, file=file.path(outputDir, paste0("cNVPanelizerResultsConsensusAndThreshold", ".xlsx")), row.names = TRUE)
  #  write.xlsx(cNVPanelizerResultsConsensusAndThresholdAndMinimumNumberOfAmplicons, file=file.path(outputDir, paste0("cNVPanelizerResultsConsensusAndThresholdAndMinimumNumberOfAmplicons", ".xlsx")), row.names = TRUE)

  return(cNVPanelizerResultsConsensusAndThresholdAndMinimumNumberOfAmplicons)
}


# TODO Could be provided as a measurement of quality for the run
# Overview <- function(
#   #  sampleIdentifier,
#   sampleReadCounts1,
#   referenceReadCounts1,
#   allSamplesReadCounts,
#   genomicRangesFromBed,
#   #  normalizationMethod = "tss",
#   normalizationMethod = "tmm",
#   outputDir = file.path(getwd(), "CNVPanelizer", "overviewKNNReference30")) {
#
#   #  sampleIdentifier <- "I:/Run Data/Run S5-56/M1_S5-56_IonXpress_074_rawlib.bam"
#
#   #   Overview(allSamplesReadCounts[, filepaths],
#   #            allSamplesReadCounts[, referenceReadCounts],
#   #            #         allSamplesReadCounts,
#   #            genomicRangesFromBed,
#   #            #                     normalizationMethod = "tss",
#   #            #                     normalizationMethod = "tmm",
#   #            outputDir = file.path(getwd(), "CNVPanelizer", "overview"))
#
#   sampleReadCounts1 <- allSamplesReadCounts[, filepaths]
#   #  referenceReadCounts1 <- allSamplesReadCounts[, referenceReadCounts]
#   referenceReadCounts1 <- allSamplesReadCounts[, referenceFilepathsKmeans]
#
#   allSamplesReadCountsNormalized <- NormalizeCounts(allSamplesReadCounts, method = normalizationMethod)
#   samplesNormalizedReadCounts = allSamplesReadCountsNormalized[, colnames(sampleReadCounts1)]
#   referenceNormalizedReadCounts = allSamplesReadCountsNormalized[, colnames(referenceReadCounts1)]
#
#   runSamples <- samplesNormalizedReadCounts
#   refSamples <- referenceNormalizedReadCounts
#   entireSetOfSamples <- allSamplesReadCountsNormalized
#
#   geneNames <- unique(genomicRangesFromBed$geneNames)
#   ampliconNames <- genomicRangesFromBed$ampliconNames
#   # TODO validate if ampliconNames are the rownames of samples matrices
#   # ampliconNames <- rownames(allSamplesNormalized)
#
#   sampleIdentifiers <- colnames(runSamples)
#   #sampleIdentifiers <- sampleIdentifiers[1]
#   for (sampleIdentifier in sampleIdentifiers) {
#     myDir <- file.path(outputDir, basename(sampleIdentifier))
#     dir.create(myDir, recursive = TRUE)
#     for (geneName in geneNames) {
#       #    geneName <- "BRAF"
#       myAmpliconsIndexes <- grep(geneName, ampliconNames)
#       png(file=file.path(myDir, paste0(geneName,'_amplicons.png')),
#           width=1000,
#           height=(ceiling(length(myAmpliconsIndexes)/2) * 500))
#
#       #    par(mfrow=c(ceiling(length(myAmpliconsIndexes)/2), 2), mar=c(3,3,5,1))
#       par(mfrow=c(ceiling(length(myAmpliconsIndexes)/2), 2), mar=c(3, 3, 5, 1), oma=c(0,0,7,0))
#
#       for (ampliconName in ampliconNames[myAmpliconsIndexes]) {
#         #  ampliconName <- "chr7_140453102_140453221_BRAF_Ex15-NM_004333"
#         tmp <- data.frame(y = c(as.vector(runSamples[ampliconName, ])
#                                 ,as.vector(refSamples[ampliconName, ])
#                                 ,as.vector(entireSetOfSamples[ampliconName, ])
#         ),
#
#         x = factor(
#           c(rep("runSamples", length(as.vector(runSamples[ampliconName, ])))
#             ,rep("refSamples",length(as.vector(refSamples[ampliconName, ])))
#             ,rep("allSamples",length(as.vector(entireSetOfSamples[ampliconName, ])))
#           ), levels = c("runSamples", "refSamples", "allSamples"))
#
#         #                         x = c(rep("runSamples", length(as.vector(runSamples[ampliconName, ])))
#         #                               ,rep("refSamples",length(as.vector(refSamples[ampliconName, ])))
#         #                               )
#         )
#
#         myValues <- runSamples[ampliconName, sampleIdentifier]
#
#         boxplot(y~x,
#                 data=tmp,
#                 ylab="Read Counts",
#                 xlab="Sample set",
#                 names = c("Run Samples", "Reference", "All Samples"),
#                 #              main = ampliconName)
#                 main = paste(myValues, ampliconName))
#         stripchart(y ~ x,
#                    vertical = TRUE,
#                    data = tmp,
#                    method = "jitter",
#                    add = TRUE,
#                    #                pch = 20,
#                    pch = 20,
#                    col = 'blue')
#
#
#         #    stripchart(list(0.5, 1), vertical=TRUE, add=TRUE, method="stack", col='red', pch="*")
#         stripchart(myValues, vertical=TRUE, add=TRUE, method="stack", col='red', pch=16, cex=2)
#       }
#       title(main = paste(geneName, "-", sampleIdentifier), font.main = 4, outer=TRUE)
#       dev.off()
#     }
#   }
# }

## TODO CHECK THIS
#
#
# #rowIndexesToMerge <- IndexGenesPositions(lungGenomicRanges$geneNames)
# #GeneMultilinePlot(normalSamples, rowIndexesToMerge, normalizationMethod = "tss")
# GeneMultilinePlot <- function(dataSamples, rowIndexesToMerge, normalizationMethod =  "tmm") {
#   dataSamples <- round(NormalizeCounts(dataSamples))
#   result <- NULL
#   for (i in 1:ncol(dataSamples)) {
#     result <- cbind(result, sapply( rowIndexesToMerge, function(x) {mean(dataSamples[x, i])} ))
#   }
#   colnames(result) <- colnames(dataSamples)
#   dataSamples <- result
#   #  dataSamples <- cbind(geneReadCounts = seq_along(rownames(dataSamples)), dataSamples)
#
#   melted <- melt(dataSamples)
#
#   ggplot(data=melted, aes(x=Var1, y=value, group=Var2, colour=Var2, linetype = Var2)) +
#     #geom_line() +
#     #    scale_colour_manual(values = c('pink','orange','white', 'red')) +
#     theme(legend.position="left") +
#     geom_point(aes(color=Var2)) +
#     geom_line(size=1) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
# }
#

# #Example for testing MultiLinePlot
# dataSamples <- read.table(text = "CN CNVPanelizer panelcn.mops ioncopy
#                  CN1 0.1 0.4 0.4
#                  CN2 0.2 0.1 0.3
#                  CN3 0.1 0.2 0.1
#                  CN4 0.1 0.4 0.4
#                  CN5 0.2 0.1 0.3
#                  CN6 0.1 0.3 0.5
#                  CN7 0.6 0.6 0.6
#                  CN8 0.6 0.3 0.2
#                  CN9 0.6 0.8 0.8
#                  CN10 1 0.9 0.8", header=TRUE)

# TODO reuse this code at GeneMultilinePlot ...
# MultiLinePlot <- function(dataSamples, idColumnIdentifier = NULL, title = "", legendPosition = "left", xLabel = "", yLabel = "", smooth = FALSE) {
#   rs <- dataSamples
#   if (is.null(idColumnIdentifier)) {
#     melted = melt(rs)
#   } else {
#     melted = melt(rs, id.vars=idColumnIdentifier)
#   }
#   melted$CN <- as.character(melted$CN)
#   melted$CN <-factor(melted$CN, levels = unique(melted$CN))
#   myPlot <- ggplot(data=melted, aes(x=CN, y=value, group=variable, colour=variable, linetype = variable)) +
#     scale_colour_discrete(name = "Method") + #  scale_colour_manual(values = c('pink','orange','white', 'red')) +
#     ggtitle(title) +
#     geom_point(aes(color=variable)) + # plots the dots..
#     theme(text = element_text(size=10),
#           axis.text.x = element_text(angle=90, hjust=1)) +
#     theme(panel.border = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black")) +
#     theme_bw() + # sets the background colour
#     theme(plot.title = element_text(hjust = 0.5)) +
#     theme(legend.position=legendPosition) + xlab(xLabel) + ylab(yLabel)
#   if(smooth) {
#     #      myPlot = myPlot + geom_smooth(method = "loess", se = FALSE, span = 1.5, show.legend = FALSE)
#     myPlot = myPlot + geom_smooth(method = "auto", se = FALSE, span = 1.5, show.legend = FALSE)
#   } else {
#     myPlot = myPlot + geom_line(size = 1, show.legend = FALSE)    # links the dots with lines..
#   }
#   return(myPlot)
# }


SampleReadCountsPools <- function(ampliconsReadCounts, numberOfPools) {
  listOfPools <- list()
  for (i in 1:numberOfPools) {
    listOfPools[[i]] <- ampliconsReadCounts[seq(from = i, to = length(ampliconsReadCounts), by = numberOfPools)]
  }
  return(listOfPools)
}


SampleReadCountsByPool <-function(ampliconReadCounts, pools) {
  # existingPools <- as.numeric(names(table(pools)))
  # if (any(is.na(existingPools)) {
  #   stop("Pools are not defined as 1, 2, ...")
  # }
  existingPools <- names(table(pools))
  if (length(pools) != length(ampliconReadCounts)) {
    stop("Length  of pools is different than ampliconReadCounts")
  }
  listOfPools <- list()
  for (i in existingPools) {
    listOfPools[[i]] <- (ampliconReadCounts[which(pools==i)])
  }
  return(listOfPools)
}


PoolsPlots <- function(sampleReadCountsByPool) {
  plotValues <- as.integer(unlist(lapply(sampleReadCountsByPool, mean)))
  barplot(plotValues, ylab = "Read Counts mean", xlab="Pools", names.arg = names(sampleReadCountsByPool), col = "dodgerblue4",
          #          ylim=c(0, trunc(max(plotValues)))
  )
}








################################################################################
# KaryotypeAberrationPlot
################################################################################
# KaryotypeAberrationPlot(bedFilepath = bedFilepath,
#                         ampliconColumn = 4,
#                         filepath = "TESTANDO_APAGAR_KARYOTYPE.png")

KaryotypeAberrationPlot <- function(bedFilepath,
                                    ampliconColumn,
                                    status = sample(c("red", "green", "blue"), 30, replace = TRUE),
                                    filepath = "KaryotypeAberrationPlot.png",
                                    width = 1200,
                                    height = 800) {
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("chromPlot")
  #library(chromPlot)
  #library(stringr)
  hg_gap <- NULL # not sure if this is good solution to avoid 'NOTE'
  data(hg_gap)
  getGeneRegionsStatusSummary <- GetGeneRegionsStatusSummary(bedFilepath = bedFilepath,
                                                             ampliconColumn = ampliconColumn,
                                                             status = sample(c("red", "green", "blue"), 30, replace=TRUE))
  geneRegionsBands <- GetGeneRegionsBand(getGeneRegionsStatusSummary)
  png(filepath, width = width, height = height)
  #  chromPlot(gaps=hg_gap)
  #  chromPlot(gaps=hg_gap, stat=geneRegions, statCol="Value")
  #  chromPlot(gaps=hg_gap, stat=geneRegions, statCol="Value", chr=unique(geneRegions$Chrom))
  #  chromPlot(gaps=hg_gap, bands=hg_cytoBandIdeo, stat=geneRegions, statCol="Value", chr=unique(geneRegions$Chrom), statName="Value", noHist=TRUE, figCols=5, cex=0.7, statTyp="n", chrSide=c(1,1,1,1,1,1,-1,1))

  # TODO TEMPORARARIO PORQUE O chromPlot nao pode ser instalado em versao 3.2.0 ...
  # karyotypeAberrationPlot <- chromPlot(gaps = hg_gap,
  #                                      bands = geneRegionsBands,
  #                                      chr = unique(getGeneRegionsStatusSummary$Chrom),
  #                                      stat = getGeneRegionsStatusSummary,
  #                                      statCol = "Value",
  #                                      statName = "Value",
  #                                      noHist = TRUE,
  #                                      figCols = 7,
  #                                      cex = 0.7,
  #                                      statTyp = "n",
  #                                      chrSide = c(1,1,1,1,1,1,-1,1))
  karyotypeAberrationPlot <- NA
  dev.off()
  return(karyotypeAberrationPlot)
}


# It Merges all regions related to the same gene and assigns it a color that relates to the aberratiosn status
# Status should be a vector of colors representing the amplifications (red for amplification, green for deletion and blue for normal)
GetGeneRegionsStatusSummary <- function(bedFilepath, ampliconColumn = 4, split = "_", status = c("blue")) {
  a <- read.csv(bedFilepath, sep="\t", skip=1, header = FALSE)
  a[, 1] <- str_replace(a$V1, "chr","")
  a[, ampliconColumn] <- vapply(strsplit(as.vector(a[, ampliconColumn]), "_", fixed = TRUE), "[", "", 1)
  a <- a[,1:4]
  colnames(a) <- c("Chrom", "Start", "End", "ID"
                   #                   , "Colors"
  )
  geneNames <- unique(a$ID)
  geneRegions <- NULL
  for (geneName in geneNames) {
    genomicRegions <- a[a$ID==geneName, ]
    chromossome <- unique(genomicRegions$Chrom)
    start <- min(genomicRegions$Start)
    end <- max(genomicRegions$End)
    gene <- unique(genomicRegions$ID)
    geneRegions <- rbind(geneRegions, c(chromossome, start, end, gene))
  }
  colnames(geneRegions) <- colnames(a)
  geneRegions <- as.data.frame(geneRegions,
                               #                               as.is = TRUE)
                               stringsAsFactors = FALSE)
  geneRegions <- cbind(geneRegions, Colors = status)

  #  geneRegions$Chrom <- as.numeric(geneRegions$Chrom)
  geneRegions$Start <- as.numeric(geneRegions$Start)
  geneRegions$End <- as.numeric(geneRegions$End)
  geneRegions$ID <- as.factor(geneRegions$ID)
  return(geneRegions)
}


# Returns objected needed for plotting.. it changes some column names and makes the regions larger to be able to be visible in the plot..
#GetGeneRegionsBand <- function(geneRegionsStatusSummary) {
GetGeneRegionsBand <- function(geneRegions) {
  geneRegionsBands <- geneRegions
  colnames(geneRegionsBands) <- c("Chrom", "Start", "End", "Name", "Colors")
  geneRegionsBands$Chrom <- paste0("chr", geneRegionsBands$Chrom)
  geneRegionsBands$Name <- as.character(geneRegionsBands$Name)
  geneRegionsBands$Colors <- as.character(geneRegionsBands$Colors)
  # since the regions are too small to be visualized.. around 1 milion times smaller than the size of chromossome..
  geneRegionsBands$End <- geneRegionsBands$End + 1E6
  return(geneRegionsBands)
}






# TODO
# Could take in consideration the size of regions and
# the whole coverage of the gene (number of amplicons * seqLength)
# > 50 bp
# panelStatistics(genomicRangesFromBed) {
#   # hist(width(ranges(genomicRangesFromBed)))
#   # genomicRangesFromBed$geneNames
#
# }



#   boxplotGeneLevel <- function(allSamplesNormalized,
#                                cohortFilepaths,
#                                referenceFilepathsKNN,
#                                referenceFilepathsPercentil,
#                                geneName,
#                                outputDir = file.path(getwd(),"batchAnalysisNoMininumMappingQuality2")) {
#     #boxplotPerAmplicon <- function(ampliconIndex) {
#     #geneName <- "BRAF"
#     #geneName <- "APC"
#     #  geneName <- "ERBB4"
#
# #    setwd(outputDir)
#
#     #  cohortFilepaths <- filepaths
#     martina <- allSamplesNormalized[, (colnames(allSamplesNormalized) %in% cohortFilepaths)]
#
#     nonMartina <- allSamplesNormalized[, !(colnames(allSamplesNormalized) %in% cohortFilepaths)]
#
#
#     for (geneName in geneNames) {
#
#
#
#       ampliconNames <- rownames(allSamplesNormalized)
#       myAmpliconsIndexes <- grep(geneName, ampliconNames)
#
#       #  pdf(file=paste0(geneName,'_plot.pdf'), width=10, height= ceiling(length(myAmpliconsIndexes)/2) * 5)
#       png(file=paste0(geneName,'_plot.png'), width=10, height= ceiling(length(myAmpliconsIndexes)/2) * 5)
#
#       par(mfrow=c(ceiling(length(myAmpliconsIndexes)/2), 2), mar=c(3,3,5,1))
#
#       for (ampliconName in ampliconNames[myAmpliconsIndexes]) {
#         #  ampliconName <- "chr7_140453102_140453221_BRAF_Ex15-NM_004333"
#         tmp <- data.frame(y = c(as.vector(nonMartina[ampliconName, ]),
#                                 as.vector(nonMartina[ampliconName, referenceFilepathsKNN]),
#                                 as.vector(nonMartina[ampliconName, referenceFilepathsPercentil]),
#                                 as.vector(martina[ampliconName, ])),
#                           x = c(rep("nonMartina",length(as.vector(nonMartina[ampliconName, ]))),
#                                 rep("referenceKNN", length(as.vector(nonMartina[ampliconName, referenceFilepathsKNN]))),
#                                 rep("referencePercentil", length(as.vector(nonMartina[ampliconName, referenceFilepathsPercentil]))),
#                                 rep("martina", length(as.vector(martina[ampliconName, ])))))
#
#         boxplot(y~x, data=tmp, ylab="Read Counts", xlab="Sample set", main = ampliconName)
#         stripchart(y ~ x,
#                    vertical = TRUE,
#                    data = tmp,
#                    method = "jitter",
#                    add = TRUE,
#                    pch = 20,
#                    col = 'blue')
#       }
#       dev.off()
#     }
#   }





# i <- 0
# for (x in 1:100) {
#   i <- i+1
#   if (i %% 10 == 0) {
#     Sys.sleep(1)
#     cat(".")
#   }
# }
#
