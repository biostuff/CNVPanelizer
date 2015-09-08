###############################################################################
# functions to read the read counts
###############################################################################

BedToGenomicRanges <- function(panelBedFilepath,
                               ampliconColumn,
                               split = ";",
                               doReduce = TRUE,
                               rangeExtend = 0,
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
  return(gr)
}


#load the files you want to analyze
ReadCountsFromBam <- function(bamFilenames,
                              sampleNames,
                              gr,
                              ampliconNames,
                              removeDup = FALSE) {

  if (missing(ampliconNames)) {
    ampliconNames = elementMetadata(gr)$ampliconNames
  }

  # Because of package check complaints
  i <- NULL
  curbam = foreach(i = seq_along(bamFilenames), .combine = cbind) %do% {
    countBamInGRanges(bamFilenames[i],
                      gr,
                      remove.dup = removeDup,
                      min.mapq = 20,
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
                                      ampliconNames = NULL) {
  allCounts <- cbind(sampleCounts, referenceCounts)
  classes <- rep(c("samples", "reference"),
                 c(ncol(sampleCounts), ncol(referenceCounts)))

  tinyValueToAvoidZeroReadCounts <- 1e-05

  bamDataRangesNorm <- tmm(allCounts,
                           refColumn = NULL,
                           k = tinyValueToAvoidZeroReadCounts)

  if (!is.null(ampliconNames)) {
    rownames(bamDataRangesNorm) = ampliconNames
  }

  normalizedSamples <- bamDataRangesNorm[, which(classes == "samples"),
                                         drop = FALSE]

  normalizedReference <- bamDataRangesNorm[, which(classes == "reference"),
                                           drop = FALSE]

  return(list(samples = normalizedSamples, reference = normalizedReference))
}

###############################################################################
# helper functions
###############################################################################

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
PercentBelowValue <- function(vector, value) {
  sum(vector < value)/length(vector)
}

# check this params with Thomas
#index the bam files if there is no index yet
IndexMultipleBams <- function(bams,
                              index_type = ".bam.bai") {
  #check if the index already exists and need to be indexed
  potentialBaiFilenames <- gsub(".bam",
                                bams,
                                replacement = index_type,
                                ignore.case = TRUE)
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

WriteListToXLSX <- function(listOfDataFrames, filepath = "list.xlsx") {
  sizeOfList <- length(listOfDataFrames)
  dataFrameNames <- names(listOfDataFrames)
  index <- 1
  write.xlsx(listOfDataFrames[[index]],
             filepath,
             sheetName = dataFrameNames[index])
  for(i in seq(sizeOfList - 1)) {
    index <- i + 1
    write.xlsx(listOfDataFrames[[index]],
               filepath,               
               sheetName = dataFrameNames[index],
               append = TRUE)
  }
}

ReadXLSXToList <- function(filepath) {
  return(list(reference = unfactorize(read.xlsx(filepath,
                                                header = FALSE,
                                                sheetName = "reference")),
              sample = unfactorize(read.xlsx(filepath,
                                             header = FALSE,
                                             sheetName = "sample"))))
}

unfactorize <- function(df) {
  for(i in which(sapply(df, class) == "factor")) {
    df[[i]]=as.character(df[[i]])
  }
  return(df)
}

###############################################################################
# functions for bootstrapping
###############################################################################

#a function that randomly samples positions from a vector using subsampling
SubsamplingPositions <- function(pos, mtry = round(sqrt(length(pos)))) {
  return(sample(pos, mtry, replace = FALSE))
}

#calculates the mean from a vector given specific positions
PositionMean <- function(positon, vector) {
  return(mean(vector[positon]))
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
  return(ratioList)
}

ReferenceWeights <- function(refmat, varianceFunction = sd) {
  variance = apply(refmat, 2, varianceFunction)
  weights = 1/variance
  return(weights)
}

# split into a function generating the bootstrap distribution 
# and a function averaging over amplicons.
# define the function to generate bootstrap distributions and
# return a list with the genewise bootstrapping
BootList <- function(geneNames, sampleMatrix, refmat, replicates) {
  # get the genes positions in the matrix as a list from a gene name vector
  genesPos <- IndexGenesPositions(geneNames)

  # a vector and matrix are not the same and for a vector iterating over
  # the column makes no sense so we have
  # to check if a mtrix or a vector was passed. ncol only works for matrix
  # not for vector
  if (class(sampleMatrix) == "matrix") {
    iterator <- 1:ncol(sampleMatrix)
  } else {
    iterator <- 1
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
      } else {
        testSample <- sampleMatrix
      }

      # for each gene subsample the amplicon positions independently

      # sample the samples using bootstrapping
      sampleBootPos <- sample(1:ncol(refmat), ncol(refmat), replace = TRUE)

      geneBootPos <- c(lapply(genesPos, SubsamplingPositions),
                       recursive = TRUE)

      # given the obtained sampling using the bootstraps calulated above
      refMatPos <- refmat[geneBootPos, sampleBootPos, drop=FALSE]

      bootRatio <- testSample[geneBootPos]/rowMeans(refMatPos)

      # after the bootstrapping the gene positions in the vector changes
      # so recalculate them
      adjustedLength <- function(a) {
        round(sqrt(length(a)))
      }

      splitClass <- rep(names(genesPos), sapply(genesPos, adjustedLength))
      newGenesPos <- split(seq_along(splitClass), splitClass)
      sapply(newGenesPos, PositionMean, vector = bootRatio)
    }

  names(bootListSamples) <- basename(colnames(sampleMatrix))
  return(bootListSamples)
}

#calculate significance 
CheckSignificance <- function(bootList) {
  # package check complains
  i <- NULL
  j <- NULL
  sigTables =  foreach(i = seq_along(bootList)) %:%
    foreach(j = 1:ncol(bootList[[i]]), .combine = rbind) %do% {
      meanNoise <- mean(bootList[[i]][,j])
      lowerBound <- quantile(bootList[[i]][,j], 0.025, type = 1)
      upperBound <- quantile(bootList[[i]][,j], 1 - 0.025, type = 1)
      bounds <- c(lowerBound, meanNoise, upperBound)
      roundedBounds <- round(bounds, digits = 2)
      maybeAmplification <- roundedBounds[1] > 1
      maybeDeletion <- roundedBounds[3] < 1
      roundedSignificant <- maybeAmplification | maybeDeletion
      c(roundedBounds, roundedSignificant)
    }
  return(sigTables)
}

SignificantGenes <- function(sigList, genesPositionsIndex) {
  # package check complains
  i <- NULL
  sigGenes = foreach(i = seq_along(sigList)) %do% {
    names(genesPositionsIndex)[sigList[[i]][, 3] == 1]
  }
  return(sigGenes)
}

# #############################################################################
# # functions for background noise estimation
# #############################################################################
NonSignificantGeneIndex <- function(sigList, genesPositionsIndex) {
  # package check complains
  i <- NULL
  genePosNonSig = foreach(i = seq_along(sigList)) %do% {
    sigGenes = names(genesPositionsIndex)[sigList[[i]][, 3] == 1]
    RemSigGenes(genesPositionsIndex, sigGenes)
  }
  return(genePosNonSig)
}

AmplProbMultipeSamples <- function(genePosNonSig) {
  # package check complains
  i <- NULL
  amplWeights = foreach(i = seq_along(genePosNonSig)) %do% {
    AmplProb(genePosNonSig[[i]])
  }
  return(amplWeights)
}

Background <- function(geneNames,
                       samplesNormalizedReadCounts,
                       referenceNormalizedReadCounts,
                       bootList,
                       replicates = 1000) {

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
  # calculate the weight for each amplicon of significant genes
  amplWeights <- AmplProbMultipeSamples(genesPosNonSig)
  uniqueAmpliconNumbers <- NumberOfUniqueAmplicons(genesPositionsIndex)

  i <- NULL
  backgroundObject <- foreach(i = seq_along(genesPosNonSig)) %do% {
    IterateAmplNum(uniqueAmpliconNumbers,
                   ratioMatrix[unlist(genesPosNonSig[[i]]), i],
                   replicates = replicates,
                   probs = amplWeights[[i]])
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
    isSig <- (sigList[[i]][, 3] < 1) | (sigList[[i]][, 1] > 1)
    backgroundUp <- backgroundReport[[i]][, 3]
    backgroundDown <- backgroundReport[[i]][, 1]
    rMatGene <- ratioMatGene[, i]
    aboveNoise <- (rMatGene > 1 & (rMatGene) > backgroundUp) |
      (ratioMatGene[, i] < 1 &
         (ratioMatGene[, i]) < backgroundDown)

    dfTemp <- data.frame(meanRatio = ratioMatGene[, i],
                         lowerBoundBootstrapRatio = sigList[[i]][, 1],
                         meanBootstrapRatio = sigList[[i]][, 2],
                         upperBoundBootstrapRatio = sigList[[i]][, 3],
                         lowerNoise = backgroundReport[[i]][, 1],
                         meanNoise = backgroundReport[[i]][, 2],
                         upperNoise = backgroundReport[[i]][, 3],
                         significant = isSig,
                         aboveNoise = aboveNoise,
                         amplNum = as.vector(table(geneNames)),
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
                       "Passed")
    dfTemp
  }
  names(reportTables) <- seq_along(reportTables)
  return(reportTables)
}

BackgroundReport <- function(background, geneNames) {
  # because package generaton complains..
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
  geneCounts = elementLengths(genesPos)
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
  randomPos = sample(seq_along(ratios),numAmpl,prob = amplWeights)
  randomlySelectedRatios = ratios[randomPos]
  randomMean = mean(randomlySelectedRatios)
  return(randomMean)
}

SampleNoiseGenes <- function(numAmpl = 2,
                             ratios,
                             replicates = 100,
                             probs = NULL) {
  # now repeat the sampling for the selected number of
  # amplicons replicates
  sampleNoiseDistribution = replicate(replicates,
                                      SampleRatio(ratios,
                                                  numAmpl,
                                                  amplWeights = probs))
  #get the upper, mean and lower values of the noise distribution
  meanNoise <- mean(sampleNoiseDistribution)
  sdNoise <- sd(sampleNoiseDistribution)

  #define a mean function that accepts an index for boot
  #mean.fun <- function(dat, idx) mean(dat[idx], na.rm = TRUE)
  #perform the bootstrapping
  #bootOut <- boot(sampleNoiseDistribution,
  #                 mean.fun,
  #                 R=replicates,
  #                 sim="ordinary")
  #calculate the confidence interval
  #ci <- boot.ci(bootOut, conf = specificityLevel,
  #              type = "bca")
  #get the lower bound
  #lowerBound <- ci$bca[4]
  #get the upper bound
  #upperBound <- ci$bca[5]
  #specificityLevel <- 2

  lowerBound <- meanNoise +  qnorm(0.025) * sdNoise
  upperBound <- meanNoise + qnorm(1 - 0.025) * sdNoise

  sampledNoise = c(lowerBound, meanNoise, upperBound)
  names(sampledNoise) = paste0(c("Lower", "Mean", "Upper"), "Noise")
  return(sampledNoise)
}

#how many unique amplicon numbers are there
NumberOfUniqueAmplicons <- function(genesPos) {
  # calcUniqueAmpliconNumbers <-function(genesPos) {
  ampliconNumbers = elementLengths(genesPos)
  uniqueAmpliconNumbers = sort(unique(ampliconNumbers))
  return(uniqueAmpliconNumbers)
}

#now calculate the background for each of the unique amplicon numbers
IterateAmplNum <- function(uniqueAmpliconNumbers,
                           ratios,
                           replicates = 100,
                           probs = NULL) {
  # Needed because of package check complaints..
  i <- NULL
  noiseResults = foreach(i = seq_along(uniqueAmpliconNumbers)) %do% {
    sampledNoise = SampleNoiseGenes(uniqueAmpliconNumbers[i],
                                    ratios = ratios,
                                    replicates = replicates,
                                    probs = probs)
    sampledNoise
  }
  names(noiseResults) = as.character(uniqueAmpliconNumbers)
  return(noiseResults)
}

PlotBootstrapDistributions  <- function(bootList,
                                        reportTables,
                                        outputFolder = getwd(),
                                        sampleNames = NULL,
                                        save = FALSE,
                                        scale = 7) {
  i <- NULL
  plotList <- foreach(i = seq_along(bootList)) %do% {
    selSample <- i
    test <- as.factor(reportTables[[selSample]][, "Passed"])
    test <- revalue(test, c(`0` = "noChange",
                            `1` = "nonReliableChange",
                            `2` = "ReliableChange"))
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
    } else {
      filename <- sampleNames[selSample] 
    }
    test <- ggplot(df, aes(x = class, y = log2(ratios), fill = testsPassed)) +
      geom_violin() + ggtitle(filename) +
      theme(plot.title = element_text(lineheight = 0.8, face = "bold"),
            text = element_text(size = 15),
            axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(name = "CNV Reliability",
                        values = c("#56B4E9", "#CC79A7" ,"#D55E00"),
                        labels = c("0" = "Foo", "1" = "Bar")) +
      coord_cartesian(ylim = ylim1) +
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=20,face="bold")) + 
      scale_x_discrete("Gene Names") +
      #       geom_hline(yintercept=log2(1.5), color="#009E73") +
      #       geom_hline(yintercept=log2(0.5), color="#009E73") +
      geom_hline(yintercept=0, color="#009E73") 

    if(save == TRUE) {
      dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
      filepath <- paste0(outputFolder, "/", filename, "_plot.pdf")
      ggsave(filename = filepath,
             plot = test,
             height = 7.42 * 1,
             width = 8.11 * 2,
             limitsize = FALSE)
    }
    return(test)
  }
  return(plotList)
}

