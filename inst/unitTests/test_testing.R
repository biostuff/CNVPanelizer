test.bootList <- function() {
    kNumberOfReferenceMatrixLines <- 4
    kNumberOfReferenceMatrixColumns <- 4
    kNumberOfElements <- kNumberOfReferenceMatrixLines * 
                        kNumberOfReferenceMatrixColumns
    referenceValuesReadCounts <- 1:kNumberOfElements
    # Reference read counts matrix
    referenceReadCounts <- matrix(referenceValuesReadCounts,
                              nrow = kNumberOfReferenceMatrixLines,
                              ncol = kNumberOfReferenceMatrixColumns,
                              byrow = FALSE)
    # TODO Should this be a single column matrix?
#     sampleReadCounts <- 1:kNumberOfReferenceMatrixLines
    
    sampleReadCounts <- matrix(1:kNumberOfReferenceMatrixLines)
    
#    colnames(referenceReadCounts) <- paste0("ref",1:kNumberOfReferenceMatrixColumns)
#     rownames(referenceReadCounts) <- paste0("c:/somefile",
#                                       1:kNumberOfReferenceMatrixLines, ".bam")
#     colnames(sampleReadCounts) <- paste0("c:/somefile",
#                                        1:kNumberOfReferenceMatrixLines, ".bam")

    colnames(sampleReadCounts) <- paste0("c:/somefile", 1, ".bam")
    
    # amplicons indexes for each gene
    genesPositionsIndex = list("gene_1" = c(1),
                               "gene_2" = c(2, 3, 4))
    
    geneNames <- c("GENE1",rep("GENE2",3))
    
    set.seed(1)
    kNumberOfReplicates <- 1
    bootList = BootList(geneNames,
                        sampleReadCounts,
                        referenceReadCounts,
                        replicates = kNumberOfReplicates)

    resultForGENE1 <- bootList[[1]]["GENE1"][[1]]
    resultForGENE2 <- bootList[[1]]["GENE2"][[1]]

    tolerance <- 0.0001
    checkEquals(resultForGENE1, 0.1428571, tolerance = tolerance)
    checkEquals(resultForGENE2, 0.3666667, tolerance = tolerance)
  
#     # This values has small differents if you do the
#   ratio of the means or vice versa..
#     expectedResult <- list(data.frame("GENE1" = 0.1428571, "GENE2" = 0.3666667))
# 
#   
# 
#     checkEquals()
# 
#     print (bootList[[1]]["GENE1"][[1]])
# #     print(bootList)
#     print("__________________________")
#     print(expectedResult)
#     
#     # TODO implement this properly..
# #     checkEquals(bootList, expectedResult)
#     checkEquals(1, 1)
}

test.examples <- function() {
  checkEquals(1, 1)
#  checkEquals(6, factorial(3))
#  checkEqualsNumeric(6, factorial(3))
#  checkIdentical(6, factorial(3))
#  checkTrue(2 + 2 == 4, 'Arithmetic works')
#  checkException(log('a'), 'Unable to take the log() of a string')
}

test.deactivation <- function() {
#   DEACTIVATED('Deactivating this test function')
  checkEquals(1, 1)
}
