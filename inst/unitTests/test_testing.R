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
    sampleReadCounts <- matrix(1:kNumberOfReferenceMatrixLines)

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
}
