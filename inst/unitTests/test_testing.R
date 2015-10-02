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

    print(paste("GENE 1: " , resultForGENE1))
    print(paste("GENE 2: " , resultForGENE2))

    tolerance <- 0.0001
    checkEquals(resultForGENE1, 0.125, tolerance = tolerance)
    checkEquals(resultForGENE2, 0.3302891, tolerance = tolerance)
}

test.ReportTablesWithSingleSample <- function() {
  
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
  kNumberOfReplicates <- 1000
  bootList = BootList(geneNames,
                      sampleReadCounts,
                      referenceReadCounts,
                      replicates = kNumberOfReplicates)
  
  normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts, referenceReadCounts)
  
  samplesNormalizedReadCounts = normalizedReadCounts["samples"][[1]]
  referenceNormalizedReadCounts = normalizedReadCounts["reference"][[1]]
  
  backgroundNoise <- Background(geneNames,
                                samplesNormalizedReadCounts,
                                referenceNormalizedReadCounts,
                                bootList,
                                replicates = kNumberOfReplicates)

  # if the test crashes it is wrong number of dimensions here..
  reportTables <- ReportTables(geneNames,
                               samplesNormalizedReadCounts,
                               referenceNormalizedReadCounts,
                               bootList,
                               backgroundNoise)
}

test.PlotBootstrap <- function() {
  referenceReadCounts <- as.matrix(read.table(header = TRUE, text = "
r1 r2 r3 r4
gene1_1 1 2 1 2
gene2_1 3 2 2 2
gene2_2 1 3 3 1
gene3_1 2 3 1 2
gene3_2 2 2 2 1
gene3_3 1 1 2 2
"))

  sampleReadCounts <- as.matrix(read.table(header = TRUE, text = "
s1 s2 s3 s4
gene1_1 2 5 2 1
gene2_1 2 2 2 4
gene2_2 2 3 2 4
gene3_1 2 4 3 6
gene3_2 2 5 2 6
gene3_3 1 6 2 6
"))

  ampliconNames <- rownames(sampleReadCounts)
  ampliconNames <- rownames(referenceReadCounts)
  listOfAmpliconNames <- strsplit(ampliconNames, split="_")
  geneNames  <- unlist(listOfAmpliconNames)[ c(TRUE, FALSE) ]

  set.seed(1)
  kNumberOfReplicates <- 10000
  bootList = BootList(geneNames,
                      sampleReadCounts,
                      referenceReadCounts,
                      replicates = kNumberOfReplicates)

  normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts, referenceReadCounts)

  samplesNormalizedReadCounts = normalizedReadCounts["samples"][[1]]
  referenceNormalizedReadCounts = normalizedReadCounts["reference"][[1]]

  backgroundNoise <- Background(geneNames,
                                samplesNormalizedReadCounts,
                                referenceNormalizedReadCounts,
                                bootList,
                                replicates = kNumberOfReplicates)

  reportTables <- ReportTables(geneNames,
                               samplesNormalizedReadCounts,
                               referenceNormalizedReadCounts,
                               bootList,
                               backgroundNoise)

  # From the first sample no gene had either nonReliable or Reliable changes
  checkEquals(reportTables[[1]][,"Passed"], c(0, 0, 0))

  # From the second sample only the first and third gene had Reliable changes, second gene had noChange
  checkEquals(reportTables[[2]][,"Passed"], c(2, 0, 2))

gene3NonReliableChange <- read.table(header = TRUE, text = "
      Signif. AboveNoise
gene1   FALSE      FALSE
gene2   FALSE      FALSE
gene3    TRUE      FALSE")
  checkEquals(reportTables[[3]][,c("Signif.","AboveNoise")], gene3NonReliableChange)

gene1NoChange_gene2NonReliableChange_gene3ReliableChange <- read.table(
header = TRUE, text = "
      Signif. AboveNoise
gene1   FALSE      FALSE
gene2    TRUE      FALSE
gene3    TRUE       TRUE")
checkEquals(reportTables[[4]][,c("Signif.","AboveNoise")],
            gene1NoChange_gene2NonReliableChange_gene3ReliableChange)

  PlotBootstrapDistributions(bootList,
                            reportTables)
}
