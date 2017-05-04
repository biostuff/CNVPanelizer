# test.descriptiveStatistics <- function() {
# 	distribution <- c(1, 2, 3, 4, 5)
# 	myDescriptiveStatistics <- DescriptiveStatistics(distribution)
# 	checkEquals(myDescriptiveStatistics$centralTendency == 3)
# 	checkEquals(myDescriptiveStatistics$variability == 1.581139)
# }
#
# test.descriptiveStatisticsRobust <- function() {
# 	distribution <- c(1, 2, 3, 4, 5)
# 	myDescriptiveStatistics <- DescriptiveStatistics(distribution, robust = TRUE)
# 	checkEquals(myDescriptiveStatistics$centralTendency == 3)
# 	checkEquals(myDescriptiveStatistics$variability == 1)
# }

test.bootListForSamplasAsNonMatrix <- function() {
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
	# TODO this is the important part of the test.. maybe is possible to reduce code duplication
	sampleReadCounts <- 1:kNumberOfReferenceMatrixLines
	# colnames(sampleReadCounts) <- paste0("c:/somefile", 1, ".bam")
		# amplicons indexes for each gene
	genesPositionsIndex = list("gene_1" = c(1),
														 "gene_2" = c(2, 3, 4))
	geneNames <- c("GENE1",rep("GENE2",3))
	set.seed(1)
	kNumberOfReplicates <- 1
	#getOption("RUnit")$silent
	checkException(BootList(geneNames,
													sampleReadCounts,
													referenceReadCounts,
													replicates = kNumberOfReplicates),
								 silent = TRUE)
}

test.bootListForSamplasAsMatrixWithNonNamedColumns <- function() {
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
	# TODO this is the important part of the test.. maybe is possible to reduce code duplication
	sampleReadCounts <- matrix(1:kNumberOfReferenceMatrixLines)
	#	colnames(sampleReadCounts) <- paste0("c:/somefile", 1, ".bam")

	# amplicons indexes for each gene
	genesPositionsIndex = list("gene_1" = c(1),
														 "gene_2" = c(2, 3, 4))

	geneNames <- c("GENE1",rep("GENE2",3))

	set.seed(1)
	kNumberOfReplicates <- 1

	#getOption("RUnit")$silent
	checkException(BootList(geneNames,
													sampleReadCounts,
													referenceReadCounts,
													replicates = kNumberOfReplicates),
								 silent = TRUE)
}





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

	#    resultForGENE1 <- bootList[[1]]["GENE1"][[1]]
	#    resultForGENE2 <- bootList[[1]]["GENE2"][[1]]

	resultForGENE1 <- bootList[[1]][1, "GENE1"]
	resultForGENE2 <- bootList[[1]][1, "GENE2"]

	print(paste("GENE 1: " , resultForGENE1))
	print(paste("GENE 2: " , resultForGENE2))

	tolerance <- 0.0001
	checkEquals(resultForGENE1, 0.125,
							checkNames = FALSE,
							tolerance = tolerance)
	checkEquals(resultForGENE2, 0.3302891,
							checkNames = FALSE,
							tolerance = tolerance)
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
	myAmpliconNames <- c("amp1", "amp2", "amp3", "amp4")
	normalizedReadCountsWithAmpliconNames <- CombinedNormalizedCounts(sampleReadCounts,
																																		referenceReadCounts,
																																		ampliconNames = myAmpliconNames)

	checkEquals(rownames(normalizedReadCountsWithAmpliconNames$samples),
							myAmpliconNames)
	checkEquals(rownames(normalizedReadCountsWithAmpliconNames$reference),
							myAmpliconNames)
}

test.RobustTRUE <- function() {
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
																replicates = kNumberOfReplicates,
																robust = TRUE)

	tolerance <- 0.00001
	checkTrue(all.equal(backgroundNoise[[1]]$`1`["LowerNoise"], c(LowerNoise=0.5401818), tolerance = tolerance))
	checkTrue(all.equal(backgroundNoise[[1]]$`1`["MeanNoise"], c(MeanNoise=0.5401818), tolerance = tolerance))
	checkTrue(all.equal(backgroundNoise[[1]]$`1`["UpperNoise"], c(UpperNoise=0.5401818), tolerance = tolerance))
	checkTrue(all.equal(backgroundNoise[[1]]$`3`["LowerNoise"], c(LowerNoise=0.7654648), tolerance = tolerance))
	checkTrue(all.equal(backgroundNoise[[1]]$`3`["MeanNoise"], c(MeanNoise=0.8365187), tolerance = tolerance))
	checkTrue(all.equal(backgroundNoise[[1]]$`3`["UpperNoise"], c(UpperNoise=0.9141681), tolerance = tolerance))
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
	# it should be the same as they should be generated from the same bed file
	#ampliconNames <- rownames(referenceReadCounts)

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
														 reportTables,
														 outputFolder = tempdir(),
														 save=TRUE)

	sampleNames <- paste0("sample", 1:4)
	PlotBootstrapDistributions(bootList,
														 reportTables,
														 outputFolder = tempdir(),
														 sampleNames = sampleNames,
														 save=TRUE)

}

#' test.BedToGenomicRanges <- function() {
#'
#' 	bedContent <- 'track name="IAD34679_Design" description="CoveredBases_AmpliSeqID_IAD34679" type=bedDetail color=77,175,74 priority=2 AmpliSeq_Version=1.2.9
#' 	chr1	27056172	27056291	GENE11	.	AMPL3573652750
#' 	chr1	27057875	27058001	GENE11	.	AMPL3498925712
#' 	chr1	27059149	27059276	GENE11	.	AMPL3573725051
#' 	chr2	27088638	27088749	GENE21	.	AMPL414842062
#' 	chr2	27100836	27100964	GENE21	.	AMPL3573753863
#' 	chr2	27101392	27101517	GENE22	.	AMPL3573628614
#' 	chr3	27105485	27105614	GENE31	.	AMPL414126459
#' 	chr3	27105839	27105963	GENE31	.	AMPL414225948
#' 	chr3	27106283	27106414	GENE31	.	AMPL4828584845
#' 	chrX	27106751	27106852	GENEX1	.	AMPL3573650239
#' 	chrX	115252187	115252313	GENEX1	.	AMPL4713167298
#' 	chrX	115256484	115256587	GENEX2	.	AMPL389446477
#' 	chrX	115258676	115258805	GENEX2	.	AMPL389117292
#' 	chrY	48023081	48023207	GENEY1	.	AMPL3562810493'
#'
#' 	bedFilepath <- file.path(tempdir(), "myBed.bed")
#' 	writeLines(bedContent, bedFilepath)
#'
#' 	genomicRangesFromBed <- BedToGenomicRanges(bedFilepath,
#' 																						 4,
#' 																						 split = "_",
#' 																						 dropChromossomes = c("chrX"))
#'
#' 	metadataFromGenomicRanges <- elementMetadata(genomicRangesFromBed)
#' 	geneNames = metadataFromGenomicRanges["geneNames"][,1]
#' 	ampliconNames = metadataFromGenomicRanges["ampliconNames"][,1]
#'
#' 	checkEquals(geneNames, c("GENE11",
#' 													 "GENE11",
#' 													 "GENE11",
#' 													 "GENE21",
#' 													 "GENE21",
#' 													 "GENE22",
#' 													 "GENE31",
#' 													 "GENE31",
#' 													 "GENE31",
#' 													 "GENEY1"))
#'
#'
#' 	#
#' 	# TODO Testing for doReduce = FALSE (should be in another test)
#' 	#
#' 	genomicRangesFromBed <- BedToGenomicRanges(bedFilepath,
#' 																						 4,
#' 																						 doReduce = FALSE,
#' 																						 split = "_",
#' 																						 dropChromossomes = c("chrX"))
#'
#' 	metadataFromGenomicRanges <- elementMetadata(genomicRangesFromBed)
#' 	geneNames = metadataFromGenomicRanges["geneNames"][,1]
#' 	ampliconNames = metadataFromGenomicRanges["ampliconNames"][,1]
#'
#' 	checkEquals(geneNames, c("GENE11",
#' 													 "GENE11",
#' 													 "GENE11",
#' 													 "GENE21",
#' 													 "GENE21",
#' 													 "GENE22",
#' 													 "GENE31",
#' 													 "GENE31",
#' 													 "GENE31",
#' 													 "GENEY1"))
#'
#' }
#'
#' # test.adjustedLength <- function() {
#' #   param <- c(1, 2, 3)
#' #   lengthParam <- length(param)
#' #   checkEquals(lengthParam, 3)
#' #   sqrtParam <- sqrt(lengthParam)
#' #   checkEquals(sqrtParam, 1.732051)
#' #   roundSqrtParam <- round(sqrtParam)
#' #   checkEquals(roundSqrtParam, 2)
#' # }
#' #
#' # test.ratiosMean <- function() {
#' #   checkEquals(ratiosMean(c(1, 2, 3)), 1.817121)
#' # }
#'
#' test.readCountsFromBam <- function() {
#' 	# https://www.biostars.org/p/150010/
#'
#' 	bedHeaderContent <- 'track name="IAD34679_Design" description="CoveredBases_AmpliSeqID_IAD34679" type=bedDetail color=77,175,74 priority=2 AmpliSeq_Version=1.2.9'
#' 	bedBodyContent <- 'chr1	1	100	GENE11	.	AMPL3573652750
#' 	chr1	300	1000	GENE11	.	AMPL3498925712
#' 	chr2	100	500	GENE21	.	AMPL414842062
#' 	chr2	1000	2000	GENE21	.	AMPL3573753863'
#' 	bedContentHasNoEmptySpaces <- !grepl(" ", bedBodyContent)
#' 	checkTrue(bedContentHasNoEmptySpaces, msg="The only alowed column separator is the TAB (\t), and white spaces were found at the bed file")
#' 	bedContent <- paste0(bedHeaderContent, "\n", bedBodyContent)
#' 	bedFilepath <- file.path(tempdir(), "myBed.bed")
#' 	writeLines(bedContent, bedFilepath)
#' 	genomicRangesFromBed <- BedToGenomicRanges(bedFilepath,
#' 																						 4,
#' 																						 split = "_")
#' 	samContent <- "@SQ	SN:chr1	LN:2000
#' 	@SQ	SN:chr2	LN:3000
#' 	r001	0	chr1	1	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r002	16	chr1	400	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r003	16	chr1	410	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r004	0	chr1	450	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r005	16	chr1	800	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r006	0	chr2	900	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r007	16	chr2	1001	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r008	0	chr2	1500	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r009	0	chr2	1900	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r010	0	chr2	2050	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK"
#' 	samContentHasNoEmptySpaces <- !grepl(" ", samContent)
#' 	checkTrue(samContentHasNoEmptySpaces,
#' 						msg="The only alowed separator is the TAB (\t), and white spaces were found at the sam file")
#'
#' 	samFilepath <- file.path(tempdir(), "mySam.sam")
#' 	writeLines(samContent, samFilepath)
#' 	bamFilepathWithoutExtension <- gsub(pattern = "\\.sam$", "", samFilepath)
#' 	bamFilepath <- paste0(bamFilepathWithoutExtension, ".bam")
#' 	asBam(samFilepath,
#' 				destination = bamFilepathWithoutExtension,
#' 				overwrite = TRUE)
#'
#' 	sampleFilename <- "example"
#' 	readCountsFromBam <- ReadCountsFromBam(c(bamFilepath),
#' 																				 sampleNames = sampleFilename,
#' 																				 genomicRangesFromBed,
#' 																				 #                    ampliconNames = NULL,
#' 																				 removeDup = FALSE)
#'
#' 	readCountsResult = matrix(c(1, 4, 0, 3),
#' 														nrow=4,
#' 														ncol=1)
#' 	colnames(readCountsResult) = sampleFilename
#' 	#TODO why cant i just to this.. checkEquals(readCountsFromBam, readCountsResult)
#' 	checkEquals(as.vector(readCountsFromBam), as.vector(readCountsResult))
#' }
#'
#' test.IndexMultipleBams <- function() {
#' 	temporaryDirectory <- tempdir()
#' 	samContent <- "@SQ	SN:chr1	LN:2000
#' 	@SQ	SN:chr2	LN:3000
#' 	r001	0	chr1	1	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r002	16	chr1	400	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r003	16	chr1	410	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r004	0	chr1	450	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r005	16	chr1	800	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r006	0	chr2	900	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r007	16	chr2	1001	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r008	0	chr2	1500	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r009	0	chr2	1900	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK
#' 	r010	0	chr2	2050	30	20M	*	0	0	NNNNNNNNNNNNNNNNNNNN	KKKKKKKKKKKKKKKKKKKK"
#' 	samContentHasNoEmptySpaces <- !grepl(" ", samContent)
#' 	checkTrue(samContentHasNoEmptySpaces,
#' 						msg="The only alowed separator is the TAB (\t), and white spaces were found at the sam file")
#'
#' 	samFilepath <- file.path(temporaryDirectory, "mySam.sam")
#' 	writeLines(samContent, samFilepath)
#' 	bamFilepathWithoutExtension <- gsub(pattern = "\\.sam$", "", samFilepath)
#' 	bamFilepath <- paste0(bamFilepathWithoutExtension, ".bam")
#' 	baiFilepath <- paste0(bamFilepathWithoutExtension, ".bam.bai")
#' 	asBam(samFilepath,
#' 				destination = bamFilepathWithoutExtension,
#' 				overwrite = TRUE,
#' 				indexDestination = FALSE)
#'
#' 	IndexMultipleBams(bamFilepath)
#'
#' 	checkTrue(file.exists(baiFilepath), msg = paste("bai file was not generated at", bamFilepath))
#' }

test.WriteListToXLSXandReadXLSXToList <- function() {
	temporaryDirectory <- tempdir()
	myDataFrame <- as.data.frame(matrix(1:9, 3, 3))
	dataFrameName <- "SomeDataFrame"
	myList <- list(myDataFrame)
	names(myList) <- dataFrameName
	filepath <- file.path(temporaryDirectory, "samples.xlsx")
	WriteListToXLSX(myList, filepath)
	checkTrue(file.exists(filepath), msg = paste("xlsx file was not generated at", filepath))
	otherList <- ReadXLSXToList(filepath)
	checkEquals(names(myList), names(otherList))
	checkTrue(all(myList[[dataFrameName]] == otherList[[dataFrameName]]))
}
