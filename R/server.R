#library(shiny)
#library(shinyFiles)
#library(CNVPanelizer)
#library(shinyjs)
#library(openxlsx)
#library(gridExtra)

# library(openxlsx)
# library(ggplot2)
# library(plyr)
# library(NOISeq)
# library(foreach)
# library(GenomicRanges)
# library(ExomeDepth)
# library(exomeCopy)

#source("./Library.R")
options(shiny.maxRequestSize=6000*1024^2)


CNVPanelizerServer <- function(input, output, session) {
#shinyServer(function(input, output, session) {

# Hides somethings nice to check from GUI at development time..
#  debug = TRUE
  debug = FALSE

  outputDirectory <- CreateSessionOutputDirectory(session)
  inputDirectory <- CreateSessionInputDirectory(session)

  print(paste("outputDirectory >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", outputDirectory))
  print(paste("inputDirectory >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", inputDirectory))


################################################################################
# Relative to the bed file
################################################################################
# TODO can be improved.. or replaced if data table were configurable..
FixDataTable <- function(x, inputDirectory, type) {
  newX <- x
  newX$datapath = UploadedFilepaths(x, inputDirectory)
  numberOfFiles <- length(newX$datapath)
  newX$type = rep(type, numberOfFiles)
  if (!is.null(newX$size)) {
    newX$size <- sapply(newX$size, utils:::format.object_size, "auto")
    if (!debug) {
      datapathColumnIndex <- 4
      newX <- newX[, -c(datapathColumnIndex)]
    }
  }
  return(newX)
}

observeEvent(input$bedFilepath, {
  if (!is.null(input$bedFilepath)) {
    fixUploadedFilesNames(input$bedFilepath, inputDirectory)
    myDataTable <- read.csv(UploadedFilepaths(input$bedFilepath, inputDirectory), header= FALSE, skip= 1, sep = "\t")
    numberOfCols <- length(colnames(myDataTable))
    colnames(myDataTable) <- c("chr", "start", "end", rep("", numberOfCols-3))
    output$bedFile <- renderDataTable(myDataTable,
                                      options = list(dom = "", searching = TRUE))
    output$importedbedFilepath <- renderDataTable(FixDataTable(input$bedFilepath, inputDirectory, type = "bed"),
                                                  options = list(dom = "", searching = FALSE))
  }
})

################################################################################
# Relative to the uploaded sample bam files
################################################################################

observeEvent(input$samplesFilepaths, {
  fixUploadedFilesNames(input$samplesFilepaths, inputDirectory)
  output$importedSamplesFilepaths <- renderDataTable(
    FixDataTable(input$samplesFilepaths, inputDirectory, "bam"),
    options = list(dom = "", searching = TRUE))
})


################################################################################
# Relative to the uploaded Reference files
################################################################################
#observeEvent(input$referenceBamFilepaths, {
observeEvent(input$referenceFilepaths, {
# if (!is.null(input$referenceBamFilepaths)) {
  fixUploadedFilesNames(input$referenceFilepaths, inputDirectory)


  output$importedReferenceFilepaths <- renderDataTable(
    #  fixUploadedFilesNames(input$referenceFilepaths),
    #  FixedBamDataStructure(input$referenceFilepaths, inputDirectory),
    FixDataTable(input$referenceFilepaths, inputDirectory, type = "bam"),
    options = list(dom = "", searching = TRUE)
  )
})



RunCNVPanelizerAndUpdateGUI <- function(sampleBamFilepaths,
                                      referenceBamFilepaths,
                                      bedFilepath,
                                      amplColumnNumber,
                                      minimumMappingQuality,
                                      numberOfBootstrapReplicates,
                                      removePcrDuplicates,
                                      robust,
                                      backgroundSignificanceLevel) {

# IT WORKS WITH THIS SAMPLES !!!!!!
# sampleBamFilepaths <- c("D:\\data\\NGS\\Reference\\CRC Diagnostik\\697 N 4_1PG-196_IonXpress_094_rawlib.bam",
#                         "D:\\data\\NGS\\Reference\\CRC Diagnostik\\589-06_1PG-134_IonXpress_001_rawlib.bam",
#                         "D:\\data\\NGS\\Reference\\CRC Diagnostik\\589-06_1PG-134_IonXpress_001_rawlib.bam",
#                         "D:\\data\\NGS\\Reference\\CRC Diagnostik\\589-06_1PG-134_IonXpress_001_rawlib.bam",
#                         "D:\\data\\NGS\\Reference\\CRC Diagnostik\\589-06_1PG-134_IonXpress_001_rawlib.bam",
#                         "D:\\data\\NGS\\Reference\\CRC Diagnostik\\589-06_1PG-134_IonXpress_001_rawlib.bam",
#                         "D:\\data\\NGS\\Reference\\CRC Diagnostik\\589-06_1PG-134_IonXpress_001_rawlib.bam",
#                         "D:\\data\\NGS\\Reference\\CRC Diagnostik\\589-06_1PG-134_IonXpress_001_rawlib.bam",
#                         "D:\\data\\NGS\\Reference\\CRC Diagnostik\\589-06_1PG-134_IonXpress_001_rawlib.bam",
#                         "D:\\data\\NGS\\Reference\\CRC Diagnostik\\697 N 4_1PG-196_IonXpress_093_rawlib.bam")
# referenceBamFilepaths <- c("D:\\data\\NGS\\Reference\\CRC Diagnostik\\589-06_1PG-134_IonXpress_001_rawlib.bam", "D:\\data\\NGS\\Reference\\CRC Diagnostik\\697 N 4_1PG-196_IonXpress_093_rawlib.bam")
# bedFilepath <- "D:\\data\\NGS\\bed-files\\panels-with-chr\\CRC_CORRECTED2.bed"


# TEMPORARY THING TO PLAY AROUND...
#sampleBamFilepaths <- input$samplesFilepaths$datapath
#bedFilepath <- input$bedFilepath$datapath


# IT WORKS FOR QUICK TESTS..
# testFilesBaseDirectory <- "D:/wip/CNVPanelizerShiny/tmp/input/testing"
# bamFilepaths <- file.path(testFilesBaseDirectory,
#                           c("589-06_1PG-134_IonXpress_001_rawlib.bam",
#                             "697 N 4_1PG-196_IonXpress_093_rawlib.bam",
#                             "697 N 4_1PG-196_IonXpress_094_rawlib.bam",
#                             "697 N 4_1PG-196_IonXpress_093_rawlib.bam",
#                             "697 N 4_1PG-196_IonXpress_094_rawlib.bam",
#                             "697 N 4_1PG-196_IonXpress_093_rawlib.bam",
#                             "697 N 4_1PG-196_IonXpress_094_rawlib.bam",
#                             "697 N 4_1PG-196_IonXpress_095_rawlib.bam"))
# sampleBamFilepaths <- bamFilepaths[1:2]
# referenceBamFilepaths <- bamFilepaths[3:4]
# bedFilepath <- file.path(testFilesBaseDirectory, "CRC_CORRECTED2.bed")
# amplColumnNumber <- 4
#
#
# minimumMappingQuality <- 20
# numberOfBootstrapReplicates <- 3
# removePcrDuplicates <- TRUE
# robust <- TRUE
# backgroundSignificanceLevel <- 0.01
# outputDirectory = "D:/wip/CNVPanelizerShiny/tmp/94f180481edaee20cab99c95d7c960ed/output/plots/"

print(sampleBamFilepaths)
print(referenceBamFilepaths)
print(bedFilepath)

  CNVPanelizerResults <- CNVPanelizer(sampleBamFilepaths = sampleBamFilepaths,
                                      referenceBamFilepaths = referenceBamFilepaths,
                                      bedFilepath = bedFilepath,
                                      amplColumnNumber = amplColumnNumber,
                                      minimumMappingQuality = minimumMappingQuality, # TODO
                                      numberOfBootstrapReplicates = numberOfBootstrapReplicates, # TODO #numberOfBootstrapReplicates,
                                      removePcrDuplicates = removePcrDuplicates,
                                      robust = robust,
                                      backgroundSignificanceLevel = backgroundSignificanceLevel,
                                      outputDir = outputDirectory)
    # TODO why create another object..
    myOutput <- list()
    myOutput$outputDirectory <- outputDirectory
    myOutput$plots <- CNVPanelizerResults@plots
    return(myOutput)
}

observeEvent(input$Run, {

  notificationDuration <- 3

  if(is.null(input$samplesFilepaths) &
    (is.null(input$referenceFilepaths)) &
    (is.null(input$bedFilepath))) {
    showNotification(paste("Please upload all required files"), duration = notificationDuration)
    return()
  }

  if (!(!is.na(input$minimumMappingQuality) && !is.null(input$minimumMappingQuality))) {
      showNotification(paste("Please fill Min. Mapping Quality"), duration = notificationDuration)
      return()
  }

  if (!(!is.na(input$numberOfBootstrapReplicates) && !is.null(input$numberOfBootstrapReplicates))) {
      showNotification(paste("Please fill Number of Bootstrap Replicates"), duration = notificationDuration)
      return()
  }

  if (!(!is.na(input$backgroundSignificanceLevel) && !is.null(input$backgroundSignificanceLevel))) {
      showNotification(paste("Please fill Background significance Level"), duration = notificationDuration)
      return()
      }

  if (!(!is.na(input$bedColumn) && !is.null(input$bedColumn))) {
      showNotification(paste("Please fill Bed Column"), duration = notificationDuration)
      return()
      }

  EmptySessionOutputDirectory(session)
  updateTabsetPanel(session, "inTabset", selected = "reportTables")
  sampleBamFilepaths <- UploadedFilepaths(input$samplesFilepaths, inputDirectory)
  referenceBamFilepaths <- UploadedFilepaths(input$referenceFilepaths, inputDirectory)
  bedFilepath <- UploadedFilepaths(input$bedFilepath, inputDirectory)

  withProgress(message = 'Analysis in progress', {
    # not working!?!?
  #  hide("plot")
    result <- RunCNVPanelizerAndUpdateGUI(sampleBamFilepaths = sampleBamFilepaths,
                                          referenceBamFilepaths = referenceBamFilepaths,
                                bedFilepath = bedFilepath,
                                amplColumnNumber = input$bedColumn,
                                minimumMappingQuality = input$minimumMappingQuality, # TODO
                                numberOfBootstrapReplicates = input$numberOfBootstrapReplicates, # TODO
                                removePcrDuplicates = input$removePcrDuplicates,
                                robust = input$robust,
                                backgroundSignificanceLevel = input$backgroundSignificanceLevel)

    output$plot <- renderPlot({multiplot(plotlist = result$plots)}, height = length(result$plots) * 300)
    Sys.sleep(5) # by some strange reason, the plot does not display imediately... :S
  # not working..
  #  show("plot")
  })
})

# downloadHandler() takes two arguments, both functions.
# The content function is passed a filename as an argument, and
#   it should write out data to that filename.
output$downloadData <- downloadHandler(
  # This function returns a string which tells the client
  # browser what name to use when saving the file.
  filename = function() {
    paste("CNVPanelizer", "zip", sep = ".")},
  content = function(fname) {
  originalWorkingDirectory <- getwd()
  setwd(outputDirectory)
  print(paste("sessionToken : ", session$token))
  print(paste("getwd() : ", getwd()))
  a <- zip(zipfile = fname, files=c("xlsx", "plots"))
  setwd(originalWorkingDirectory)
}, contentType = "application/zip"
)
}
#)

