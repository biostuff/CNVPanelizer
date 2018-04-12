multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# TODO improve ( A little fragile... to say the least..)
MoveFileTo <- function(from, to) {
  print(paste("from : ", from))
  print(paste("to : ", to))
#  todir <- dirname(to)
  todir <- unique(dirname(to))

  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
#  file.rename(from = from,  to = to)

  for (i in seq_along(from)) {
    file.rename(from = from[i],  to = to[i])
  }
}

UploadedFilepaths <- function(x, inputDirectory) {
  return(file.path(inputDirectory, x$name))
}

fixUploadedFilesNames <- function(x, inputDirectory = "") {
  if (is.null(x)) {
    return()
  }
  oldNames = x$datapath
  newNames = UploadedFilepaths(x, inputDirectory)
  MoveFileTo(from = oldNames, to = newNames)
  x$datapath <- newNames
  x
}

GetSessionOutputDirectory <- function(session) {
  outputDirectory <- file.path(getwd(), "tmp", session$token, "output")
  return(outputDirectory)
}

CreateSessionOutputDirectory <- function(session) {
  outputDirectory <- GetSessionOutputDirectory(session)
  dir.create(outputDirectory, recursive = TRUE)
  return(outputDirectory)
}

GetSessionInputDirectory <- function(session) {
  intputDirectory <- file.path(getwd(), "tmp", session$token, "input")
  return(intputDirectory)
}

CreateSessionInputDirectory <- function(session) {
  intputDirectory <- GetSessionInputDirectory(session)
  dir.create(intputDirectory, recursive = TRUE)
  return(intputDirectory)
}

EmptySessionOutputDirectory <- function(session) {
  print("temporary directory Contents before deleting...")
  print(list.files(GetSessionOutputDirectory(session)))
  unlink(GetSessionOutputDirectory(session), recursive = TRUE)
  CreateSessionOutputDirectory(session)
  print("temporary directory Contents after deleting...")
  print(list.files(GetSessionOutputDirectory(session)))
}


#runCNVPanelizerShiny <- function(dataset=NULL, port=8100,...) {
RunCNVPanelizerShiny <- function(port=8100) {
  options(shiny.maxRequestSize=1024^3) #allows uploading 1Gb file.
  addResourcePath(prefix="www", directoryPath=system.file("extdata/www/", package="CNVPanelizer"))

  # if( !is.null(dataset)) {
  #   if (! is(dataset, "rTResult") ) {
  #     stop("The dataset object must be a rTANDEM result object of class='rTResult'.")
  #   }
  # }

  #redefine the env of shinyTandemServer to the calling environment so that
  #shinyTandemServer() has access to the 'dataset' object.
  environment(shinyServer) <- environment()

#
#  app <- list(
##    ui= shinyUI(),
##    serverFunct=shinyTandemServer,
#    ui= shinyUI(),
#    server = shinyServer(),
#    port=port)
#


#  runApp(app)
#runApp()
shiny::shinyApp(ui = CNVPanelizerUI(), server = function(input,output,session) {CNVPanelizerServer(input, output, session)}, options = list(port = port, launch.browser = TRUE))

}

