
# HIDE TABS TODO
# https://stackoverflow.com/questions/31703241/activate-tabpanel-from-another-tabpanel/31719425#31719425

# CHANGE Files
#https://stackoverflow.com/questions/36144499/shiny-fileinput-does-not-keep-file-name

#library(shiny)
#library(data.table)
#library(shinyFiles)
#library(CNVPanelizer)
#library(shinyjs)


CNVPanelizerUI <- function() {
shinyUI(
  fluidPage(
    tags$head(
      tags$link(rel = "shortcut icon", href = "favicon.png"),
      tags$style(HTML("hr {border-top: 1px solid #000000;}")),
      tags$style(HTML("div#plot img {width:  100%; height:  100%;}"))
    ),
    useShinyjs(),   # https://stackoverflow.com/questions/42795524/r-shiny-shinyjs-remove-plot-and-draw-it-again-if-button-is-clicked
    titlePanel(div("CNVPanelizer")),
#    titlePanel(img(src='CNVPanelizerLogo.png', height="60px"), "CNVPanelizer"),
#  textOutput("txt_file"),
    sidebarLayout(
      sidebarPanel(
        width = 2,
        fileInput("samplesFilepaths",
                  label="Samples Bam files",
                  multiple = TRUE,
                  accept = c(
                ".bam")),
        #https://stackoverflow.com/questions/19470426/r-shiny-add-tabpanel-to-tabsetpanel-dynamically-with-the-use-of-renderui
        tabPanel("Test",
          div(style="display: inline-block;font-weight: bold; width: 180px;",
                     p("Bed amplicon column")),
                 div(style="display: inline-block;width: 65px;",
                   numericInput(inputId = "bedColumn", label= "", value = 4, min = 1))),
        fileInput("bedFilepath", label="Bed file", multiple = FALSE, accept = c(".bed")),

        tabPanel("minimumMappingQuality",
           div(style="display: inline-block;font-weight: bold; width: 180px;",
               p("Min. Mapping Quality")),
           div(style="display: inline-block;width: 65px;",
            numericInput(inputId = "minimumMappingQuality", label= "", value = 20, min = 0))),
        fileInput("referenceFilepaths", label="Reference Bam files", multiple = TRUE, accept = c(".bam")),
        numericInput("numberOfBootstrapReplicates", "Bootstrap Replicates", 2, min = 1, max = 100000),     # TODO increase to 10000 by default...
        numericInput("backgroundSignificanceLevel", "Significance Level:", 0.05, min = 0.01, max = 0.1),
        selectInput("multipleComparisonCorrection", "Comparison Correction:", c("Bonferroni" = "bonferroni")),
        checkboxInput("removePcrDuplicates", "Remove Pcr Duplicates", FALSE),
        checkboxInput("robust", "Robust", TRUE),
      div(actionButton('Run', label = "OK", icon = icon("ok", lib = "glyphicon")), style="float:right;")
      ,withTags( # hack to fix margin issue..
        div(style = "margin-top: 80px")
      )),
      mainPanel(
        width = 10,
        tabsetPanel(id = "inTabset",
          tabPanel('Samples',
                   div(style = "margin-top: 30px"),
                   dataTableOutput("importedSamplesFilepaths")),
          tabPanel('Bed',
                   div(style = "margin-top: 30px"),
                   dataTableOutput("importedbedFilepath"),
                   div(style = "margin-top: 60px"),
                   dataTableOutput("bedFile")),
          tabPanel('Reference',
                   div(style = "margin-top: 30px"),
                   dataTableOutput("importedReferenceFilepaths")),
          tabPanel('Results', value = "analysis",
                  textOutput("espaco"),
          withTags(
            div(class = "someSpace",
                h3("  "))
          ),
        div(downloadButton('downloadData', 'Download'), style="float:right;"),
        withTags(
          div(style = "margin-top: 10px")
        ), plotOutput("plot", width = "100%", height="100%")))
      )
    )
  )
)


}
