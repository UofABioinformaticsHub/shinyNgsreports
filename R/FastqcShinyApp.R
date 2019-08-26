#' Run fastQC shiny app
#'
#' @description Returns a shiny app
#' interface to parse many fastQC
#' objects using the package ngsReports
#'
#' @details Currently some plots can
#' take a while to render if the
#' \code{FastqcDataList} passed to
#' \code{fastqcInput} has many elements
#'
#' @param fastqcInput can be a \code{FastqcFileList},
#'  \code{fastqcDataList},
#' or simply a \code{character} vector
#'  of paths to fastqc files.
#'
#'
#' @return UI data for fastQC shiny.
#'
#' @import shinydashboard
#' @import ngsReports
#' @importFrom magrittr %>%
#' @importFrom plotly layout
#' @importFrom plotly plotlyOutput
#' @importFrom plotly renderPlotly
#' @importFrom plotly event_data
#' @importFrom shiny br
#' @importFrom shiny h1
#' @importFrom shiny h5
#' @importFrom shiny observe
#' @importFrom shiny radioButtons
#' @importFrom shiny checkboxInput
#' @importFrom shiny column
#' @importFrom shiny reactiveValues
#' @importFrom shiny observeEvent
#' @importFrom shiny icon
#' @importFrom shiny htmlOutput
#' @importFrom shiny selectInput
#' @importFrom shiny plotOutput
#' @importFrom shiny renderUI
#' @importFrom shiny renderPrint
#' @importFrom shiny runApp
#' @importFrom shiny textOutput
#' @importFrom shiny renderText
#' @importFrom shiny reactive
#' @importFrom shiny withProgress
#' @importFrom shinyDirectoryInput choose.dir
#' @importFrom shinyDirectoryInput directoryInput
#' @importFrom shinyDirectoryInput readDirectoryInput
#' @importFrom shinyDirectoryInput updateDirectoryInput
#' @importFrom shinyFiles shinyFilesButton
#' @importFrom shinyFiles shinyFileChoose
#' @importFrom shinyFiles parseFilePaths
#' @importFrom shinyFiles parseDirPath
#' @importFrom shinyFiles getVolumes
#'
#' @examples
#' \dontrun{
#' ## Get the files included with the package ngsReports
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fl <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- FastqcDataList(fl)
#'
#' # Run the Shiny app
#' fastqcShiny(fdl)
#' 
#' ## Or can simply be run using and input files into the shiny app
#' fastqcShiny()
#' 
#' }
#'
#' @export
#' @rdname fastqcShiny
#'

fastqcShiny <- function(fastqcInput = NULL) {
  ## check if the initial value of fastqcInput is null and check class of
  ## object passed along with length
  if (!is.null(fastqcInput)) {
    stopifnot(length(fastqcInput) > 1)
    stopifnot(
      any(is(fastqcInput, "character"), is(fastqcInput, "FastqcDataList"))
    )
  }
  
  ## start rendering the dashboard gui
  body <- dashboardBody(
    tabItems(
      tabItem(
        tabName = "BS",
        column(
          width = 2,
          box(
            h5("Choose FastQC Report:"),
            shinyFilesButton(
              id = "files",
              label = "Choose files",
              multiple = TRUE,
              title = ""
            ),
            br(),
            textOutput("report"),
            br(),
            checkboxInput("Sumcluster", "Cluster", value = TRUE),
            collapsible = TRUE,
            width = NULL,
            title = "Options"
          )
        ),
        box(
          h1("Summary of fastQC Flags"),
          h5(
            "Heatmap of fastQC flags (PASS/WARN/FAIL) for each fastQC report"
          ),
          plotlyOutput("SummaryFlags", height = 700),
          width = 10
        )
      ),
      tabItem(
        tabName = "TS",
        box(
          checkboxInput("showDup", "Show Duplicated?", value = FALSE),
          collapsible = TRUE,
          width = 2,
          title = "Options"
        ),
        box(
          h1("Total Sequences"),
          h5("Total number of unique and duplicated reads in each sample"),
          plotlyOutput("ReadTotals"),
          width = 10
        )
      ),
      tabItem(
        tabName = "BQ",
        column(
          width = 2,
          box(
            radioButtons(
              inputId = "BQplotValue",
              label = "Base Quality",
              choices = c("Mean", "Median"),
              selected = "Mean"
            ),
            checkboxInput("BQcluster", "Cluster", value = TRUE),
            collapsible = TRUE,
            width = NULL,
            title = "Options"
          ),
          valueBoxOutput("BQboxP", width = NULL),
          valueBoxOutput("BQboxW", width = NULL),
          valueBoxOutput("BQboxF", width = NULL)
        ),
        box(
          h1("Per Base Sequence Quality"),
          h5(
            paste(
              "Per base sequence quality in each sample can either view",  
              "mean or median for each cycle"
            )
          ),
          h5("Click sidebar on heatmap to change line plots"),
          plotlyOutput("baseQualHeatmap"),
          br(),
          plotlyOutput("BaseQualitiesSingle"),
          width = 10
        )
      ),
      tabItem(
        tabName = "SQ",
        column(
          width = 2,
          box(
            radioButtons(
              inputId = "SQType",
              label = "Sequence Quality",
              choices = c("Frequency", "Counts"),
              selected = "Frequency"
            ),
            checkboxInput("SQcluster", "Cluster", value = TRUE),
            collapsible = TRUE,
            width = NULL,
            title = "Options"
          ),
          valueBoxOutput("SQboxP", width = NULL),
          valueBoxOutput("SQboxW", width = NULL),
          valueBoxOutput("SQboxF", width = NULL)
        ),
        box(
          h1("Per Sequence Quality Scores"),
          h5(
            paste(
              "Per base sequence quality in each sample, can either",
              "view mean or median for each cycle"
            )
          ),
          h5("Click sidebar on heatmap to change line plots"),
          plotlyOutput("seqQualHeatmap"),
          br(),
          plotlyOutput("SeqQualitiesSingle"),
          width = 10
        )
      ),
      tabItem(
        tabName = "SC",
        column(
          width = 2,
          box(
            checkboxInput("SCcluster", "Cluster", value = TRUE),
            collapsible = TRUE,
            width = NULL,
            title = "Options"
          ),
          valueBoxOutput("SCboxP", width = NULL),
          valueBoxOutput("SCboxW", width = NULL),
          valueBoxOutput("SCboxF", width = NULL)
        ),
        box(
          h1("Per Base Sequence Content"),
          h5(
            paste(
              "Per base sequence content in each sample.", 
              "Colours at each base indicate sequence bias"
            )
          ),
          h5("G = Black, A = Green, T = Red, C = Blue"),
          plotlyOutput("SCHeatmap"),
          br(),
          plotlyOutput("SCsingle"),
          width = 10
        )
      ),
      tabItem(
        tabName = "GC",
        column(
          width = 2,
          box(
            checkboxInput("GCcluster", "Cluster", value = TRUE),
            checkboxInput("theoreticalGC",
                          "Normalize To Theoretical GC",
                          value = FALSE),
            htmlOutput("theoreticalGC"),
            htmlOutput("GCspecies"),
            collapsible = TRUE,
            width = NULL,
            title = "Options"
          ),
          valueBoxOutput("GCboxP", width = NULL),
          valueBoxOutput("GCboxW", width = NULL),
          valueBoxOutput("GCboxF", width = NULL)
        ),
        box(
          h1("Per Sequence GC Content"),
          h5(
            paste(
              "GC content (%) in sample can either view total count", 
              "or frequency"
            )
          ),
          h5("Click sidebar on heatmap to change line plots"),
          plotlyOutput("GCheatmap"),
          br(),
          plotlyOutput("GCSingle"),
          width = 10
        )
      ),
      tabItem(
        tabName = "NC",
        column(
          width = 2,
          box(
            checkboxInput("Ncluster", "Cluster", value = TRUE),
            collapsible = TRUE,
            width = NULL,
            title = "Options"
          ),
          valueBoxOutput("NCboxP", width = NULL),
          valueBoxOutput("NCboxW", width = NULL),
          valueBoxOutput("NCboxF", width = NULL)
        ),
        box(
          h1("Per base N content"),
          h5("N content (%) in sample"),
          h5("If dendrogram is truncated double click on dendrogram to resize"),
          plotlyOutput("NCheatmap"),
          br(),
          plotlyOutput("NCsingle"),
          width = 10
        )
      ),
      tabItem(
        tabName = "SLD",
        column(
          width = 2,
          box(
            radioButtons(
              inputId = "SLType",
              label = "Value to plot",
              choices = c("Frequency", "Counts"),
              selected = "Frequency"
            ),
            checkboxInput("SLcluster", "Cluster", value = TRUE),
            collapsible = TRUE,
            width = NULL,
            title = "Options"
          ),
          valueBoxOutput("SLDboxP", width = NULL),
          valueBoxOutput("SLDboxW", width = NULL),
          valueBoxOutput("SLDboxF", width = NULL)
        ),
        box(
          h1("Sequence Length Distribution"),
          h5(
            paste(
            "Sequence length distribution in each sample, can either view", 
            "total count or frequency"
            )
          ),
          h5("Click sidebar on heatmap to change line plots"),
          plotlyOutput("SLHeatmap"),
          br(),
          plotlyOutput("SLSingle"),
          width = 10
        )
      ),
      tabItem(
        tabName = "SDL",
        column(
          width = 2,
          box(
            checkboxInput("Dupcluster", "Cluster", value = TRUE),
            collapsible = TRUE,
            width = NULL,
            title = "Options"
          ),
          valueBoxOutput("SDLboxP", width = NULL),
          valueBoxOutput("SDLboxW", width = NULL),
          valueBoxOutput("SDLboxF", width = NULL)
        ),
        box(
          h1("Sequence Duplication Levels"),
          h5("Sequence duplication in each sample"),
          h5("Click sidebar on heatmap to change line plots"),
          plotlyOutput("DupHeatmap"),
          br(),
          plotlyOutput("DupSingle"),
          width = 10
        )
      ),
      tabItem(
        tabName = "OS",
        column(
          width = 2,
          box(
            checkboxInput("OScluster", "Cluster", value = TRUE),
            h5("Export Overrepresented Sequences"),
            directoryInput('ORdirs', label = 'Select a directory', value = '~'),
            collapsible = TRUE,
            width = NULL,
            title = "Options"
          ),
          valueBoxOutput("OSboxP", width = NULL),
          valueBoxOutput("OSboxW", width = NULL),
          valueBoxOutput("OSboxF", width = NULL)
        ),
        box(
          h1("Overrepresented Sequences"),
          h5("Origin of Overrepresented sequences within each sample"),
          plotlyOutput("OSummary"),
          br(),
          plotlyOutput("OSsingle"),
          width = 10
        )
      ),
      tabItem(
        tabName = "AC",
        column(
          width = 2,
          box(
            textInput(
              "ACtype",
              "Regular expression matching the adapter(s) to be plotted",
              value = "Total"
            ),
            checkboxInput("ACcluster", "Cluster", value = TRUE),
            collapsible = TRUE,
            width = NULL,
            title = "Options"
          ),
          valueBoxOutput("ACboxP", width = NULL),
          valueBoxOutput("ACboxW", width = NULL),
          valueBoxOutput("ACboxF", width = NULL)
        ),
        box(
          h1("Adapter content"),
          h5("Adapter content (%) across all reads"),
          plotlyOutput("ACheatmap"),
          br(),
          plotlyOutput("ACsingle"),
          width = 10
        )
      ),
      tabItem(
        tabName = "KC",
        column(
          width = 2,
          box(
            checkboxInput("KMcluster", "Cluster", value = TRUE),
            collapsible = TRUE,
            width = NULL,
            title = "Options"
          ),
          valueBoxOutput("KCboxP", width = NULL),
          valueBoxOutput("KCboxW", width = NULL),
          valueBoxOutput("KCboxF", width = NULL)
        ),
        box(
          h1("Kmer Content"),
          h5(
            paste(
              "Total Identified Kmer Count by Position.",
              "\nPlease load a file to see the top 6 Kmers."
            )
          ),
          plotlyOutput("Kheatmap"),
          br(),
          plotlyOutput("Ksingle"),
          width = 10
        )
      ),
      tabItem(
        tabName = "HTML",
        column(
          width = 2,
          box(
            radioButtons(
              inputId = "omicsType",
              label = "-omic Type",
              choices = c("Genome", "Transcriptome"),
              selected = "Genome"
            ),
            htmlOutput("sequencedSpecies"),
            h5("Output report for files"),
            directoryInput('dirs', label = 'Select a directory', value = '~'),
            collapsible = TRUE,
            width = NULL,
            title = "Options"
          )
        ),
        box(
          h1("Output HTML Report Using the Default Template "),
          h5(
            paste(
            "Select the type of omic data used in your study",
            "and the most suitable organism"
            )
          ),
          h5(
            "from the dropdown list. Upon
            selecting the applicable omic
            and species, select the directory"
          ),
          h5(
            "containing the FASTQC files
            you wish to make the log for.
            Currently even if files have
            been loaded"
          ),
          h5(
            "into the Shiny app the folder
            containing the data must still
            be selected."
          ),
          width = 10
        )
      )
    )
  )
  
  
  ### ui page
  ui <- dashboardPage(
    dashboardHeader(title = "fastqcRShiny"),
    dashboardSidebar(
      #### gender the menu items
      sidebarMenu(
        menuItem(text = "Summary", tabName = "BS"),
        menuItem(text = "Total Sequences", tabName = "TS"),
        menuItemOutput("BQflag"),
        menuItemOutput("SQflag"),
        menuItemOutput("SCflag"),
        menuItemOutput("GCflag"),
        menuItemOutput("NCflag"),
        menuItemOutput("SLDflag"),
        menuItemOutput("SDLflag"),
        menuItemOutput("OSflag"),
        menuItemOutput("ACflag"),
        menuItemOutput("KCflag"),
        menuItem(text = "Output HTML Report", tabName = "HTML")
        
        
      ),
      width = 250
    ),
    body
  )
  
  
  
  ### backend
  server <- function(input, output, session) {
    # set up reactives
    values <- reactiveValues()
    
    observeEvent(input$files, {
        fastqcInput <<- NULL
    })
    
    ## wait timer
    #rective function repsonsible for loading in the selected files or just using the fdl supplied
    data <- reactive({
        
        
      input$files
        
      data <- fastqcInput
      #check that input$files is empty (ie no files are selected)
      
      withProgress(
          min = 0,
          max = 1,
          value = 0.8,
          message = "Reading in files.",
          {
      if(!length(fastqcInput)){
        # get the volumes when the shiny button is clicked
         
        volumes <- getVolumes()
        #choose the files from the root directory of the current volume and show only show zip files in the loadout
        shinyFileChoose(
          input,
          "files",
          roots = volumes,
          session = session,
          filetypes = "zip",
          defaultPath = getwd()
        )
        #get the metadata for the file that is selected, files is the element bade in the body script
        fileSelected <- parseFilePaths(volumes, input$files)
        # make the selected file a character vector
        data <- as.character(fileSelected$datapath)
        # import the selected file(s)
        #selectedData <- fileSelected)
        
      }
      else{
          fastqcInput <- NULL
        #if a character string is provided, it will import the data for you
        if (class(data) == "character") {
            FastqcDataList(data)
        } else data
      }
          }
      )
    })
    
    
    
    #render the UI for gcTheoretical
    output$sequencedSpecies <- renderUI({
      selectInput(
        "omicSpecies",
        "Select species",
        choices = gcAvail(ngsReports::gcTheoretical, type = input$omicsType)$Name,
        selected = "Hsapiens"
      )
      
    })
    
    
    species <- observeEvent(input$omicSpecies, {
      values$omicSpecies <- input$omicSpecies
    })
    
    #export overrepresented
    
    observeEvent(
        ignoreNULL = TRUE,
        eventExpr = {
            input$ORdirs
        },
        handlerExpr = {
            if (input$ORdirs > 0) {
                
                path <- choose.dir(default = readDirectoryInput(session, "ORdirs"))
                withProgress(
                    min = 0,
                    max = 1,
                    value = 0.8,
                    message = "Exporting Overrepresented Sequences",
                    {
                        exportOverrepresented(
                            data(),
                            path = paste0(path, "/OverrepSequences", "-", Sys.Date()),
                            n = 10,
                            noAdapters = TRUE
                        )
                    }
                )
                
                # update the widget value
                updateDirectoryInput(session, "ORdirs", value = path)
            }
        }
    )
    
    # expOS <- reactive({
    #   volumes <- shinyFiles::getVolumes()
    #   shinyFiles::shinyDirChoose(input, "dirOS",
    #                              roots = volumes, session = session)
    #   dirSelected <- shinyFiles::parseDirPath(volumes, input$dirOS)
    #   as.character(dirSelected)
    # })
    # 
    # observe({
    #   if (length(expOS())) {
    #     exportOverrepresented(
    #       data(),
    #       path = paste0(expOS(), "/OverrepSequences", "-", Sys.Date()),
    #       n = 10,
    #       noAdapters = TRUE
    #     )
    #   }
    # })
    
    # # export HTML
    # dir <- reactive({
    #   input$dirs
    #   volumes <- getVolumes()
    #   shinyFiles::shinyDirChoose(
    #     input = input, id = "dirs", roots = volumes, session = session, defaultPath = "~/"
    #   )
    #   dirSelected <- shinyFiles::parseDirPath(volumes, input$dirs)
    #   as.character(dirSelected)
    # })
    # 
    # observe({
    #   dir()
    #   if (length(dir())) {
    #     withProgress(
    #       min = 0,
    #       max = 1,
    #       value = 0.8,
    #       message = "Writing report",
    #       {
    #         writeHtmlReport(dir(),
    #                         species = values$omicSpecies,
    #                         gcType = input$omicsType, overwrite = TRUE)
    #       }
    #     )
    #     output$report2 <- renderText("Done!")
    #   }
    # })
    # 
    
    #### new export html
    observeEvent(
        ignoreNULL = TRUE,
        eventExpr = {
            input$dirs
        },
        handlerExpr = {
            if (input$dirs > 0) {
                
                path <- choose.dir(default = readDirectoryInput(session, "dirs"))
                withProgress(
                    min = 0,
                    max = 1,
                    value = 0.8,
                    message = "Writing report",
                    {
                        writeHtmlReport(path,
                                        species = values$omicSpecies,
                                        gcType = input$omicsType, overwrite = TRUE)
                    }
                )
                output$report2 <- renderText("Done!")
                
                
                # update the widget value
                updateDirectoryInput(session, "dirs", value = path)
            }
        }
    )
    
    
    #### dynamic tabs to show pass warn Fail ########
    
    # render the menu item for the Summary tab
    
    output$BQflag <- renderMenu({
      if (!is.null(fastqcInput) | length(input$files) > 1) {
        flags <- getSummary(data())
        Category <- c()
        flags <- subset(flags, Category == "Per base sequence quality")
        
        items <- .menuItemLogic(flags = flags)
        
        values$BQflag <- items[[1]]
        values$BQcolour <- items[[2]]
        values$BQcountF <- items[[3]]
        values$BQcountW <- items[[4]]
        values$BQcountP <- items[[5]]
        
        menuItem(
          text = "Base Quality",
          tabName = "BQ",
          badgeLabel = values$BQflag,
          badgeColor = values$BQcolour
        )
        
      }
      else{
        menuItem(text = "Base Quality", tabName = "BQ")
      }
    })
    
    
    output$SQflag <- renderMenu({
      if (!is.null(fastqcInput) | length(input$files) > 1) {
        flags <- getSummary(data())
        flags <- subset(flags, Category == "Per sequence quality scores")
        
        items <- .menuItemLogic(flags = flags)
        
        values$SQflag <- items[[1]]
        values$SQcolour <- items[[2]]
        values$SQcountF <- items[[3]]
        values$SQcountW <- items[[4]]
        values$SQcountP <- items[[5]]
        
        menuItem(
          text = "Sequence Quality",
          tabName = "SQ",
          badgeLabel = values$SQflag,
          badgeColor = values$SQcolour
        )
        
      }
      else{
        menuItem(text = "Sequence Quality", tabName = "SQ")
      }
    })
    
    output$SCflag <- renderMenu({
      if (!is.null(fastqcInput) | length(input$files) > 1) {
        flags <- getSummary(data())
        flags <- subset(flags, Category == "Per base sequence content")
        
        items <- .menuItemLogic(flags = flags)
        
        values$SCflag <- items[[1]]
        values$SCcolour <- items[[2]]
        values$SCcountF <- items[[3]]
        values$SCcountW <- items[[4]]
        values$SCcountP <- items[[5]]
        
        menuItem(
          text = "Sequence Content",
          tabName = "SC",
          badgeLabel = values$SCflag,
          badgeColor = values$SCcolour
        )
        
      }
      else{
        menuItem(text = "Sequence Content", tabName = "SC")
      }
    })
    
    output$GCflag <- renderMenu({
      if (!is.null(fastqcInput) | length(input$files) > 1) {
        flags <- getSummary(data())
        flags <- subset(flags, Category == "Per sequence GC content")
        
        items <- .menuItemLogic(flags = flags)
        
        values$GCflag <- items[[1]]
        values$GCcolour <- items[[2]]
        values$GCcountF <- items[[3]]
        values$GCcountW <- items[[4]]
        values$GCcountP <- items[[5]]
        
        menuItem(
          text = "GC Content",
          tabName = "GC",
          badgeLabel = values$GCflag,
          badgeColor = values$GCcolour
        )
        
      }
      else{
        menuItem(text = "GC Content", tabName = "GC")
      }
    })
    
    output$NCflag <- renderMenu({
      if (!is.null(fastqcInput) | length(input$files) > 1) {
        flags <- getSummary(data())
        flags <- subset(flags, Category == "Per base N content")
        
        items <- .menuItemLogic(flags = flags)
        
        values$NCflag <- items[[1]]
        values$NCcolour <- items[[2]]
        values$NCcountF <- items[[3]]
        values$NCcountW <- items[[4]]
        values$NCcountP <- items[[5]]
        
        menuItem(
          text = "N Content",
          tabName = "NC",
          badgeLabel = values$NCflag,
          badgeColor = values$NCcolour
        )
        
      }
      else{
        menuItem(text = "N Content", tabName = "NC")
      }
    })
    
    output$SLDflag <- renderMenu({
      if (!is.null(fastqcInput) | length(input$files) > 1) {
        flags <- getSummary(data())
        flags <- subset(flags, Category == "Sequence Length Distribution")
        
        items <- .menuItemLogic(flags = flags)
        
        values$SLDflag <- items[[1]]
        values$SLDcolour <- items[[2]]
        values$SLDcountF <- items[[3]]
        values$SLDcountW <- items[[4]]
        values$SLDcountP <- items[[5]]
        
        menuItem(
          text = "Sequence Length Distribution",
          tabName = "SLD",
          badgeLabel = values$SLDflag,
          badgeColor = values$SLDcolour
        )
        
      }
      else{
        menuItem(text = "Sequence Length Distribution", tabName = "SLD")
      }
    })
    
    output$SDLflag <- renderMenu({
      if (!is.null(fastqcInput) | length(input$files) > 1) {
        flags <- getSummary(data())
        flags <- subset(flags, Category == "Sequence Duplication Levels")
        
        items <- .menuItemLogic(flags = flags)
        
        values$SDLflag <- items[[1]]
        values$SDLcolour <- items[[2]]
        values$SDLcountF <- items[[3]]
        values$SDLcountW <- items[[4]]
        values$SDLcountP <- items[[5]]
        
        menuItem(
          text = "Sequence Duplication Levels",
          tabName = "SDL",
          badgeLabel = values$SDLflag,
          badgeColor = values$SDLcolour
        )
        
      }
      else{
        menuItem(text = "Sequence Duplicaiton Levels", tabName = "SDL")
      }
    })
    
    output$OSflag <- renderMenu({
      if (!is.null(fastqcInput) | length(input$files) > 1) {
        flags <- getSummary(data())
        flags <- subset(flags, Category == "Overrepresented sequences")
        
        items <- .menuItemLogic(flags = flags)
        
        values$OSflag <- items[[1]]
        values$OScolour <- items[[2]]
        values$OScountF <- items[[3]]
        values$OScountW <- items[[4]]
        values$OScountP <- items[[5]]
        
        menuItem(
          text = "Overrepresented Sequences",
          tabName = "OS",
          badgeLabel = values$OSflag,
          badgeColor = values$OScolour
        )
        
      }
      else{
        menuItem(text = "Overrepresented Sequences", tabName = "OS")
      }
    })
    
    
    output$ACflag <- renderMenu({
      if (!is.null(fastqcInput) | length(input$files) > 1) {
        flags <- getSummary(data())
        flags <- subset(flags, Category == "Adapter Content")
        
        items <- .menuItemLogic(flags = flags)
        
        values$ACflag <- items[[1]]
        values$ACcolour <- items[[2]]
        values$ACcountF <- items[[3]]
        values$ACcountW <- items[[4]]
        values$ACcountP <- items[[5]]
        
        menuItem(
          text = "Adapter Content",
          tabName = "AC",
          badgeLabel = values$ACflag,
          badgeColor = values$ACcolour
        )
        
      }
      else{
        menuItem(text = "Adapter Content", tabName = "AC")
      }
    })
    
    output$KCflag <- renderMenu({
      if (!is.null(fastqcInput) | length(input$files) > 1) {
        flags <- getSummary(data())
        flags <- subset(flags, Category == "Kmer Content")
        
        items <- .menuItemLogic(flags = flags)
        
        values$KCflag <- items[[1]]
        values$KCcolour <- items[[2]]
        values$KCcountF <- items[[3]]
        values$KCcountW <- items[[4]]
        values$KCcountP <- items[[5]]
        
        menuItem(
          text = "K-mer Content",
          tabName = "KC",
          badgeLabel = values$KCflag,
          badgeColor = values$KCcolour
        )
        
      }
      else{
        menuItem(text = "K-mer Content", tabName = "KC")
      }
    })
    
    observe({
      input$files
      
      output$"BQboxF" <- .renderValBox(
        count = values$BQcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
      
      output$SQboxF <- .renderValBox(
        count = values$SQcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
      
      output$SCboxF <- .renderValBox(
        count = values$SCcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
      
      output$GCboxF <- .renderValBox(
        count = values$GCcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
      
      output$NCboxF <- .renderValBox(
        count = values$NCcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
      
      output$SLDboxF <- .renderValBox(
        count = values$SLDcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
      
      output$SDLboxF <- .renderValBox(
        count = values$SDLcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
      
      output$OSboxF <- .renderValBox(
        count = values$OScountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
      
      output$ACboxF <- .renderValBox(
        count = values$ACcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
      
      output$KCboxF <- .renderValBox(
        count = values$KCcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
      
      #render warn value boxes
      
      output$BQboxW <- .renderValBox(
        count = values$BQcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
      
      output$SQboxW <- .renderValBox(
        count = values$SQcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
      
      output$SCboxW <- .renderValBox(
        count = values$SCcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
      
      output$GCboxW <- .renderValBox(
        count = values$GCcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
      
      output$NCboxW <- .renderValBox(
        count = values$NCcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
      
      output$SLDboxW <- .renderValBox(
        count = values$SLDcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
      
      output$SDLboxW <- .renderValBox(
        count = values$SDLcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
      
      output$OSboxW <- .renderValBox(
        count = values$OScountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
      
      output$ACboxW <- .renderValBox(
        count = values$ACcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
      
      output$KCboxW <- .renderValBox(
        count = values$KCcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
      
      #render Pass value boxes
      
      output$BQboxP <- .renderValBox(
        count = values$BQcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
      
      output$SQboxP <- .renderValBox(
        count = values$SQcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
      
      output$SCboxP <- .renderValBox(
        count = values$SCcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
      
      output$GCboxP <- .renderValBox(
        count = values$GCcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
      
      output$NCboxP <- .renderValBox(
        count = values$NCcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
      
      output$SLDboxP <- .renderValBox(
        count = values$SLDcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
      
      output$SDLboxP <- .renderValBox(
        count = values$SDLcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
      
      output$OSboxP <- .renderValBox(
        count = values$OScountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
      
      output$ACboxP <- .renderValBox(
        count = values$ACcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
      
      output$KCboxP <- .renderValBox(
        count = values$KCcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
    })
    #render fail value boxes
    
    
   
    
    # render plots
    
    ####################
    # Summary
    ####################
    
    output$SummaryFlags <- renderPlotly({
      if (!length(data())) {
        stop("Please load data to display plot.")
      }
      else{
          withProgress(
              min = 0,
              max = 1,
              value = 0.8,
              message = "Rendering plot.",
              {
        plotSummary(
          data(),
          usePlotly = TRUE,
          cluster = input$Sumcluster,
          dendrogram = input$Sumcluster
        ) %>%
          layout(margin = list(r = 200))}
          )
      }
    })
    
    ####################
    # Read Totals
    ####################
    
    output$ReadTotals <- renderPlotly({
      if (!length(data())) {
        stop("Please load data to display plot.")
      }
      else{
          withProgress(
              min = 0,
              max = 1,
              value = 0.8,
              message = "Rendering plot.",
              {
        plotReadTotals(data(),
                       usePlotly = TRUE,
                       duplicated = input$showDup) %>%
          layout(margin = list(l = 100, r = 200))}
          )
      }
    })
    
    ####################
    # Base Quality
    ####################
    
    output$baseQualHeatmap <- renderPlotly({
      if (!length(data())) {
        stop("Please load data to display plot.")
      }
      else{
          withProgress(
              min = 0,
              max = 1,
              value = 0.8,
              message = "Rendering plot.",
              {
        plotBaseQuals(
          data(),
          usePlotly = TRUE,
          plotType = "heatmap",
          plotValue = input$BQplotValue,
          cluster = input$BQcluster,
          dendrogram = input$BQcluster
        ) %>%
          layout(margin = list(r = 200))   
              }
          )
      }
    })
    
    
    output$BaseQualitiesSingle <- renderPlotly({
        if (!length(data())) {
            stop("Please load data to display plot.")
        }
        else{
            if(class(data()) == "character") names <- fqName(FastqcDataList(data()))
            else names <- fqName(data())
            
            if (is.null(event_data("plotly_click")$key[[1]])) num <- 1
            else {
                click <- event_data("plotly_click")
                if(!event_data("plotly_click")$key[[1]] %in% names) num <- 1
                else num <- which(names == click$key[[1]])
            }
            sub_fdl <- data()[[num]]
            plotBaseQuals(sub_fdl, usePlotly = TRUE,
                          plotValue = input$BQplotValue) %>%
                layout(margin = list(r = 200, b = 50))
        }
    })
    
    ####################
    # Sequence Quality
    ####################
    
    output$seqQualHeatmap <- renderPlotly({
      if (!length(data())) {
        stop("Please load data to display plot.")
      }
      else{withProgress(
          min = 0,
          max = 1,
          value = 0.8,
          message = "Rendering plot.",
          {
        plotSeqQuals(
          data(),
          cluster = input$SQcluster,
          counts = input$SQType == "Counts",
          dendrogram = input$SQcluster,
          usePlotly = TRUE
        ) %>% layout(margin = list(r = 200))}
      )
      }
    })
    
    output$SeqQualitiesSingle <- renderPlotly({
        if (!length(data())) {
            stop("Please load data to display plot.")
        }
        else{
            if(class(data()) == "character") names <- fqName(FastqcDataList(data()))
            else names <- fqName(data())
            
            if (is.null(event_data("plotly_click")$key[[1]])) num <- 1
            else {
                click <- event_data("plotly_click")
                if(!event_data("plotly_click")$key[[1]] %in% names) num <- 1
                else num <- which(names == click$key[[1]])
            }
        sub_fdl <- data()[[num]]
        qualPlot <-
          plotSeqQuals(sub_fdl, usePlotly = TRUE,
                       counts = input$SQType == "Counts") %>%
          layout(margin = list(r = 200),
                 legend = list(orientation = 'h', title = ""))
      }
    })
    
    ####################
    # Sequence Content
    ####################
    
    output$SCHeatmap <- renderPlotly({
      if (!length(data())) {
        stop("Please load data to display plot.")
      }
      else{withProgress(
          min = 0,
          max = 1,
          value = 0.8,
          message = "Rendering plot.",
          {
        plotSeqContent(
          data(),
          cluster = input$SCcluster,
          dendrogram = input$SCcluster,
          usePlotly = TRUE
        ) %>%
          layout(margin = list(r = 200))}
      )
      }
    })
    
    
    output$SCsingle <- renderPlotly({
        if (!length(data())) {
            stop("Please load data to display plot.")
        }
        else{
            if(class(data()) == "character") names <- fqName(FastqcDataList(data()))
            else names <- fqName(data())
            
            if (is.null(event_data("plotly_click")$key[[1]])) num <- 1
            else {
                click <- event_data("plotly_click")
                if(!event_data("plotly_click")$key[[1]] %in% names) num <- 1
                else num <- which(names == click$key[[1]])
            }
        sub_fdl <- data()[[num]]
        plotSeqContent(sub_fdl, usePlotly = TRUE) %>%
          layout(margin = list(r = 200))
      }
    })
    
    
    ####################
    # GC Content
    ####################
    
    output$theoreticalGC <- renderUI({
      if (input$theoreticalGC) {
        radioButtons(
          inputId = "theoreticalType",
          label = "What type of data?",
          choices = c("Genome", "Transcriptome"),
          selected = "Genome"
        )
      }
    })
    
    output$GCspecies <- renderUI({
      if (!is.null(input$theoreticalGC)) {
        if (input$theoreticalGC) {
          selectInput(
            "GCspecies",
            "Select species",
            choices = gcAvail(ngsReports::gcTheoretical, type = input$theoreticalType)$Name,
            selected = "Hsapiens")
          
        }
      }
    })
    
    
    
    output$GCheatmap <- renderPlotly({
      if (!length(data())) {
        stop("Please load data to display plot.")
      }
      else{
        if (is.null(input$GCspecies)) {
          GCspecies <- FALSE
        } else
          GCspecies <- input$GCspecies
        
        GCtype <- input$GCheatType == "Count"
        withProgress(
            min = 0,
            max = 1,
            value = 0.8,
            message = "Rendering plot.",
            {
        if (is.null(input$theoreticalGC)) {
          plotGcContent(
            data(),
            cluster = input$GCcluster,
            plotType = "heatmap",
            theoreticalType = input$theoreticalType,
            dendrogram = input$GCcluster,
            usePlotly = TRUE
          ) %>%
            layout(margin = list(r = 200))
        } else{
          plotGcContent(
            data(),
            cluster = input$GCcluster,
            plotType = "heatmap",
            theoreticalType = input$theoreticalType,
            theoreticalGC = input$theoreticalGC,
            dendrogram = input$GCcluster,
            species = GCspecies,
            usePlotly = TRUE
          ) %>%
            layout(margin = list(r = 200))
        }}
        )
        
      }
    })
    
    
   
    output$GCSingle <- renderPlotly({
      if (!length(data())) {
        stop("Please load data to display plot.")
      }
      else{
        if (is.null(input$GCspecies)) {
          GCspecies <- FALSE
        } else
          GCspecies <- input$GCspecies
        
        if(class(data()) == "character") names <- fqName(FastqcDataList(data()))
        else names <- fqName(data())
        
        if (is.null(event_data("plotly_click")$key[[1]])) num <- 1
        else {
            click <- event_data("plotly_click")
            if(!event_data("plotly_click")$key[[1]] %in% names) num <- 1
            else num <- which(names == click$key[[1]])
        }
        sub_fdl <- data()[[num]]
        if (is.null(input$theoreticalGC)) {
          GCSingle <- plotGcContent(sub_fdl, usePlotly = TRUE,
                                    counts = FALSE)
        } else{
          GCSingle <- plotGcContent(
            sub_fdl,
            usePlotly = TRUE,
            counts = FALSE,
            theoreticalGC = input$theoreticalGC,
            theoreticalType = input$theoreticalType,
            species = GCspecies
          )
        }
        GCSingle %>%
          layout(margin = list(r = 200),
                 legend = list(orientation = 'h', title = ""))
      }
    })
    
    
    ####################
    # N Content
    ####################
    
    output$NCheatmap <- renderPlotly({
      if (!length(data())) {
        stop("Please load data to display plot.")
      }
      else{
          withProgress(
              min = 0,
              max = 1,
              value = 0.8,
              message = "Rendering plot.",
              {
        plotNContent(
          data(),
          cluster = input$Ncluster,
          dendrogram = input$Ncluster,
          usePlotly = TRUE
        )}
          )
      }
    })
    
    
    
    
    # N Content single plot
    
    output$NCsingle <- renderPlotly({
        if (!length(data())) {
            stop("Please load data to display plot.")
        }
        else{
            if(class(data()) == "character") names <- fqName(FastqcDataList(data()))
            else names <- fqName(data())
            
            if (is.null(event_data("plotly_click")$key[[1]])) num <- 1
            else {
                click <- event_data("plotly_click")
                if(!event_data("plotly_click")$key[[1]] %in% names) num <- 1
                else num <- which(names == click$key[[1]])
            }
        sub_fdl <- data()[[num]]
        plotNContent(sub_fdl, usePlotly = TRUE) %>%
          layout(margin = list(r = 200, b = 50))
      }
    })
    
    
    ####################
    # Sequence Length Distribution
    ####################
    
    output$SLHeatmap <- renderPlotly({
      if (!length(data())) {
        stop("Please load data to display plot.")
      }
      else{withProgress(
          min = 0,
          max = 1,
          value = 0.8,
          message = "Rendering plot.",
          {
        plotSeqLengthDistn(
          data(),
          cluster = input$SLcluster,
          dendrogram = input$SLcluster,
          counts = input$SLType == "Counts",
          usePlotly = TRUE
        ) %>%
          layout(margin = list(r = 200))}
      )
      }
    })
    
    
    output$SLSingle <- renderPlotly({
        if (!length(data())) {
            stop("Please load data to display plot.")
        }
        else{
            if(class(data()) == "character") names <- fqName(FastqcDataList(data()))
            else names <- fqName(data())
            
            if (is.null(event_data("plotly_click")$key[[1]])) num <- 1
            else {
                click <- event_data("plotly_click")
                if(!event_data("plotly_click")$key[[1]] %in% names) num <- 1
                else num <- which(names == click$key[[1]])
            }
        sub_fdl <- data()[[num]]
        plotSeqLengthDistn(sub_fdl,
                           usePlotly = TRUE,
                           plotType = "line",
                           counts = input$SLType == "Counts") %>%
          layout(margin = list(r = 200))
      }
    })
    
    ####################
    # Sequence Duplicaiton Levels
    ####################
    
    output$DupHeatmap <- renderPlotly({
      if (!length(data())) {
        stop("Please load data to display plot.")
      }
      else{
          withProgress(
              min = 0,
              max = 1,
              value = 0.8,
              message = "Rendering plot.",
              {
        plotDupLevels(
          data(),
          cluster = input$Dupcluster,
          dendrogram = input$Dupcluster,
          usePlotly = TRUE
        ) %>%
          layout(margin = list(r = 200))}
          )
      }
    })
    
    output$DupSingle <- renderPlotly({
        if (!length(data())) {
            stop("Please load data to display plot.")
        }
        else{
            if(class(data()) == "character") names <- fqName(FastqcDataList(data()))
            else names <- fqName(data())
            
            if (is.null(event_data("plotly_click")$key[[1]])) num <- 1
            else {
                click <- event_data("plotly_click")
                if(!event_data("plotly_click")$key[[1]] %in% names) num <- 1
                else num <- which(names == click$key[[1]])
            }
        sub_fdl <- data()[[num]]
        plotDupLevels(sub_fdl, usePlotly = TRUE) %>%
          layout(margin = list(r = 200))
      }
    })
    
    ####################
    # Overrepresented sequences
    ####################
    
    output$OSummary <- renderPlotly({
      if (!length(data())) {
        stop("Please load data to display plot.")
      }
      else{
          withProgress(
              min = 0,
              max = 1,
              value = 0.8,
              message = "Rendering plot.",
              {
        plotOverrep(
          data(),
          usePlotly = TRUE,
          cluster = input$OScluster,
          dendrogram = input$OScluster
        ) %>%
          layout(margin = list(r = 200))}
          )
      }
      
    })
    
    
    output$OSsingle <- renderPlotly({
        if (!length(data())) {
            stop("Please load data to display plot.")
        }
        else{
            if(class(data()) == "character") names <- fqName(FastqcDataList(data()))
            else names <- fqName(data())
            
            if (is.null(event_data("plotly_click")$key[[1]])) num <- 1
            else {
                click <- event_data("plotly_click")
                if(!event_data("plotly_click")$key[[1]] %in% names) num <- 1
                else num <- which(names == click$key[[1]])
            }
        sub_fdl <- data()[[num]]
        plotOverrep(sub_fdl, usePlotly = TRUE) %>%
          layout(margin = list(r = 200))
      }
    })
    
    ####################
    # Adapter Content
    ####################
    
    
    output$ACheatmap <- renderPlotly({
      if (!length(data())) {
        stop("Please load data to display plot.")
      }
      else{withProgress(
          min = 0,
          max = 1,
          value = 0.8,
          message = "Rendering plot.",
          {
        ACplot <- plotAdapterContent(
          data(),
          adapterType = input$ACtype,
          usePlotly = TRUE,
          dendrogram = input$ACcluster,
          cluster = input$ACcluster
        )}
      )
        if (!is.null(ACplot))
          ACplot %>% layout(margin = list(r = 200))
        else
          stop(
            paste(
              "Sequences did not contain any",
              input$ACtype,
              "content, Please load another."
            )
          )
      }
    })
    
    output$ACsingle <- renderPlotly({
        if (!length(data())) {
            stop("Please load data to display plot.")
        }
        else{
            if(class(data()) == "character") names <- fqName(FastqcDataList(data()))
            else names <- fqName(data())
            
            if (is.null(event_data("plotly_click")$key[[1]])) num <- 1
            else {
                click <- event_data("plotly_click")
                if(!event_data("plotly_click")$key[[1]] %in% names) num <- 1
                else num <- which(names == click$key[[1]])
            }
        sub_fdl <- data()[[num]]
        ACsing <- plotAdapterContent(sub_fdl, usePlotly = TRUE)
        
        if (!is.null(ACsing))
          ACsing %>%
          layout(margin = list(r = 200),
                 legend = list(orientation = 'h', title = ""))
        else
          stop(
            paste(
              "Sequences did not contain any",
              input$ACtype,
              "content, Please load another."
            )
          )
      }
    })
    
    
    ####################
    # k-mer Content
    ####################
    
    output$Kheatmap <- renderPlotly({
      if (!length(data())) {
        stop("Please load data to display plot.")
      }
      else{
          Kplot <- NULL
          withProgress(
              min = 0,
              max = 1,
              value = 0.8,
              message = "Rendering plot.",
              {
        Kplot <- plotKmers(
          data(),
          usePlotly = TRUE,
          cluster = input$KMcluster,
          dendrogram = input$KMcluster
        )}
          )
        if (!is.null(Kplot))
          Kplot %>% layout(margin = list(r = 200))
        else
          stop(paste("Samples have no Kmer content"))
      }
    })
    
    output$Ksingle <- renderPlotly({
        if (!length(data())) {
            stop("Please load data to display plot.")
        }
        else{
            if(class(data()) == "character") names <- fqName(FastqcDataList(data()))
            else names <- fqName(data())
            
            if (is.null(event_data("plotly_click")$key[[1]])) num <- 1
            else {
                click <- event_data("plotly_click")
                if(!event_data("plotly_click")$key[[1]] %in% names) num <- 1
                else num <- which(names == click$key[[1]])
            }
        sub_fdl <- data()[[num]]
        Ksing <- plotKmers(sub_fdl, usePlotly = TRUE)
        
        if (!is.null(Ksing))
          Ksing %>%
          layout(margin = list(r = 200, b = 50))
        else
          stop(paste(
            "Library did not contain any identified
            Kmers please load another."
          ))
      }
    })
  }
  
  runApp(list(ui = ui, server = server), launch.browser = TRUE)
  
}

