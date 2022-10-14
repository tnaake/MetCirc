#' @name shinyCircos
#'
#' @title Interactive visualisation of similarity and navigation of MS/MS
#' features
#'
#' @description
#' Visualise the similarity of MS/MS features in a reactive 
#' context. See \code{Details} the vignette for further descriptions on how to 
#' use \code{shinyCircos}.
#'
#' @param similarityMatrix \code{matrix}, \code{similarityMatrix} 
#' contains pair-wise similarity coefficients which give information about the
#' similarity between MS/MS features
#' @param sps \code{Spectra}, \code{sps} will be used to display information 
#' about the selected feature and will store information of annotation
#' @param condition \code{character} vector, specifies which conditions/samples
#' are displayed
#' @param ... further arguments passed to \code{shinyCircos}, e.g.
#' \code{cexFeatureNames} to pass to \code{plotCircos} to set font size in
#' \code{plotCircos} of feature names
#'
#' @details
#' The function is based on the \code{shiny} and \code{circlize} package.
#' The user can choose interactively thresholds, type of links (between or
#' within groups), display information about MS/MS features, permanently select
#' MS/MS features and export selected precursors. The \code{Spectra} object
#' stores annotation information about the MS/MS features. Names of features
#' within the \code{similarityMatrix} have to be found as entries
#' in \code{Spectra}. \code{sps$name} are used as identifiers and
#' \code{colnames}/\code{rownames} from \code{similarityMatrix} are cleaved
#' by the group identifier (separated by \code{"_"}). Annotation information 
#' is taken from \code{spectra} from the columns "names", "information", 
#' "classes" and "adduct" in the slot \code{metadata} of \code{spectra}. 
#' After exiting the application, the annotation will be written to the 
#' respective columns in the slot \code{metadata}. If one or several of these 
#' columns is already present in \code{metadata}, the column(s) will be used as 
#' the source of annotation information.
#'
#' @return \code{character}, \code{shinyCircos} returns a \code{character} 
#' vector with the permanently selected precursors and an object with the 
#' \code{MSpectra} object containing the annotation.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples 
#' data("spectra", package = "MetCirc")
#' similarityMat <- Spectra::compareSpectra(sps_tissue[1:10],
#'     FUN = MsCoreUtils::ndotproduct, ppm = 10, m = 0.5, n = 2)
#' rownames(similarityMat) <- colnames(similarityMat) <- sps_tissue$name[1:10]
#' \dontrun{
#' shinyCircos(similarityMatrix = similarityMat,
#'     sps = sps_tissue, condition = c("SPL", "LIM", "ANT", "STY"))
#' }
#' 
#' @importFrom circlize circos.par circlize circos.initialize 
#' @importFrom circlize circos.trackPlotRegion circos.clear
#' @importFrom MsCoreUtils ndotproduct
#' @importFrom shiny fluidPage tags column fluidRow tabsetPanel tabPanel HTML
#' @importFrom shiny wellPanel radioButtons sliderInput uiOutput actionButton
#' @importFrom shiny checkboxInput plotOutput conditionalPanel textInput
#' @importFrom shiny verbatimTextOutput renderUI htmlOutput reactiveValues 
#' @importFrom shiny isolate actionButton reactive eventReactive observe 
#' @importFrom shiny renderPlot renderText stopApp runApp req updateSelectInput
#' 
#' @export
shinyCircos <- function(similarityMatrix, sps, condition, ...) {

    ## check if names, information, classes or adduct is present in
    ## @metadata and add respectively, these slots will be queried
    if (is.null(sps@metadata$names))
        sps@metadata[["names"]] <- "Unknown"

    if (is.null(sps@metadata$information))
        sps@metadata[["information"]] <- "Unknown"

    if (is.null(sps@metadata$classes))
        sps@metadata[["classes"]] <- "Unknown"

    if (is.null(sps@metadata$adduct))
        sps@metadata[["adduct"]] <- "Unknown"

    ## circlize parameters
    circlize::circos.par(gap.degree = 0, cell.padding = c(0, 0, 0, 0),
        track.margin = c(0.0, 0))

    ## order sps per condition according to rt, mz, cluster
    rt <- sps$rtime
    prec_mz <- sps$precursorMz
    groupname <- rownames(similarityMatrix)

    ## create plots and assign to objects by recordPlot
    ## rt
    RT <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMatrix,
        sps = sps, type = "retentionTime", condition = condition)
    link0dfRT <- RT[["link0df"]]
    rt_match <- RT[["type_match"]]

    ## get group and name from rt_match
    ## rt_match is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    groupRT <- lapply(strsplit(rt_match, split = "_"), "[", 1) |>
        unlist()
    nameRT <- lapply(strsplit(rt_match, split = "_"), function(x) x[length(x)]) |>
        unlist()

    ## plot filled and get degree features
    fillRT <- MetCirc:::recordPlotFill_degreeFeatures(rt_match, ...)
    degFeatRT <- fillRT[["degreeFeatures"]]
    fillRT <- fillRT[["plotFill"]]

    ## plot highlight
    highlightRT <- MetCirc:::recordPlotHighlight(rt_match, ...)

    ## mz
    MZ <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMatrix,
        sps = sps, type = "mz", condition = condition)
    link0dfMZ <- MZ[["link0df"]]
    mz_match <- MZ[["type_match"]]

    ## get group and name from mz_match
    ## mz_match is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    groupMZ <- lapply(strsplit(mz_match, split = "_"), "[", 1) |>
        unlist()
    nameMZ <- lapply(strsplit(mz_match, split = "_"), function(x) x[length(x)]) |>
        unlist()

    ## plot filled and get degree features
    fillMZ <- MetCirc:::recordPlotFill_degreeFeatures(mz_match, ...)
    degFeatMZ <- fillMZ[["degreeFeatures"]]
    fillMZ <- fillMZ[["plotFill"]]

    ## plot highlight
    highlightMZ <- MetCirc:::recordPlotHighlight(mz_match, ...)

    ## clustering
    Clust <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMatrix,
        sps = sps, type = "clustering", condition = condition)
    link0dfClust <- Clust[["link0df"]]
    clust_match <- Clust[["type_match"]]

    ## get group and name from clust_match
    ## clust_match is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    groupClust <- lapply(strsplit(clust_match, split = "_"), "[", 1) |>
        unlist()
    nameClust <- lapply(strsplit(clust_match, split = "_"), function(x) x[length(x)]) |>
        unlist()

    ## plot filled and get degree features
    fillClust <- MetCirc:::recordPlotFill_degreeFeatures(clust_match, ...)
    degFeatClust <- fillClust[["degreeFeatures"]]
    fillClust <- fillClust[["plotFill"]]

    ## plot highlight
    highlightClust <- MetCirc:::recordPlotHighlight(clust_match, ...)

    ## create list with all plots
    plot_l <- list(highlightMz = highlightMZ, fillMz = fillMZ,
        highlightRT = highlightRT, fillRT = fillRT,
        highlightClust = highlightClust, fillClust = fillClust)

    ui <- shiny::fluidPage(
        shiny::tags$head(shiny::tags$script('
            $(document).on("shiny:connected", function(e) {
            Shiny.onInputChange("innerWidth", window.innerWidth);
            });
            $(window).resize(function(e) {
            Shiny.onInputChange("innerWidth", window.innerWidth);
            });
            '
        )),
        shiny::column(4,
            shiny::fluidRow(
                shiny::tabsetPanel(id = "tabs",
                    shiny::tabPanel("Main", id = "Main", shiny::wellPanel(
                        shiny::radioButtons("choiceLinks", "choose type of links",
                            choices = c("all" = "all", "inter-class links" = "inter",
                                "intra-class links" = "intra"), selected = "all"),
                        shiny::sliderInput("threshold",
                            "Threshold for similarity to display",
                                   min = 0, max = 1, value = c(0.8, 1)),
                        shiny::radioButtons("order", "order within groups",
                            choices = c("clustering" = "clustering", "m/z" = "mz",
                              "retention time" = "retentionTime"), selected = "mz"),
                        shiny::uiOutput("annotationName"),
                        shiny::uiOutput("annotationClass"),
                        shiny::uiOutput("annotationInformation"),
                        shiny::uiOutput("annotationAdduct"),
                        shiny::uiOutput("annotationButton"),
                        shiny::actionButton("resetClickIndices", "Unselect features"),
                        shiny::actionButton("stop", "Stop and export \n selected features")
                    )),
                    shiny::tabPanel("Spectra", id = "Spectra", shiny::wellPanel(
                        shiny::selectInput(inputId = "subjectSpectra", 
                            label = "Choose MS/MS spectra #1:",
                            choices = sps$name, selected = sps$name[1]),
                        shiny::selectInput(inputId = "querySpectra", 
                            label = "Choose MS/MS spectra #2:",
                            choices = sps$name, selected = sps$name[1])
                    )),
                    shiny::tabPanel("Appearance", id = "Appearance", shiny::wellPanel(
                        shiny::sliderInput("plotSize", "plot size",
                           min = 0.5, max = 1.5, value = 1),
                        shiny::sliderInput("precision", "precision of numbers",
                           value = 2, min = 0, max = 5, step = 1),
                        shiny::checkboxInput("legend", "legend", value = FALSE)
                    ))
                )
            ),
            shiny::plotOutput("circosLegend", height = "300")
        ),
        shiny::column(8,
            shiny::conditionalPanel(condition = "input.tabs  == 'Spectra'",
                shiny::plotOutput("plotSpec")),
            shiny::fluidRow(shiny::plotOutput("circos",
                dblclick = "circosDbl", click = "circosSgl", 
                width = "80%", height = "80%")),
            shiny::fluidRow(
                shiny::verbatimTextOutput("dimension_display"),
                shiny::textOutput("sglConnectedFeature"),
                shiny::verbatimTextOutput("dblClickFeature")
            )
        )
    )

    server <- function(input, output, session) {

        ## reactiveValues for sps, stores annotation information
        spe <- shiny::reactiveValues(sps = sps)
        
        ## single click: which is the current sector?
        indSgl <- shiny::reactiveValues(ind = NULL)
        shiny::observe({
            if (!is.null(input$circosSgl$x))
                indSgl$ind <- minFragCart2Polar(input$circosSgl$x,
                    input$circosSgl$y, degFeat())
        })
        
        ## reactiveValues for click coordinates
        coordsNewSgl <- shiny::reactiveValues(X = 0, Y = 0)
        coordsOldSgl <- shiny::reactiveValues(X = 0, Y = 0)

        shiny::observe({
            if (!is.null(input$circosSgl$x)) {
                coordsNewSgl$X <- input$circosSgl$x
                coordsNewSgl$Y <- input$circosSgl$y
                coordsOldSgl$X <- coordsNewSgl$X
                coordsOldSgl$Y <- coordsNewSgl$Y
            } else {
                coordsNewSgl$X <- coordsOldSgl$X
                coordsNewSgl$Y <- coordsOldSgl$Y
            }
        })

        ## double click: which is the current sector?
        coordsNewDbl <- shiny::reactiveValues(X = 0, Y = 0)
        coordsOldDbl <- shiny::reactiveValues(X = 0, Y = 0)

        shiny::observe({
            if (!is.null(input$circosDbl$x)) {
                coordsNewDbl$X <- input$circosDbl$x
                coordsNewDbl$Y <- input$circosDbl$y
                coordsOldDbl$X <- coordsNewDbl$X
                coordsOldDbl$Y <- coordsNewDbl$Y
            } else {
                coordsNewDbl$X <- coordsOldDbl$X
                coordsNewDbl$Y <- coordsOldDbl$Y
            }
        })
        
        ## reactive value which stores double clicked indices (inds = storage,
        ## new = new indices)
        indDbl <- shiny::reactiveValues(ind = NULL, new = NULL)
        shiny::observe({
            if (!is.null(input$circosDbl$x)) {

                minInd <- minFragCart2Polar(input$circosDbl$x,
                        input$circosDbl$y, degFeat())
                if (!is.na(minInd)) {
                    GNselect <- GN()[minInd]
                    selected <- strsplit(GNselect, split = "_")[[1]]
                    groupSelected <- selected[1]
                    nameSelected <- selected[2:length(selected)]
                    newNG <- paste(groupSelected, nameSelected, sep = "_")
                    ## write truncated name to indDbl$new
                    indDbl$new <- newNG
                } else  indDbl$new <- NULL
            }
        })


        ## write double-clicked (truncated) names to indDblMZ, indDblRT,
        ## indDblCluster
        indDblMZ <- shiny::reactiveValues(ind = NULL)
        shiny::observe({
            if (!is.null(input$circosDbl$x)) {
                if (!is.null(indDbl$new)) {
                    newMZ <- paste(groupMZ, nameMZ, sep = "_")
                    newIndMZ <- match(indDbl$new, newMZ)
                    if (shiny::isolate(newIndMZ %in% indDblMZ$ind)) {
                        indDblMZ$ind <- shiny::isolate(indDblMZ$ind[-which(newIndMZ == indDblMZ$ind)])
                    } else {
                        indDblMZ$ind <- shiny::isolate(c(indDblMZ$ind, newIndMZ))
                    }
                }
            }
        })

        indDblRT <- shiny::reactiveValues(ind = NULL)
        shiny::observe({
            if (!is.null(input$circosDbl$x)) {
                if (!is.null(indDbl$new)) {
                    newRT <- paste(groupRT, nameRT, sep = "_")
                    newIndRT <- match(indDbl$new, newRT)
                    if (shiny::isolate(newIndRT %in% indDblRT$ind)) {
                        indDblRT$ind <- isolate(indDblRT$ind[-which(newIndRT == indDblRT$ind)])
                    } else {
                        indDblRT$ind <- shiny::isolate(c(indDblRT$ind, newIndRT))
                    }
                }
            }
        })

        indDblCluster <- shiny::reactiveValues(ind = NULL)
        shiny::observe({
            if (!is.null(input$circosDbl$x)) {
                if(!is.null(indDbl$new)) {
                    newCl <- paste(groupClust, nameClust, sep = "_")
                    newIndCl <- match(indDbl$new, newCl)
                    if (isolate(newIndCl %in% indDblCluster$ind)) {
                        indDblCluster$ind <- shiny::isolate(indDblCluster$ind[-which(newIndCl == indDblCluster$ind)])
                    } else {
                        indDblCluster$ind <- shiny::isolate(c(indDblCluster$ind, newIndCl))
                    }
                }
            }
        })
        
        shiny::observe({
            input$resetClickIndices
            shiny::isolate(indDblMZ$ind <- NULL)
            shiny::isolate(indDblRT$ind <- NULL)
            shiny::isolate(indDblCluster$ind <- NULL)
            shiny::isolate(indDbl$new <- NULL)
            shiny::isolate(onCircle$is <- FALSE)
            shiny::isolate(indSgl$ind <- NULL)
        })

        ## reset indSgl when changing radio button order
        shiny::observe({
            input$order
            shiny::isolate(onCircle$is <- FALSE)
            shiny::isolate(indSgl$ind <- NULL)
        })

        
        ## ordering of features, use predefined groupname object
        GN <- shiny::reactive({
            MetCirc:::select(input$order, mz_match, rt_match, clust_match)
        })

        ## get degree of features
        degFeat <- shiny::reactive({
            MetCirc:::select(input$order, degFeatMZ, degFeatRT, degFeatClust)
        })
        
        ind <- reactive({
            if (length(indSgl$ind) > 0) {
                
                nameSgl <- GN()[indSgl$ind] |>
                    strsplit(split = "_") |>
                    lapply(function(i) i[[-1]]) |>
                    lapply(function(x) paste(x, collapse = "_")) |>
                    unlist()
                which(nameSgl == sps$name)
            }
        })
        
        ## is mouse over the track 1?
        onCircle <- shiny::reactiveValues(is = NULL)
        shiny::observe({
            if (!is.null(coordsNewSgl$X)) {
                .dist <- sqrt(coordsOldSgl$X^2 + coordsOldSgl$Y^2)
                if (.dist >= 0.8 & .dist <= 1) {
                    onCircle$is <- TRUE
                } else {
                    onCircle$is <- FALSE
                }
            } else onCircle$is <- FALSE
        })

        ## annotation
        output$annotationName <- shiny::renderUI({
            if (length(indSgl$ind) > 0 && onCircle$is) {
                shiny::textInput("names", label = "name",
                    value = shiny::isolate(spe$sps@metadata$names[ind()]))
            } else NULL
        })

        output$annotationClass <- shiny::renderUI({
            if (length(indSgl$ind) > 0 && onCircle$is) {
                shiny::textInput("classes", label = "class",
                    value = shiny::isolate(spe$sps@metadata$classes[ind()]))
            } else NULL
        })

        output$annotationInformation <- shiny::renderUI({
            if (length(indSgl$ind) > 0 && onCircle$is) {
                shiny::textInput("information", label = "information",
                    value = shiny::isolate(spe$sps@metadata$information[ind()]))
            } else NULL
        })

        output$annotationAdduct <- shiny::renderUI({
            if (length(indSgl$ind) > 0 && onCircle$is) {
                shiny::textInput("adduct", label = "adduct",
                    value = shiny::isolate(spe$sps@metadata$adduct[ind()]))
            } else NULL
        })

        output$annotationButton <- shiny::renderUI({
            if (length(indSgl$ind) > 0 && onCircle$is)
                shiny::actionButton(inputId = "annotate", label = "update annotation")
        })

        ## eventReactive for annotation ind
        indAnn <- shiny::eventReactive(input$annotate, {ind()})

        ## eventReactive for input$name (records "character" for annotation)
        annotateNames <- shiny::eventReactive(input$annotate, {
            input$names
        })

        ## eventReactive for input$classes
        annotateClasses <- shiny::eventReactive(input$annotate, {
            input$classes
        })

        ## eventReactive for input$information
        annotateInformation <- shiny::eventReactive(input$annotate, {
            input$information
        })

        ## eventReactive for input$adduct
        annotateAdduct <- shiny::eventReactive(input$annotate, {
            input$adduct
        })

        ## observe annotations and write to respective columns in slot
        ## elementMetadata
        shiny::observe({
            spe$sps@metadata$names[indAnn()] <- annotateNames()
        })
        shiny::observe({
            spe$sps@metadata$classes[indAnn()] <- annotateClasses()
        })
        shiny::observe({
            spe$sps@metadata$information[indAnn()] <- annotateInformation()
        })
        shiny::observe({
            spe$sps@metadata$adduct[indAnn()] <- annotateAdduct()
        })
        ## end annotation


        ## get link0df
        link0df <- shiny::reactive({
            MetCirc:::select(input$order, link0dfMZ, link0dfRT, link0dfClust)
        })

        ## create reactive expression for linkDf which is cut according to
        ## set radioButton (input$choiceLinks)
        linkDf_cut <- shiny::reactive(
            cutLinkDf(link0df(), type = input$choiceLinks))

        ## threshold linkDf_cut
        linkDf_threshold <- shiny::reactive(
            thresholdLinkDf(linkDf_cut(),
                input$threshold[1], input$threshold[2]))

        ## plotting
        initializePlot <- shiny::reactive({
            circlize::circos.initialize(factor(GN(), levels = GN()),
                xlim = matrix(rep(c(0,1), length(mz_match)), ncol = 2, byrow = TRUE))
            circlize::circos.trackPlotRegion(factor(GN(), levels = GN()), ylim = c(0,1))
        })

        ## assign to output$circos: actual plotting
        output$circos <- shiny::renderPlot({
            indDblMZ$ind
            initializePlot()
            MetCirc:::replayPlotOrder(orderMatch = input$order,
                onCircle = onCircle$is, plot_l = plot_l, ind = indDblMZ$ind)
            MetCirc:::replayPlotAdd(orderMatch = input$order,
                onCircle = onCircle$is, linkDf = linkDf_threshold(),
                mz_match = mz_match, rt_match = rt_match,
                clust_match = clust_match, ind = indSgl$ind,
                indMz = indDblMZ$ind, indRT = indDblRT$ind,
                indCluster = indDblCluster$ind)
        }, 
            width = shiny::reactive(
                ifelse(is.null(input$innerWidth), "50%", 
                    input$innerWidth * 0.6)),
            height = shiny::reactive(
                ifelse(is.null(input$innerWidth), "50%", 
                    input$innerWidth * 0.6))
        )

        output$circosLegend <- shiny::renderPlot({
            if (!is.null(input$legend)) if(input$legend)
                circosLegend(rt_match, highlight = TRUE)
        })

        ## show when Clicking the feature which connects to it
        linkDfIndsSgl <- shiny::reactive({
            getLinkDfIndices(GN()[indSgl$ind], linkDf_threshold())
        })

        output$sglConnectedFeature <- shiny::renderText({
            if (!is.null(onCircle$is) & onCircle$is) {
                    printInformationSelect(select = GN()[indSgl$ind],
                        sps = spe$sps, linkDfInd = linkDfIndsSgl(),
                        linkDf = linkDf_threshold(),
                        similarityMatrix = similarityMatrix,
                        roundDigits = input$precision)
            }
        })

        output$dblClickFeature <- shiny::renderText({
            if (length(indDblMZ$ind) > 0)
                c("(permanently) selected features: ", mz_match[indDblMZ$ind])
            else "no features (permanently) selected"
        })

        ## Tab Spectra
        shiny::observe({
            shiny::updateSelectInput(session = session,
                inputId = "subjectSpectra", choices = GN(),
                selected = GN()[indSgl$ind])
        })
        shiny::observe({
            shiny::req(GN())
            shiny::updateSelectInput(session = session,
                inputId = "querySpectra", choices = GN())
        })

        output$plotSpec <- shiny::renderPlot({
            if (!is.null(input$subject) && !is.null(input$query)) {
                plotSpectra(sps, subject = input$subject, query = input$query)
            }
        })

        ## on exit
        shiny::observe({
            if (input$stop == 0) {return()
            } else {
                circlize::circos.clear()
                selFeat <- as.character(paste(
                    groupMZ[indDblMZ$ind], nameMZ[indDblMZ$ind],
                    sep = "_"))
                shiny::stopApp(list(sps = spe$sps, selectedFeatures = selFeat))
            }
        })

    }

    app <- list(ui = ui, server = server)
    shiny::runApp(app)
}

#' @name printInformationSelect
#' 
#' @title Display information on connected features of selected features
#' 
#' @description
#' Displays information on connected features of selected features.
#'
#' @param select \code{character}, obtained from \code{groupname}, 
#' \code{character} of selected feature
#' @param sps \code{Spectra} object containing spectra that are compared in 
#' \code{similarityMatrix}
#' @param linkDfInd \code{numeric} indices of selected features
#' @param linkDf \code{data.frame} that contains information of linked features 
#' for given thresholds
#' @param similarityMatrix \code{matrix} that is used to get information on the 
#' degree of similarity, \code{similarityMatrix} is an ordered version of a 
#' similarity matrix, see \code{?orderSimilarityMatrix}
#' @param roundDigits \code{numeric(1)},  how many digits should be displayed?
#'
#' @details
#' \code{printInformationSelect} is for internal use.
#'
#' @return
#' \code{character} that is in HTML format
#'
#' @examples
#' data("spectra", package = "MetCirc")
#' similarityMat <- Spectra::compareSpectra(sps_tissue[1:10],
#'     FUN = MsCoreUtils::ndotproduct, ppm = 10, m = 0.5, n = 2)
#' rownames(similarityMat) <- colnames(similarityMat) <- sps_tissue$name[1:10]
#' linkDf <- createLinkDf(similarityMatrix = similarityMat,
#'     sps = sps_tissue[1:10], 
#'     condition = c("SPL", "LIM", "ANT", "STY"), lower = 0.5, upper = 1)
#'     
#' ## cut link data.frame (here: only display links between groups)
#' linkDf_cut <- cutLinkDf(linkDf, type = "inter")
#' groupname <- c(as.character(linkDf_cut[, "spectrum1"]),
#'             as.character(linkDf_cut[, "spectrum2"]))
#' groupname <- unique(groupname)
#' 
#' ## arbitrarily select a feature
#' ind <- 2
#' linkDfInds <- getLinkDfIndices(groupname[ind], linkDf_cut)
#' MetCirc:::printInformationSelect(select = groupname[ind], 
#'     sps = sps_tissue[1:10], linkDfInd = linkDfInds, 
#'     linkDf = linkDf_cut, similarityMatrix = similarityMat)
#' 
#' @importFrom MsCoreUtils ndotproduct
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
printInformationSelect <- function(select, sps = NULL, 
    linkDfInd, linkDf, similarityMatrix, roundDigits = 2) {

    ## connected features: find
    connect <- linkDf[linkDfInd, c("spectrum1", "spectrum2")] |>
        unlist() |>
        unique() |>
        as.character()
    
    ## remove duplicated hovFeat in connect
    connect <- connect[connect != select]
    
    ## truncate group
    connect_cut <- lapply(strsplit(connect, split = "_"), "[", -1) |>
        lapply(function(x) paste(x, collapse = "_")) |> 
        unlist()
    select_cut <- lapply(strsplit(select, split = "_"), "[", -1) |>
        lapply(function(x) paste(x, collapse = "_")) |>
        unlist()
    
    ## truncate sps for select and connecting features
    select_sps <- sps[sps$name %in% select_cut]
    connect_sps <- sps[sps$name %in% connect_cut]

    ## define columns to collapse
    cols <- c("names", "information", "classes", "adduct")
    
    if (length(connect) == 0) {
        
        res <- paste0(select, " (",
            round(select_sps$precursorMz, roundDigits), ", ",
            round(select_sps$rtime, roundDigits), ", ",
            paste(select_sps@metadata[select_cut, cols], collapse = ", "),
            ") does not connect to any feature ")
        
    } else { ## if length > 0

        connChar <- character()
        degreeSimilarity <- similarityMatrix[select_cut, ]
  
        for (i in seq_along(connect)) {

            connect_i <- connect_cut[i]
            degreeSimilarityI <- round(degreeSimilarity[[connect_i]], 3)
            sps_i <- connect_sps[connect_sps$name == connect_i]

            connChar <- c(connChar, 
                paste0(connect[i], " (", degreeSimilarityI, ", ",
                    round(sps_i$precursorMz, roundDigits), ", ",
                    round(sps_i$rtime, roundDigits), ", ",
                    paste(connect_sps@metadata[connect_i, cols], collapse = ", "), 
                    ")", "<br/>")
            )
        }
 
        connChar <- paste(connChar, collapse = " ")

        res <- paste0(select, " (", 
            round(select_sps$precursorMz, roundDigits), ", ",
            round(select_sps$rtime, roundDigits), ", ",
            paste(select_sps@metadata[select_cut, cols], collapse = ", "),
            ") connects to ", " <br/>", connChar)
    }
    
    ## return 
    res
}

#' @name replayPlotOrder
#' 
#' @title Wrapper for \code{replayPlot}
#' 
#' @description
#' \code{replayPlotOrder} will call \code{replayPlot} from \code{grDevices} with
#' a \code{recordedplot} object based on \code{orderMatch}.
#'
#' @param orderMatch \code{character}, either \code{"mz"}, 
#' \code{"retentionTime"} or \code{"clustering"}
#' @param plot_l \code{list} with plots
#' @param onCircle \code{logical}, are coordinates on circle. If \code{FALSE} 
#' and no features are selected (\code{length(ind) == 0}), then filled plots 
#' are replayed, otherwise highlighted plots are replayed.
#' @param ind \code{numeric}, indices of clicked features
#'
#' @details
#' Helper function for \code{shinyCircos}.
#'
#' @return \code{replayedplot}
#'
#' @examples 
#' type_match <- c("a_1", "a_2", "a_3", "b_1", "b_2", "b_3", "c_1", "c_2")
#' plotCircos(type_match, NULL, initialize = TRUE, featureNames = TRUE,
#'     groupSector = TRUE, groupName = FALSE, links = FALSE,
#'     highlight = TRUE)
#' p <- recordPlot()
#' plot.new()
#' plot_l <- list(highlightMz = p)
#' MetCirc:::replayPlotOrder(orderMatch = "mz", onCircle = TRUE,
#'     plot_l = plot_l, ind = NULL)
#'
#' @importFrom grDevices replayPlot
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
replayPlotOrder <- function(orderMatch = "mz", onCircle = FALSE, plot_l, ind) {

    if (!onCircle & length(ind) == 0) {

        ## get plot based on orderMatch 
        p <- MetCirc:::select(orderMatch, mz = plot_l[["fillMz"]],
            rt = plot_l[["fillRT"]], clust = plot_l[["fillClust"]])
        
        ## plot
        grDevices::replayPlot(p)
    } else {

        ## get plot based on orderMatch
        p <- MetCirc:::select(orderMatch, mz = plot_l[["highlightMz"]],
            rt = plot_l[["highlightRT"]], clust = plot_l[["highlightClust"]])

        ## plot
        grDevices::replayPlot(p)
    }
}

#' @name replayPlotAdd
#'
#' @title Plot plotCircos or highlight
#'
#' @description
#' \code{replayPlotAdd} plots additional plots on a plot, either
#' plots \code{plotCircos} or \code{highlight}.
#'
#' @param orderMatch \code{character(1)}, either \code{"mz"}, 
#' \code{"retentionTime"}, or \code{"clustering"}
#' @param onCircle \code{logical}, are coordinates on circle. If \code{FALSE} 
#' and no features are selected (\code{length(ind) == 0}), then filled plots are 
#' replayed, otherwise highlighted plots are replayed.
#' @param linkDf \code{data.frame} that contains information of linked 
#' features for given thresholds
#' @param mz_match \code{character}, ordered vector according to m/z
#' @param rt_match \code{character}, ordered vector according to retention time
#' @param clust_match \code{character}, ordered vector according to clustering
#' @param ind \code{numeric}, indices of clicked features
#' @param indMz \code{numeric}, indices of clicked features for \code{"mz"}
#' ordering
#' @param indRT \code{numeric}, indices of clicked features for 
#' \code{"retentionTime"} ordering
#' @param indCluster \code{numeric}, indices of clicked features for 
#' \code{"clustering"} ordering
#'
#' @details
#' Helper function for \code{shinyCircos}.
#'
#' @return
#' Depending on \code{onCircle} and \code{indMz} either returns 
#' \code{plotCircos} or \code{highlight}
#'
#' @examples
#' data("spectra", package = "MetCirc")
#' similarityMat <- Spectra::compareSpectra(sps_tissue[1:10],
#'     FUN = MsCoreUtils::ndotproduct, ppm = 10, m = 0.5, n = 2)
#' rownames(similarityMat) <- colnames(similarityMat) <- sps_tissue$name[1:10]
#'   
#' ## order according to m/z 
#' mz_match <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMat, 
#'     sps = sps_tissue, type = "mz", 
#'     condition = c("SPL", "LIM", "ANT", "STY"))
#' linkDf <- mz_match[["link0df"]]
#' mz_match <- mz_match[["type_match"]]
#' 
#' ## order according to retention time 
#' rt_match <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMat, 
#'     sps = sps_tissue, type = "retentionTime", 
#'     condition = c("SPL", "LIM", "ANT", "STY"))
#' rt_match <- rt_match[["type_match"]]
#'
#' ## order according to clustering
#' clust_match <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMat, 
#'     sps = sps_tissue, type = "clustering", 
#'     condition = c("SPL", "LIM", "ANT", "STY"))
#' clust_match <- clust_match[["type_match"]]
#' circos.initialize(mz_match,##, levels  =  mz_match),
#'     xlim = matrix(rep(c(0,1), length(mz_match)), ncol = 2, byrow = TRUE))
#' #circos.trackPlotRegion(factor(mz_match, levels = mz_match), ylim = c(0,1))  
#' MetCirc:::replayPlotAdd(orderMatch = "mz", onCircle = FALSE, linkDf = linkDf, 
#'     mz_match = mz_match, rt_match = rt_match, clust_match = clust_match, 
#'     ind = 1, indMz = NULL, indRT = NULL, indCluster = NULL)
#'
#' @importFrom MsCoreUtils ndotproduct
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
replayPlotAdd <- function(orderMatch = "mz", onCircle = FALSE, linkDf,
        mz_match, rt_match, clust_match, ind, indMz, indRT, indCluster) {

    ## get type_match based on orderMatch
    type_match <- MetCirc:::select(orderMatch, mz_match, rt_match, clust_match)
    inds <- MetCirc:::select(orderMatch, indMz, indRT, indCluster)

    if (!onCircle & length(inds) == 0) {
        ## plot
        plotCircos(type_match, linkDf, initialize = FALSE,
            featureNames = FALSE, groupSector = FALSE, groupName = FALSE,
            links = TRUE, highlight = FALSE)

    } else {
        ## get inds
        inds <- if (onCircle) c(ind, inds) else inds

        ## plot
        highlight(type_match, inds, linkDf)
    }
}


#' @name recordPlotFill_degreeFeatures
#'
#' @title Record a plot of filled features and the degree of features
#'
#' @description
#' \code{recordPlotFill_degreeFeatures} records a plot of filled 
#' features and returns the degree of features.
#'
#' @param type_match \code{character}, ordered vector according to type
#' @param ... further arguments passed to \code{plotCircos}
#'
#' @details
#' Helper function for \code{shinyCircos}.
#'
#' @return
#' \code{list} of length 2, entry \code{plotFill} is of \code{recordedplot} and 
#' entry \code{degreeFeatures} is a \code{list} of vectors of \code{numeric(1)}
#'
#' @examples
#' type_match <- c("a_1", "a_2", "a_3", "b_1", "b_2", "b_3", "c_1", "c_2")
#' MetCirc:::recordPlotFill_degreeFeatures(type_match)
#'
#' @importFrom grDevices recordPlot
#' @importFrom graphics plot.new
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
recordPlotFill_degreeFeatures <- function(type_match, ...) {
    
    plotCircos(type_match, NULL, initialize = TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = FALSE, ...)

    fill <- grDevices::recordPlot()

    ## get degree of features
    degree <- lapply(type_match, function(x) {
        mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")])
    })

    graphics::plot.new()
    
    ## return the list
    list(plotFill = fill, degreeFeatures = degree)

}

#' @name recordPlotHighlight
#'
#' @title Return a \code{recordedplot} of \code{plotCircos} plot with 
#' \code{highlight = TRUE}
#'
#' @description
#' \code{recordPlotHighlight} returns a \code{recordedplot} object of 
#' \code{plotCircos} with \code{highlight = TRUE}
#'
#' @param type_match \code{character}, ordered vector according to type
#' @param ... further arguments passed to \code{plotCircos}
#'
#' @details
#' Helper function for \code{shinyCircos}.
#'
#' @return \code{recordedplot}
#'
#' @examples 
#' type_match <- c("a_1", "a_2", "a_3", "b_1", "b_2", "b_3", "c_1", "c_2")
#' MetCirc:::recordPlotHighlight(type_match)
#'
#' @importFrom graphics plot.new
#' @importFrom grDevices recordPlot
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
recordPlotHighlight <- function(type_match, ...) {
    
    ## use plotCircos and record the plt
    plotCircos(type_match, NULL, initialize = TRUE, featureNames = TRUE,
               groupSector = TRUE, groupName = FALSE, links = FALSE,
               highlight = TRUE, ...)
    highlight <- grDevices::recordPlot()
    graphics::plot.new()
    
    ## return
    highlight
}

#' @name typeMatch_link0
#'
#' @title Get typeMatch and link0 data frame
#'
#' @description
#' \code{typeMatch_link0} returns a list with accessors \code{"link0df"} and
#' \code{"type_match"}
#'
#' @param similarityMatrix \code{matrix} with pairwise similarity values
#' @param sps \code{Spectra} object
#' @param type \code{character(1)}, either \code{"mz"}, \code{"retentionTime"}, 
#' \code{"clustering"}
#' @param condition \code{character}
#'
#' @details Helper function for \code{shinyCircos}. 
#'
#' @return
#' \code{list} of length 2,  entry \code{link0df} is a \code{data.frame} and 
#' entry \code{type_match} is a \code{character} vector
#'
#' @examples 
#' data("spectra", package = "MetCirc")
#' similarityMat <- Spectra::compareSpectra(sps_tissue[1:10],
#'     FUN = MsCoreUtils::ndotproduct, ppm = 10, m = 0.5, n = 2)
#' rownames(similarityMat) <- colnames(similarityMat) <- sps_tissue$name[1:10]
#' 
#' ## order according to retention time 
#' MetCirc:::typeMatch_link0(similarityMatrix = similarityMat, 
#'     sps = sps_tissue, type = "mz", 
#'     condition = c("SPL", "LIM", "ANT", "STY"))
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
typeMatch_link0 <- function(similarityMatrix, sps, type, condition) {

    simMat <- orderSimilarityMatrix(similarityMatrix, sps = sps, type = type, 
        group = FALSE)

    ## get names of sps per condition
    inds <- MetCirc:::spectraCondition(sps, condition)

    link0df <- createLink0df(simMat, sps = sps, condition = condition)
    groupname <- rownames(simMat)
    
    type_match <- lapply(inds, function(inds_i) {
        type_match <- match(groupname, inds_i)
        type_match <- type_match[!is.na(type_match)]; inds_i[type_match]
    })
    type_match <- lapply(seq_along(type_match), function(i) {
        if (length(type_match[[i]]) > 0) {
            paste(condition[i], type_match[[i]], sep = "_")
        } else character()
    })
    type_match <- unique(unlist(type_match))

    ## return list
    list(link0df = link0df, type_match = type_match)
}

#' @name select
#'
#' @title Select variable based on condition
#'
#' @description
#' \code{select} returns \code{mz}, \code{rt} or \code{clust} depending on
#' \code{condition}.
#'
#' @param condition \code{character(1)}, either \code{"mz"},
#' \code{"retentionTime"}, or \code{"clustering"}
#' @param mz object to return if \code{condition == "mz"}
#' @param rt object to return if \code{condition == "retentionTime"}
#' @param clust object to return if \code{condition == "clustering"}
#'
#' @details
#' Helper function for \code{shinyCircos}, \code{replayPlotOrder} and 
#' \code{replayPlotAdd}.
#'
#' @return \code{mz}, \code{rt} or \code{clust} depending on \code{condition}
#'
#' @examples 
#' mz <- 1
#' rt <- 2
#' clust <- 3
#' MetCirc:::select(condition = "mz", mz = mz, rt = rt, clust = clust)
#'
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
select <- function(condition, mz, rt, clust) {

    condition <- match.arg(condition, c("mz", "retentionTime", "clustering"))

    if (condition == "mz") res <- mz
    if (condition == "retentionTime") res <- rt
    if (condition == "clustering") res <- clust
    
    res
}

