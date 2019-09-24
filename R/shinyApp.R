#' @import grDevices
#' @import graphics

#' @name shinyCircos
#' @title Interactive visualisation of similarity and navigation of MS/MS features
#' @description Visualise the similarity of MS/MS features in a reactive 
#'  context. See `Details` the vignette for further descriptions on how to use 
#'  `shinyCircos`.
#' @usage shinyCircos(similarityMatrix, spectra, condition, ...)
#' @param similarityMatrix `matrix`, `similarityMatrix` contains 
#' pair-wise similarity coefficients which give information about the similarity 
#' between MS/MS features
#' @param spectra an S4 object of class `Spectra`, the 
#'  `Spectra` object will be used to display information about the selected 
#'  feature and will store information of annotation
#' @param condition `character` vector, specifies which condtions/samples
#' are displayed
#' @param ... further arguments passed to `shinyCircos`, e.g. 
#' `cexFeatureNames` to pass to `plotCircos` to set font size in 
#' `plotCircos` of feature names
#' @details The function is based on the `shiny` and `circlize` 
#' package. 
#' The user can choose interactively thresholds, type of links (between or 
#' within groups), display information about MS/MS features, permanently select 
#' MS/MS features and export selected precursors. The `Spectra` object
#' stores annotation information about the MS/MS features. Names of features 
#' within the `similarityMatrix` have to be found as entries 
#' in `Spectra`. `names(Spectra)` are used as identifiers and 
#' `colnames`/`rownames` from `similarityMatrix` are cleaved 
#' by the group identifier, separated by "_"). Annotation information is taken
#' from `spectra` from the columns "names", "information", "classes" and
#' "adduct" in the slot elementMetadata of `spectra`. After exiting 
#' the application, the annotation will be written to the respective columns
#' in the slot elementMetadata. If one or several of these columns is 
#' already present in elementMetadata, the column(s) will be used as the 
#' source of annotation information.
#' @return `character`, `shinyCircos` returns a `character` vector with the 
#' permanently selected precursors and an object with the `Spectra`
#' object containing the annotation.
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples data("spectra", package = "MetCirc")
#' similarityMat <- compare_Spectra(spectra_tissue[1:10],
#'     fun = normalizeddotproduct, binSize = 0.01)
#' \dontrun{shinyCircos(similarityMatrix = similarityMat, 
#'     spectra = spectra_tissue, condition = c("SPL", "LIM", "ANT", "STY))}
#' @export
shinyCircos <- function(similarityMatrix, spectra, condition, ...) {
    
    ## check if names, information, classes or adduct is present in
    ## @elementMetadata and add respectively, these slots will be querid
    if (is.null(spectra@elementMetadata$names))
        spectra@elementMetadata <- cbind(spectra@elementMetadata,
                                                    names = "Unknown")
    if (is.null(spectra@elementMetadata$information))
        spectra@elementMetadata <- cbind(spectra@elementMetadata,
                                                    information = "Unknown")
    if (is.null(spectra@elementMetadata$classes))
        spectra@elementMetadata <- cbind(spectra@elementMetadata,
                                                    classes = "Unknown")
    if (is.null(spectra@elementMetadata$adduct))
        spectra@elementMetadata <- cbind(spectra@elementMetadata,
                                                    adduct = "Unknown")
    
    ## circlize parameters
    circos.par(gap.degree = 0, cell.padding = c(0, 0, 0, 0),
                                                    track.margin = c(0.0, 0))
    
    ## order spectra per condition according to rt, mz, cluster
    rt <- unlist(lapply(spectra@listData, function(x) x@rt))
    prec_mz <- unlist(lapply(spectra@listData, function(x) x@precursorMz))
    
    groupname <- rownames(similarityMatrix)
    
    ## create plots and assign to objects by recordPlot
    ## rt
    RT <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMatrix,
        spectra = spectra, type = "retentionTime", condition = condition)
    link0dfRT <- RT[["link0df"]]
    rt_match <- RT[["type_match"]]
    
    ## get group and name from rt_match
    ## rt_match is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    groupRT <- lapply(strsplit(rt_match, split = "_"), "[", 1)
    groupRT <- unlist(groupRT)
    nameRT <- lapply(strsplit(rt_match, split = "_"), function (x) x[length(x)])
    nameRT <- unlist(nameRT)
    
    ## plot filled and get degree features
    fillRT <- MetCirc:::recordPlotFill_degreeFeatures(rt_match, ...)
    degFeatRT <- fillRT[["degreeFeatures"]]
    fillRT <- fillRT[["plotFill"]]
    
    ## plot highlight
    highlightRT <- MetCirc:::recordPlotHighlight(rt_match, ...)

    ## mz
    MZ <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMatrix,
        spectra = spectra, type = "mz", condition = condition)
    link0dfMZ <- MZ[["link0df"]]
    mz_match <- MZ[["type_match"]]
    
    ## get group and name from mz_match
    ## mz_match is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    groupMZ <- lapply(strsplit(mz_match, split = "_"), "[", 1)
    groupMZ <- unlist(groupMZ)
    nameMZ <- lapply(strsplit(mz_match, split = "_"), function (x) x[length(x)])
    nameMZ <- unlist(nameMZ)
    
    ## plot filled and get degree features
    fillMZ <- MetCirc:::recordPlotFill_degreeFeatures(mz_match, ...)
    degFeatMZ <- fillMZ[["degreeFeatures"]]
    fillMZ <- fillMZ[["plotFill"]]
    
    ## plot highlight
    highlightMZ <- MetCirc:::recordPlotHighlight(mz_match, ...)
    
    ## clustering
    Clust <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMatrix,
        spectra = spectra, type = "clustering", condition = condition)
    link0dfClust <- Clust[["link0df"]]
    clust_match <- Clust[["type_match"]]
    
    ## get group and name from clust_match
    ## clust_match is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    groupClust <- lapply(strsplit(clust_match, split = "_"), "[", 1)
    groupClust <- unlist(groupClust)
    nameClust <- lapply(strsplit(clust_match, split = "_"), function (x) x[length(x)])
    nameClust <- unlist(nameClust)
    
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
    
    
    ui <- fluidPage(
        tags$head(tags$script('
                $(document).on("shiny:connected", function(e) {
                Shiny.onInputChange("innerWidth", window.innerWidth);
                });
                $(window).resize(function(e) {
                Shiny.onInputChange("innerWidth", window.innerWidth);
                });
                '
        )),
        column(4,
            fluidRow(
                tabsetPanel(id = "tabs",
                    tabPanel("Main", id = "Main", wellPanel(
                        radioButtons("choiceLinks", "choose type of links",
                            choices = c("all" = "all", "inter-class links" = "inter",
                                "intra-class links" = "intra"), selected = "all"),
                        sliderInput("threshold",
                            "Threshold for similarity to display",
                                   min = 0, max = 1, value = c(0.8, 1)),
                        radioButtons("order", "order within groups",
                            choices = c("clustering" = "clustering", "m/z" = "mz",
                              "retention time" = "retentionTime"), selected = "mz"),
                        uiOutput("annotationName"),
                        uiOutput("annotationClass"),
                        uiOutput("annotationInformation"),
                        uiOutput("annotationAdduct"),
                        uiOutput("annotationButton"),
                        actionButton("resetClickIndices", "Unselect features"),
                        actionButton("stop", "Stop and export \n selected features")
                    )),
                    tabPanel("Spectra", id = "Spectra", wellPanel(
                        uiOutput("subjectSpectra"),
                        uiOutput("querySpectra")
                    )),
                    tabPanel("Appearance", id = "Appearance", wellPanel(
                        sliderInput("plotSize", "plot size", 
                           min = 0.5, max = 1.5, value = 1),
                        sliderInput("precision", "precision of numbers",
                           value = 2, min = 0, max = 5, step = 1),
                        checkboxInput("legend", "legend", value = FALSE)
                    ))
                               
                )
            ),
            plotOutput("circosLegend", height = "300")
        ),
        column(8,
            conditionalPanel(condition = "input.tabs  == 'Spectra'",
                plotOutput("plotSpec")), 
            fluidRow(uiOutput("sized_plot")),
            fluidRow(
                verbatimTextOutput("dimension_display"),
                htmlOutput("sglConnectedFeature"),
                verbatimTextOutput("dblClickFeature")
            )
        )
    )
    
    server <- function(input, output, session) {
        
        ## annotation
        ## reactiveValues for spectra, stores annotation information
        spe <- reactiveValues(spectra = spectra)
        
        output$annotationName <- renderUI({
            if (length(indSgl$ind) > 0 && onCircle$is) {
                textInput("names", label = "name", 
                    value = isolate(spe$spectra@elementMetadata$names[ind()]))
            } else NULL  
        })
        
        output$annotationClass <- renderUI({
            if (length(indSgl$ind) > 0 && onCircle$is) {
                textInput("classes", label = "class",
                    value = isolate(spe$spectra@elementMetadata$classes[ind()]))
            } else NULL 
        })
        
        output$annotationInformation <- renderUI({
            if (length(indSgl$ind) > 0 && onCircle$is) {
                #if (onCircle$is) {
                textInput("information", label = "information",
                    value = isolate(spe$spectra@elementMetadata$information[ind()]))
            } else NULL  
        })
        output$annotationAdduct <- renderUI({
            if (length(indSgl$ind) > 0 && onCircle$is) {
                textInput("adduct", label = "adduct",
                    value = isolate(spe$spectra@elementMetadata$adduct[ind()]))
            } else NULL  
        })
        
        
        output$annotationButton <- renderUI({
            if (length(indSgl$ind) > 0 && onCircle$is)
                actionButton(inputId = "annotate", label = "update annotation")
        })
        
        ind <- reactive({
            if (length(indSgl$ind) > 0) {
                nameSgl <- GN()[indSgl$ind]
                nameSgl <- strsplit(nameSgl, split = "_")
                nameSgl <- lapply(nameSgl, "[", -1)
                nameSgl <- lapply(nameSgl, function(x) paste(x, collapse = "_"))
                nameSgl <- unlist(nameSgl)
                which(nameSgl == names(spectra))
            }
        })
        
        ## eventReactive for annotation ind
        indAnn <- eventReactive(input$annotate, {ind()})
        
        ## eventReactive for input$name (records "character" for annotation)
        annotateNames <- eventReactive(input$annotate, {
            input$names ##as.character(input$names)
        })
        ## eventReactive for input$classes
        annotateClasses <- eventReactive(input$annotate, {input$classes})
        
        ## eventReactive for input$information
        annotateInformation <- eventReactive(input$annotate, {input$information})
        
        ## eventReactive for input$adduct
        annotateAdduct <- eventReactive(input$annotate, {input$adduct})
        
        ## observe annotations and write to respective columns in slot
        ## elementMetadata
        observe({
            spe$spectra@elementMetadata$names[indAnn()] <- annotateNames()
        })
        observe({
            spe$spectra@elementMetadata$classes[indAnn()] <- annotateClasses()
        })
        observe({
            spe$spectra@elementMetadata$information[indAnn()] <- annotateInformation()
        })
        observe({
            spe$spectra@elementMetadata$adduct[indAnn()] <- annotateAdduct()
        })
        ## end annotation
        
        ## ordering of features, use predefined groupname object
        GN <- reactive({
            MetCirc:::select(input$order, mz_match, rt_match, clust_match)
        })
        
        ## get degree of features
        degFeat <- reactive({
            MetCirc:::select(input$order, degFeatMZ, degFeatRT, degFeatClust)
        })
        
        ## get link0df
        link0df <- reactive({
            MetCirc:::select(input$order, link0dfMZ, link0dfRT, link0dfClust)
        })
        
        ## create reactive expression for linkDf which is cut according to
        ## set radioButton (input$choiceLinks)
        linkDf_cut <- reactive(
            cutLinkDf(link0df(), type = input$choiceLinks))
        
        ## threshold linkDf_cut
        linkDf_threshold <- reactive(
            thresholdLinkDf(linkDf_cut(),
                            input$threshold[1], input$threshold[2]))
        
        ## reactiveValues for click coordinates
        coordsNewSgl <- reactiveValues(X = 0, Y = 0)
        coordsOldSgl <- reactiveValues(X = 0, Y = 0)
        
        observe({
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
        
        ## is mouse over the track 1?
        onCircle <- reactiveValues(is = NULL)
        observe({
            if (!is.null(coordsNewSgl$X)) {
                .dist <- sqrt(coordsOldSgl$X^2 + coordsOldSgl$Y^2)
                if (.dist >= 0.8 & .dist <= 1) {
                    onCircle$is <- TRUE
                } else {
                    onCircle$is <- FALSE
                }
            } else onCircle$is <- FALSE
        })
        
        ## single click: which is the current sector?
        indSgl <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosSgl$x))
                indSgl$ind <- minFragCart2Polar(input$circosSgl$x,
                                        input$circosSgl$y, degFeat())
        })
        
        ## double click: which is the current sector?
        coordsNewDbl <- reactiveValues(X = 0, Y = 0)
        coordsOldDbl <- reactiveValues(X = 0, Y = 0)
        
        observe({
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
        indDbl <- reactiveValues(ind = NULL, new = NULL)
        observe({
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
        
        observe({
            input$resetClickIndices
            isolate(indDblMZ$ind <- NULL)
            isolate(indDblRT$ind <- NULL)
            isolate(indDblCluster$ind <- NULL)
            isolate(indDbl$new <- NULL)
            isolate(onCircle$is <- FALSE)
            isolate(indSgl$ind <- NULL)
        })
        
        ## reset indSgl when changing radio button order
        observe({
            input$order
            isolate(onCircle$is <- FALSE)
            isolate(indSgl$ind <- NULL)
        })
        
        
        ## write double-clicked (truncated) names to indDblMZ, indDblRT,
        ## indDblCluster
        indDblMZ <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosDbl$x)) {
                if (!is.null(indDbl$new)) {
                    newMZ <- paste(groupMZ, nameMZ, sep = "_")
                    newIndMZ <- match(indDbl$new, newMZ)
                    if (isolate(newIndMZ %in% indDblMZ$ind)) {
                        indDblMZ$ind <- isolate(indDblMZ$ind[-which(newIndMZ == indDblMZ$ind)])
                    } else {indDblMZ$ind <- isolate(c(indDblMZ$ind, newIndMZ))}
                }
            }
        })
        
        indDblRT <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosDbl$x)) {
                if (!is.null(indDbl$new)) {
                    newRT <- paste(groupRT, nameRT, sep = "_")
                    newIndRT <- match(indDbl$new, newRT)
                    if (isolate(newIndRT %in% indDblRT$ind)) {
                        indDblRT$ind <- isolate(indDblRT$ind[-which(newIndRT == indDblRT$ind)])
                    } else {indDblRT$ind <- isolate(c(indDblRT$ind, newIndRT))}
                }
            }
        })
        
        indDblCluster <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosDbl$x)) {
                if(!is.null(indDbl$new)) {
                    newCl <- paste(groupClust, nameClust, sep = "_")
                    newIndCl <- match(indDbl$new, newCl)
                    
                    if (isolate(newIndCl %in% indDblCluster$ind)) {
                        indDblCluster$ind <- isolate(indDblCluster$ind[-which(newIndCl == indDblCluster$ind)])
                    } else {indDblCluster$ind <- isolate(c(indDblCluster$ind, newIndCl))}
                }
            }
        })
        
        ## plotting
        initializePlot <- reactive({
            circos.initialize(factor(GN(), levels = GN()),
                xlim = matrix(rep(c(0,1), length(mz_match)), ncol = 2, byrow = TRUE))
            circos.trackPlotRegion(factor(GN(), levels = GN()), ylim = c(0,1))  
        })
        
        ## assign to output$circos: actual plotting
        output$circos <- renderPlot({
            indDblMZ$ind
            initializePlot()
            MetCirc:::replayPlotOrder(orderMatch = input$order,
                onCircle = onCircle$is, plot_l = plot_l, ind = indDblMZ$ind)
            MetCirc:::replayPlotAdd(orderMatch = input$order, onCircle = onCircle$is,
                linkDf = linkDf_threshold(), mz_match = mz_match, rt_match = rt_match,
                clust_match = clust_match, ind = indSgl$ind, indMz = indDblMZ$ind,
                indRT = indDblRT$ind, indCluster = indDblCluster$ind)
        })

        
        output$sized_plot <- renderUI({
            plotOutput("circos",
                dblclick = "circosDbl", click = "circosSgl",
                width = ifelse(is.null(input$innerWidth), 0, input$innerWidth * 0.5 * input$plotSize),
                height = ifelse(is.null(input$innerWidth), 0, input$innerWidth * 0.5 * input$plotSize))
        })
        
        output$circosLegend <- renderPlot({
            if (!is.null(input$legend)) if(input$legend)
                circosLegend(rt_match, highlight = TRUE)
        })
        
        ## show when Clicking the feature which connects to it
        linkDfIndsSgl <- reactive({
            getLinkDfIndices(GN()[indSgl$ind], linkDf_threshold())
        })
        
        output$sglConnectedFeature <- renderUI({
            if (!is.null(onCircle$is)) {
                if (onCircle$is)
                    HTML(printInformationSelect(select = GN()[indSgl$ind],
                            spectra = spe$spectra, linkDfInd = linkDfIndsSgl(),
                            linkDf = linkDf_threshold(),
                            similarityMatrix = similarityMatrix,
                            roundDigits = input$precision))
            }
        })
        
        output$dblClickFeature <- renderText({
            if (length(indDblMZ$ind) > 0) 
                c("(permanently) selected features: ", mz_match[indDblMZ$ind])
            else "no features (permanently) selected"
        })
        
        ## Tab Spectra
        output$querySpectra <- renderUI({
            selectInput("query", "Choose MS/MS spectra #2:", choices = GN())
        })
        
        output$subjectSpectra <- renderUI({
            selectInput("subject", "Choose MS/MS spectra #1:", choices = GN(), selected = GN()[indSgl$ind])
        })
        
        output$plotSpec <- renderPlot({
            if (!is.null(input$subject) && !is.null(input$query)) {
                plotSpectra(spectra, input$subject, input$query)
            }
        })
        
        ## on exit
        observe({
            if (input$stop == 0) {return()
            } else {
                circos.clear()
                selFeat <- as.character(paste(
                    groupMZ[indDblMZ$ind], nameMZ[indDblMZ$ind],
                    sep = "_"))
                stopApp(list(spectra = spe$spectra, selectedFeatures = selFeat))
            }
        })
        
        
    }
    
    app <- list(ui = ui, server = server)
    runApp(app)
}

#' @name printInformationSelect
#' @title Display information on connected features of selected features
#' @description Displays information on connected features of selected features.
#' @usage printInformationSelect(select, spectra = NULL,
#'     linkDfInd, linkDf, similarityMatrix, roundDigits = 2) 
#' @param select `character`, obtained from groupname, `character` of 
#'     selected feature
#' @param spectra `Spectra` object containing spectra that are compared
#' in `similarityMatrix`
#' @param linkDfInd `numeric` indices of selected features
#' @param linkDf `data.frame` that contains information of linked 
#'  features for given thresholds
#' @param similarityMatrix `matrix` that is used to get information on the 
#' degree of similarity, `similarityMat` is an ordered version of a 
#' similarity matrix, see `?createOrderedSimMat`
#' @param roundDigits `numeric`,  how many digits should be displayed?
#' @details `printInformationSelect` is for internal use. 
#' @return `character` that is in HTML format
#' @examples
#' data("spectra", package = "MetCirc")
#' similarityMat <- compare_Spectra(spectra_tissue[1:10],
#'     fun = normalizeddotproduct, binSize = 0.01)
#' linkDf <- createLinkDf(similarityMatrix = similarityMat,
#'     spectra = spectra_tissue[1:10], 
#'     condition = c("SPL", "LIM", "ANT", "STY"), lower = 0.5, upper = 1)
#' ## cut link data.frame (here: only display links between groups)
#' linkDf_cut <- cutLinkDf(linkDf, type = "inter")
#' groupname <- c(as.character(linkDf_cut[, "spectrum1"]),
#'             as.character(linkDf_cut[, "spectrum2"]))
#' groupname <- unique(groupname)
#' ## arbitrarily select a feature
#' ind <- 2
#' linkDfInds <- getLinkDfIndices(groupname[ind], linkDf_cut)
#' MetCirc:::printInformationSelect(groupname[ind], spectra = spectra_tissue[1:10],
#'     linkDfInd = linkDfInds, linkDf = linkDf_cut, similarityMatrix = similarityMat)
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
printInformationSelect <- function(select, spectra = NULL, 
                linkDfInd, linkDf, similarityMatrix, roundDigits = 2) {

    ## connected features: find
    connect <- as.character(unique(unlist(linkDf[linkDfInd, c("spectrum1", "spectrum2")])))
    
    ## remove duplicated hovFeat in connect
    if (select %in% connect) connect <- connect[-which(connect == select)]
    
    ## truncate group
    connect_cut <- lapply(strsplit(connect, split = "_"), "[", -1)
    connect_cut <- unlist(lapply(connect_cut, function(x) paste(x, collapse = "_")))
    select_cut <- lapply(strsplit(select, split = "_"), "[", -1)
    select_cut <- unlist(lapply(select_cut, function(x) paste(x, collapse = "_")))
    
    ## truncate spectra for select and connecting features
    select_spectra <- spectra[select_cut, ]
    connect_spectra <- spectra[connect_cut, ]
    
    if (length(linkDfInd) == 0) {
        return(paste0(select, " (",
            round(select_spectra[[1]]@precursorMz, roundDigits), ", ",
            round(select_spectra[[1]]@rt, roundDigits), ", ",
            select_spectra@elementMetadata$names, ", ",
            select_spectra@elementMetadata$information, ", ",
            select_spectra@elementMetadata$classes, ", ",
            select_spectra@elementMetadata$adduct,
            ") does not connect to any feature "))
    } else { ## if length > 0
        
        connChar <- character()
        degreeSimilarity <- similarityMatrix[select_cut, ]
            
        for (i in 1:length(connect_cut)) {
        
            connect_i <- connect_cut[i]
            degreeSimilarityI <- round(degreeSimilarity[connect_i], 3)
            degreeSimilarityI <- as.numeric(degreeSimilarityI)
            
            spectra_i <- connect_spectra[[i]]
            
            newFeat <- paste0(connect[i], " (", degreeSimilarityI, ", ",
                round(spectra_i@precursorMz, roundDigits), ", ",
                round(spectra_i@rt, roundDigits), ", ",
                connect_spectra@elementMetadata$names[i], ", ",
                connect_spectra@elementMetadata$information[i], ", ",
                connect_spectra@elementMetadata$classes[i], ", ",
                connect_spectra@elementMetadata$adduct[i], ")", "<br/>")
                
            connChar <- c(connChar, newFeat)
        }
            
        connChar <- paste(connChar, collapse = " ")

        return(paste0(select, " (", 
            round(select_spectra[[1]]@precursorMz, roundDigits), ", ",
            round(select_spectra[[1]]@rt, roundDigits), ", ",
            select_spectra@elementMetadata$names, ", ",
            select_spectra@elementMetadata$information, ", ",
            select_spectra@elementMetadata$classes, ", ",
            select_spectra@elementMetadata$adduct, ") connects to ",
            " <br/>", connChar))
    }
}


#' @name replayPlotOrder
#' @title Wrapper for `replayPlot`
#' @description `replayPlotOrder` will call `replayPlot` from `grDevices` with 
#' a `recordedplot` object based on `orderMatch`.
#' @usage replayPlotOrder(orderMatch = "mz", onCircle = FALSE, plot_l, ind)
#' @param orderMatch `character`, either `"mz"`, `"retentionTime"` or 
#' `"clustering"`
#' @param plot_l `list` with plots
#' @param onCircle `logical`, are coordinates on circle. If FALSE and
#' no features are selected (`length(ind) == 0`), then filled plots are 
#' replayed, otherwise highlighted plots are replayed.
#' @param ind `numeric`, indices of clicked features
#' @details Helper function for `shinyCircos`.
#' @return `replayedplot`
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
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
replayPlotOrder <- function(orderMatch = "mz", onCircle = FALSE, plot_l, ind) {
    
    if (!onCircle & length(ind) == 0) {
        
        ## get plot based on orderMatch 
        p <- MetCirc:::select(orderMatch, plot_l[["fillMz"]],
                    plot_l[["fillRT"]], plot_l[["fillClust"]])
        ## plot
        replayPlot(p)
    } else {
         
        ## get plot based on orderMatch
        p <- MetCirc:::select(orderMatch, plot_l[["highlightMz"]],
                    plot_l[["highlightRT"]], plot_l[["highlightClust"]])
        
        ## plot
        replayPlot(p)
    }
}

#' @name replayPlotAdd
#' @title Plot plotCircos or highlight
#' @description `replayPlotAdd` plots additional plots on a plot, either
#' plots `plotCircos` or `highlight`.
#' @param orderMatch orderMatch `character`, either `"mz"`, `"retentionTime"` 
#' or `"clustering"`
#' @param onCircle `logical`, are coordinates on circle. If FALSE and
#' no features are selected (`length(ind) == 0`), then filled plots are 
#' replayed, otherwise highlighted plots are replayed.
#' @param linkDf `data.frame` that contains information of linked 
#'  features for given thresholds
#' @param mz_match `character`, ordered vector according to m/z
#' @param rt_match `character`, ordered vector according to retention time
#' @param clust_match `character`, ordered vector according to clustering
#' @param ind  `numeric`, indices of clicked features
#' @param indMz `numeric`, indices of clicked features for `"mz"` ordering
#' @param indRT `numeric`, indices of clicked features for `"retentionTime"` 
#' ordering
#' @param indCluster `numeric`, indices of clicked features for `"clustering"` 
#' ordering
#' @details Helper function for `shinyCircos`.
#' @return Depending on `onCircle` and `indMz` either returns `plotCircos` or 
#' `highlight`
#' @examples
#' data("spectra", package = "MetCirc")
#' similarityMat <- compare_Spectra(spectra_tissue[1:10], 
#'     fun = normalizeddotproduct, binSize = 0.01)
#' ## order according to retention time 
#' mz_match <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMat, 
#'     spectra = spectra_tissue, type = "mz", 
#'     condition = c("SPL", "LIM", "ANT", "STY"))
#' linkDf <- mz_match[["link0df"]]
#' mz_match <- mz_match[["type_match"]]
#' rt_match <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMat, 
#'     spectra = spectra_tissue, type = "retentionTime", 
#'     condition = c("SPL", "LIM", "ANT", "STY"))
#' rt_match <- rt_match[["type_match"]]
#' clust_match <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMat, 
#'     spectra = spectra_tissue, type = "clustering", 
#'     condition = c("SPL", "LIM", "ANT", "STY"))
#' clust_match <- clust_match[["type_match"]]
#' circos.initialize(mz_match,##, levels  =  mz_match),
#'     xlim = matrix(rep(c(0,1), length(mz_match)), ncol = 2, byrow = TRUE))
#' #circos.trackPlotRegion(factor(mz_match, levels = mz_match), ylim = c(0,1))  
#' MetCirc:::replayPlotAdd(orderMatch = "mz", onCircle = FALSE, linkDf = linkDf, 
#'     mz_match = mz_match, rt_match = rt_match, clust_match = clust_match, 
#'     ind = 1, indMz = NULL, indRT = NULL, indCluster = NULL)
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
replayPlotAdd <- function(orderMatch = "mz", onCircle = FALSE, linkDf,
        mz_match, rt_match, clust_match, ind, indMz, indRT, indCluster) {
    
    ## get type_match based on orderMatch
    type_match <- MetCirc:::select(orderMatch, mz_match, rt_match, clust_match)
    inds <- MetCirc:::select(orderMatch, indMz, indRT, indCluster)
    
    if (!onCircle & length(inds) == 0) {
        ## plot
        return(plotCircos(type_match, linkDf, initialize = FALSE,
            featureNames = FALSE, groupSector = FALSE, groupName = FALSE,
            links = TRUE, highlight = FALSE))

    } else {
        ## get inds
        inds <- if (onCircle) c(ind, inds) else inds
        
        ## plot
        return(highlight(type_match, inds, linkDf))
    }
}


#' @name recordPlotFill_degreeFeatures
#' @title Record a plot of filled features and the degree of features
#' @description  `recordPlotFill_degreeFeatures` records a plot of filled 
#' features and returns the degree of features
#' @usage recordPlotFill_degreeFeatures(type_match, ...)
#' @param type_match `character`, ordered vector according to type
#' @param ... further arguments passed to `plotCircos`
#' @details Helper function for `shinyCircos`.
#' @return `list` of length 2, entry `plotFill` is of `recordedplot` and 
#' entry `degreeFeatures` that is a `list` of vectors of `numeric(1)` 
#' @examples 
#' type_match <- c("a_1", "a_2", "a_3", "b_1", "b_2", "b_3", "c_1", "c_2")
#' MetCirc:::recordPlotFill_degreeFeatures(type_match)
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
recordPlotFill_degreeFeatures <- function(type_match, ...) {
    plotCircos(type_match, NULL, initialize = TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = FALSE, ...)
    fill <- recordPlot()
    ## get degree of features
    degree <- lapply(type_match, function(x) {
        mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")])
    })
    plot.new()
    return(list(plotFill = fill, degreeFeatures = degree))
    
}

#' @name recordPlotHighlight
#' @title Return a `recordedplot` of `plotCircos` plot with `highlight = TRUE`
#' @description `recordPlotHighlight` returns a `recordedplot` object of 
#' `plotCircos` with `highlight = TRUE`
#' @usage recordPlotHighlight(type_match, ...)
#' @param type_match `character`, ordered vector according to type
#' @param ... further arguments passed to `plotCircos`
#' @details Helper function for `shinyCircos`.
#' @return `recordedplot`
#' @examples 
#' type_match <- c("a_1", "a_2", "a_3", "b_1", "b_2", "b_3", "c_1", "c_2")
#' MetCirc:::recordPlotHighlight(type_match)
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
recordPlotHighlight <- function(type_match, ...) {
    ## use plotCircos
    plotCircos(type_match, NULL, initialize = TRUE, featureNames = TRUE,
               groupSector = TRUE, groupName = FALSE, links = FALSE,
               highlight = TRUE, ...)
    highlight <- recordPlot()
    plot.new()
    return(highlight)
}

#' @name typeMatch_link0
#' @title Get typeMatch and link0 data frame
#' @description `typeMatch_link0` returns a list with accessors `"link0df"` and
#' `"type_match"`
#' @usage typeMatch_link0(similarityMatrix, spectra, type, condition)
#' @param similarityMatrix `matrix` with pair-wise similarity values
#' @param spectra `Spectra` object
#' @param type `character`, either `"mz"`, `"retentionTime"`, `"clustering"`
#' @param condition `character`, tissue
#' @details Helper function for `shinyCircos`. 
#' @return `list` of length 2,  entry `link0df` is of `data.frame` and 
#' entry `type_match` that is a `character` vector
#' @examples 
#' data("spectra", package = "MetCirc")
#' similarityMat <- compare_Spectra(spectra_tissue[1:10], 
#'     fun = normalizeddotproduct, binSize = 0.01)
#' ## order according to retention time 
#' MetCirc:::typeMatch_link0(similarityMatrix = similarityMat, 
#'     spectra = spectra_tissue, type = "mz", 
#'     condition = c("SPL", "LIM", "ANT", "STY"))
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
typeMatch_link0 <- function(similarityMatrix, spectra, type, condition) {
    simMat <- orderSimilarityMatrix(similarityMatrix, spectra = spectra, 
                                    type = type, group = FALSE)
    ## get names of spectra per condition
    inds <- MetCirc:::spectraCond(spectra, condition)
    
    link0df <- createLink0df(simMat, spectra = spectra, condition = condition)
    groupname <- rownames(simMat)
    type_match <- lapply(inds, function(x) {type_match <- match(groupname, x)
    type_match <- type_match[!is.na(type_match)]; x[type_match]})
    type_match <- lapply(seq_along(type_match), function(x) {
        if (length(type_match[[x]]) > 0) {
            paste(condition[x], type_match[[x]], sep = "_")
        } else character()
    })
    type_match <- unique(unlist(type_match))
    return(list(link0df = link0df, type_match = type_match))
}

#' @name select
#' @title Select variable based on condition
#' @description `select` returns `mz`, `rt` or `clust` depending on `condition`.
#' @usage select(condition, mz, rt, clust)
#' @param condition `character`, either `"mz"`, `"retentionTime"` or 
#' `"clustering"`
#' @param mz object to return if `condition == "mz"`
#' @param rt object to return if `condition == "retentionTime"`
#' @param clust object to return if `condition == "clustering"`
#' @details Helper function for `shinyCircos`, `replayPlotOrder` and 
#' `replayPlotAdd`-
#' @return `mz`, `rt` or `clust` depending on condition
#' @examples 
#' mz <- 1
#' rt <- 2
#' clust <- 3
#' MetCirc:::select(condition = "mz", mz = mz, rt = rt, clust = clust)
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
select <- function(condition, mz, rt, clust) {
    if (!condition %in% c("mz", "retentionTime", "clustering"))
        stop("condition not equal to 'mz', 'retentionTime' or 'clustering'")
    if (!is.null(condition)) {
        if (condition == "mz") res <- mz
        if (condition == "retentionTime") res <- rt
        if (condition == "clustering") res <- clust
    } else {
        res <- NULL
    }
    return(res)
}

