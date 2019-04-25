#' @import grDevices
#' @import graphics

#' @name shinyCircos
#' @title Interactive visualisation of similarity and navigation of MS/MS features
#' @description Visualise the similarity of MS/MS features in a reactive 
#'  context. See \code{Details} the vignette for further descriptions on how to use 
#'  \code{shinyCircos}.
#' @usage shinyCircos(similarityMatrix, msp = NULL, ...)
#' @param similarityMatrix \code{matrix}, \code{similarityMatrix} contains 
#' pair-wise similarity coefficients which give information about the similarity 
#' between MS/MS features
#' @param spectra an S4 object of class \code{Spectra}, the 
#'  \code{Spectra} object will be used to display information about the selected 
#'  feature and will store information of annotation
#' @param ... further arguments passed to \code{shinyCircos}, e.g. 
#' \code{cexFeatureNames} to pass to \code{plotCircos} to set font size in 
#' \code{plotCircos} of feature names
#' @details The function is based on the \code{shiny} and \code{circlize} 
#' package. 
#' The user can choose interactively thresholds, type of links (between or 
#' within groups), display information about MS/MS features, permanently select 
#' MS/MS features and export selected precursors. The \code{Spectra} object
#' stores annotation information about the MS/MS features. Names of features 
#' within the \code{similarityMatrix} have to be found as entries 
#' in \code{Spectra}. \code{names(Spectra)} are used as identifiers and 
#' \code{colnames}/\code{rownames} from \code{similarityMatrix} are cleaved 
#' by the group identifier, separated by "_"). Annotation information is taken
#' from \code{spectra} from the columns "names", "information", "classes" and
#' "adduct" in the slot elementMetadata of \code{spectra}. After exiting 
#' the application, the annotation will be written to the respective columns
#' in the slot elementMetadata. If one or several of these columns is 
#' already present in elementMetadata, the column(s) will be used as the 
#' source of annotation information. 
#' @return \code{shinyCircos} returns a \code{character} vector with the 
#' permanently selected precursors and an object with the \code{Spectra}
#' object containing the annotation. 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("idMSMStoMSP", package = "MetCirc")
#' ## truncate files
#' finalMSP <- finalMSP[c(1:20, 29:48, 113:132, 240:259)]
#' data("binnedMSP", package = "MetCirc")
#' binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' \dontrun{shinyCircos(similarityMatrix = similarityMat, msp = finalMSP)}
#' @export
shinyCircos <- function(similarityMatrix, spectra, condition, ...) {
    
    #cond_inds <- lapply(condition, function(x) which(spectra@elementMetadata@listData[[x]] == 1))
    inds <- MetCirc:::spectraCond(spectra, condition)
    
    ## check if names, information, classes or adduct is present in 
    ## @elementMetadata and add respectively, these slots will be querid 
    if (is.null(spectra@elementMetadata$names)) 
        spectra@elementMetadata <- cbind(spectra@elementMetadata, names="Unknown")
    if (is.null(spectra@elementMetadata$information)) 
        spectra@elementMetadata <- cbind(spectra@elementMetadata, information="Unknown")
    if (is.null(spectra@elementMetadata$classes)) 
        spectra@elementMetadata <- cbind(spectra@elementMetadata, classes="Unknown")
    if (is.null(spectra@elementMetadata$adduct)) 
        spectra@elementMetadata <- cbind(spectra@elementMetadata, adduct="Unknown")
    
    ## circlize parameters
    circos.par(gap.degree=0, cell.padding=c(0, 0, 0, 0), track.margin=c(0.0, 0))
    
    ## order spectra per condition according to rt, mz, cluster
    rt <- unlist(lapply(spectra@listData, function(x) x@rt))
    prec_mz <- unlist(lapply(spectra@listData, function(x) x@precursorMz))
    
    groupname <- rownames(similarityMatrix)
    ## create plots and assign to objects by recordPlot
    ## rt
    simMatRT <- orderSimilarityMatrix(similarityMatrix, spectra=spectra, 
        type="retentionTime", group=FALSE)
    link0dfRT <- createLink0df(simMatRT, spectra=spectra, condition=condition)
    groupnameRT <- rownames(simMatRT)
    rt_match <- lapply(inds, function(x) {rt_match <- match(groupnameRT, x)
    rt_match <- rt_match[!is.na(rt_match)]; x[rt_match]})
    rt_match <- lapply(seq_along(rt_match), function(x) {
        if (length(rt_match[[x]]) > 0) {
            paste(condition[x], rt_match[[x]], sep="_")    
        } else character()
    })
    rt_match <- unique(unlist(rt_match))
    
    plotCircos(rt_match, NULL, initialize=TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = FALSE, ...)
    PlotFilledRT <- recordPlot()
    ## get group and name from rt_match
    ## rt_match is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    groupRT <- lapply(strsplit(rt_match, split = "_"), "[", 1)
    groupRT <- unlist(groupRT)
    nameRT <- lapply(strsplit(rt_match, split = "_"), function (x) x[length(x)])
    nameRT <- unlist(nameRT)
    ## get degree of features
    degreeFeaturesRT <- lapply(rt_match, 
        function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
    plot.new()
    
    plotCircos(rt_match, NULL, initialize=TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = TRUE, ...)
    PlotHighlightRT <- recordPlot()
    plot.new()
    
    ## mz
    simMatMZ <- orderSimilarityMatrix(similarityMatrix, spectra=spectra, 
        type="mz", group=FALSE)
    link0dfMZ <- createLink0df(simMatMZ, spectra=spectra, condition=condition)
    groupnameMZ <- rownames(simMatMZ)
    mz_match <- lapply(inds, function(x) {mz_match <- match(groupnameMZ, x)
    mz_match <- mz_match[!is.na(mz_match)]; x[mz_match]})
    mz_match <- lapply(seq_along(mz_match), function(x) {
        if (length(mz_match[[x]]) > 0) {
            paste(condition[x], mz_match[[x]], sep="_")    
        } else character()
    })
    mz_match <- unique(unlist(mz_match))
    plotCircos(mz_match, NULL, initialize=TRUE, featureNames=TRUE, 
               groupSector=TRUE, groupName=FALSE, links=FALSE, highlight=FALSE, ...)
    PlotFilledMZ <- recordPlot()
    ## get group and name from mz_match
    ## mz_match is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    groupMZ <- lapply(strsplit(mz_match, split = "_"), "[", 1)
    groupMZ <- unlist(groupMZ)
    nameMZ <- lapply(strsplit(mz_match, split = "_"), function (x) x[length(x)])
    nameMZ <- unlist(nameMZ)
    ## get degree of features
    degreeFeaturesMZ <- lapply(mz_match, 
        function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
    
    plot.new()
    plotCircos(mz_match, NULL, initialize=TRUE, featureNames = TRUE, 
               groupSector = TRUE, groupName = FALSE, links = FALSE, 
               highlight = TRUE, ...)
    PlotHighlightMZ <- recordPlot()
    plot.new()
    
    ## clustering
    simMatClust <- orderSimilarityMatrix(similarityMatrix, spectra=spectra, 
        type="clust", group=FALSE)
    link0dfClust <- createLink0df(simMatClust, spectra=spectra, condition=condition)
    groupnameClust <- rownames(simMatClust)
    clust_match <- lapply(inds, function(x) {clust_match <- match(groupnameClust, x)
    clust_match <- clust_match[!is.na(clust_match)]; x[clust_match]})
    clust_match <- lapply(seq_along(clust_match), function(x) {
        if (length(clust_match[[x]]) > 0) {
            paste(condition[x], clust_match[[x]], sep="_")    
        } else character()
    })
    clust_match <- unique(unlist(clust_match))
    plotCircos(clust_match, NULL, initialize=TRUE, featureNames=TRUE, 
               groupSector=TRUE, groupName=FALSE, links=FALSE, highlight=FALSE, ...)
    PlotFilledClust <- recordPlot()
    ## get group and name from clust_match
    ## clust_match is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    groupClust <- lapply(strsplit(clust_match, split = "_"), "[", 1)
    groupClust <- unlist(groupClust)
    nameClust <- lapply(strsplit(clust_match, split = "_"), function (x) x[length(x)])
    nameClust <- unlist(nameClust)
    ## get degree of features
    degreeFeaturesClust <- lapply(clust_match,
                                  function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
    plot.new()
    
    plotCircos(clust_match, NULL, initialize=TRUE, featureNames=TRUE, 
               groupSector=TRUE, groupName=FALSE, links=FALSE, highlight=TRUE, ...)
    PlotHighlightClust <- recordPlot()
    plot.new()
    
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
                tabsetPanel(id="tabs",
                    tabPanel("Main", id="Main", wellPanel(
                        radioButtons("choiceLinks", "choose type of links",
                            choices = c("all"="all", "inter-class links"="inter",
                                "intra-class links"="intra"), selected="all"),
                        sliderInput("threshold",
                            "Threshold for similarity to display",
                                   min=0, max=1, value=c(0.8, 1)),
                        radioButtons("order", "order within groups",
                            choices=c("clustering"="clustering", "m/z"="mz", 
                              "retention time"="retentionTime"), selected="mz"),
                        uiOutput("annotationName"),
                        uiOutput("annotationClass"),
                        uiOutput("annotationInformation"),
                        uiOutput("annotationAdduct"),
                        uiOutput("annotationButton"),
                        actionButton("resetClickIndices", "Unselect features"),
                        actionButton("stop", "Stop and export \n selected features")
                    )),
                    tabPanel("Spectra", id="Spectra", wellPanel(
                        uiOutput("subjectSpectra"),
                        uiOutput("querySpectra")
                    )),
                    tabPanel("Appearance", id="Appearance", wellPanel(
                        sliderInput("plotSize", "plot size", 
                           min=0.5, max=1.5, value=1),
                        sliderInput("precision", "precision of numbers", 
                           value=2, min=0, max=5, step=1),
                        checkboxInput("legend", "legend", value=FALSE)
                    ))
                               
                )
            ),
            plotOutput("circosLegend", height="300")
        ),
        column(8,
            conditionalPanel(condition="input.tabs == 'Spectra'", 
                plotOutput("plotSpec")), 
            fluidRow(uiOutput("sized_plot")),
            fluidRow(
                verbatimTextOutput("dimension_display"),
                htmlOutput("clickConnectedFeature"),
                verbatimTextOutput("dblClickFeature")
            )
        )
    )
    
    server <- function(input, output, session) {
        
        ## annotation
        ## reactiveValues for spectra, stores annotation information
        spe <- reactiveValues(spectra=spectra) 
        
        output$annotationName <- renderUI({
            if (length(indClick$ind) > 0 && onCircle$is) {
                textInput("names", label="name", 
                          value=isolate(spe$spectra@elementMetadata$names[ind()]))
            } else NULL  
        })
        
        output$annotationClass <- renderUI({
            if (length(indClick$ind) > 0 && onCircle$is) {
                textInput("classes", label="class", 
                          value=isolate(spe$spectra@elementMetadata$classes[ind()]))
            } else NULL 
        })
        
        output$annotationInformation <- renderUI({
            if (length(indClick$ind) > 0 && onCircle$is) {
                #if (onCircle$is) {
                textInput("information", label="information", 
                          value=isolate(spe$spectra@elementMetadata$information[ind()])) 
            } else NULL  
        })
        output$annotationAdduct <- renderUI({
            if (length(indClick$ind) > 0 && onCircle$is) {
                textInput("adduct", label="adduct", 
                          value=isolate(spe$spectra@elementMetadata$adduct[ind()])) 
            } else NULL  
        })
        
        
        output$annotationButton <- renderUI({
            if (length(indClick$ind) > 0 && onCircle$is) 
                actionButton(inputId="annotate", label="update annotation")
        })
        
        ind <- reactive({
            if (length(indClick$ind) > 0) {
                nameClick <- GN()[indClick$ind]
                nameClick <- strsplit(nameClick, split="_")
                nameClick <- lapply(nameClick, "[", -1)
                nameClick <- lapply(nameClick, function(x) paste(x, collapse="_"))
                nameClick <- unlist(nameClick)
                which(nameClick == names(spectra))
            }
        })
        
        ## eventReactive for annotation ind 
        indAnn <- eventReactive(input$annotate, {
            ind()
        })
        
        ## eventReactive for input$name (records "character" for annotation)
        annotateNames <- eventReactive(input$annotate, {
            as.character(input$names) 
        })
        ## eventReactive for input$classes
        annotateClasses <- eventReactive(input$annotate, {
            as.character(input$classes)
        })
        ## eventReactive for input$information
        annotateInformation <- eventReactive(input$annotate, {
            as.character(input$information)
        })
        ## eventReactive for input$adduct
        annotateAdduct <- eventReactive(input$annotate, {
            as.character(input$adduct)
        })
        
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
            if (input$order == "mz") GN <- mz_match
            if (input$order == "retentionTime") GN <- rt_match
            if (input$order == "clustering") GN <- clust_match
            GN
            
        })
        
        ## get degree of features
        degreeFeatures <- reactive({
            if (input$order == "mz") degFeatures <- degreeFeaturesMZ
            if (input$order == "retentionTime") degFeatures <- degreeFeaturesRT
            if (input$order == "clustering") degFeatures <- degreeFeaturesClust
            degFeatures
        })
        
        ## calculateLink0Matrix
        link0df <- reactive({
            if (!is.null(input$order)) {
                if (input$order == "mz") link0df <- link0dfMZ
                if (input$order == "retentionTime") link0df <- link0dfRT
                if (input$order == "clustering") link0df <- link0dfClust
                link0df
                
            }
        })
        
        ## create reactive expression for LinkDF which is cut according to 
        ## set radioButton (input$choiceLinks)
        LinkDf_cut <- reactive(
            cutLinkDf(link0df(), type=input$choiceLinks))
        
        ## threshold linkDf_cut
        LinkDf_threshold <- reactive(
            thresholdLinkDf(LinkDf_cut(), 
                            input$threshold[1], input$threshold[2]))
        
        ## reactiveValues for click Coordinates
        CoordinatesNewClick <- reactiveValues(X = 0, Y = 0)
        CoordinatesOldClick <- reactiveValues(X = 0, Y = 0)
        
        observe({
            if (!is.null(input$circosClick$x)) {
                CoordinatesNewClick$X <- input$circosClick$x
                CoordinatesNewClick$Y <- input$circosClick$y
                CoordinatesOldClick$X <- CoordinatesNewClick$X
                CoordinatesOldClick$Y <- CoordinatesNewClick$Y
            } else {
                CoordinatesNewClick$X <- CoordinatesOldClick$X
                CoordinatesNewClick$Y <- CoordinatesOldClick$Y
            }
        })
        
        ## is mouse over the track 1?
        onCircle <- reactiveValues(is = NULL)
        observe({
            if (!is.null(CoordinatesNewClick$X)) {
                .dist <- sqrt(CoordinatesOldClick$X^2 + CoordinatesOldClick$Y^2)
                if (.dist >= 0.8 & .dist <= 1) {
                    onCircle$is <- TRUE 
                } else {
                    onCircle$is <- FALSE
                }
            } else onCircle$is <- FALSE
        })
        
        ## Click: which is the current sector?
        indClick <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosClick$x))
                indClick$ind <- minFragCart2Polar(input$circosClick$x,
                                                  input$circosClick$y, degreeFeatures())
        })
        
        
        
        ## double click: which is the current sector?
        CoordinatesNewDblClick <- reactiveValues(X = 0, Y = 0)
        CoordinatesOldDblClick <- reactiveValues(X = 0, Y = 0)
        
        observe({
            if (!is.null(input$circosDblClick$x)) {
                CoordinatesNewDblClick$X <- input$circosDblClick$x
                CoordinatesNewDblClick$Y <- input$circosDblClick$y
                CoordinatesOldDblClick$X <- CoordinatesNewDblClick$X
                CoordinatesOldDblClick$Y <- CoordinatesNewDblClick$Y
            } else {
                CoordinatesNewDblClick$X <- CoordinatesOldDblClick$X
                CoordinatesNewDblClick$Y <- CoordinatesOldDblClick$Y
            }
        })
        
        ## reactive value which stores double clicked indices (inds = storage, 
        ## new = new indices)
        indDblClick <- reactiveValues(ind = NULL, new = NULL)
        observe({
            if (!is.null(input$circosDblClick$x)) {
                
                minInd <- minFragCart2Polar(input$circosDblClick$x,
                                            input$circosDblClick$y,
                                            degreeFeatures())
                if (!is.na(minInd)) {
                    GNselect <- GN()[minInd]
                    selected <- strsplit(GNselect, split = "_")[[1]]
                    groupSelected <- selected[1]
                    nameSelected <- selected[2:length(selected)]
                    
                    newNG <- paste(groupSelected, nameSelected, sep = "_")
                    ## write truncated name to indDblClick$new
                    indDblClick$new <- newNG
                } else  indDblClick$new <- NULL
            } 
        })
        
        observe({
            input$resetClickIndices
            isolate(indDblClickMZ$ind <- NULL)
            isolate(indDblClickRT$ind <- NULL)
            isolate(indDblClickCluster$ind <- NULL)
            isolate(indDblClick$new <- NULL)
            isolate(onCircle$is <- FALSE)
            isolate(indClick$ind <- NULL)
        })
        
        ## reset indClick when changing radio button order
        observe({
            input$order
            isolate(onCircle$is <- FALSE)
            isolate(indClick$ind <- NULL)
        })
        
        
        ## write double-clicked (truncated) names to indDblClickMZ, indDblClickRT, 
        ## indDblClickCluster
        indDblClickMZ <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosDblClick$x)) {
                if (!is.null(indDblClick$new)) {
                    
                    newMZ <- paste(groupMZ, nameMZ, sep = "_")
                    newIndMZ <- match(indDblClick$new, newMZ)
                    
                    if (isolate(newIndMZ %in% indDblClickMZ$ind)) {
                        indDblClickMZ$ind <- isolate(indDblClickMZ$ind[-which(newIndMZ == indDblClickMZ$ind)])
                    } else {indDblClickMZ$ind <- isolate(c(indDblClickMZ$ind, newIndMZ))}
                }
            }
        })
        
        indDblClickRT <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosDblClick$x)) {
                if (!is.null(indDblClick$new)) {
                    
                    newRT <- paste(groupRT, nameRT, sep = "_")
                    newIndRT <- match(indDblClick$new, newRT)
                    
                    if (isolate(newIndRT %in% indDblClickRT$ind)) {
                        indDblClickRT$ind <- isolate(indDblClickRT$ind[-which(newIndRT == indDblClickRT$ind)])
                    } else {indDblClickRT$ind <- isolate(c(indDblClickRT$ind, newIndRT))}
                }
            }
        })
        
        indDblClickCluster <- reactiveValues(ind = NULL)
        observe({
            if (!is.null(input$circosDblClick$x)) {
                if(!is.null(indDblClick$new)) {
                    newCl <- paste(groupClust, nameClust, sep = "_")
                    newIndCl <- match(indDblClick$new, newCl)
                    
                    if (isolate(newIndCl %in% indDblClickCluster$ind)) {
                        indDblClickCluster$ind <- isolate(indDblClickCluster$ind[-which(newIndCl == indDblClickCluster$ind)])
                    } else {indDblClickCluster$ind <- isolate(c(indDblClickCluster$ind, newIndCl))}
                }
            }
        })
        
        ## plotting
        initializePlot <- reactive({
            circos.initialize(factor(GN(), levels = GN()),
                              xlim=matrix(rep(c(0,1), length(mz_match)), ncol=2, byrow=TRUE))
            circos.trackPlotRegion(factor(GN(), levels=GN()), ylim=c(0,1))  
        })
        
        ## assign to output$circos
        output$circos <- renderPlot({
            initializePlot()
            if (onCircle$is) {
                if (input$order == "mz") {
                    replayPlot(PlotHighlightMZ)
                    highlight(mz_match, c(indClick$ind, indDblClickMZ$ind), 
                              LinkDf_threshold())  
                }
                
                if (input$order == "retentionTime") {
                    replayPlot(PlotHighlightRT)
                    highlight(rt_match, c(indClick$ind, indDblClickRT$ind), 
                              LinkDf_threshold())    
                }
                
                if (input$order == "clustering") {
                    replayPlot(PlotHighlightClust)
                    highlight(clust_match, c(indClick$ind, 
                                    indDblClickCluster$ind), LinkDf_threshold())  
                }
            } else { ## if not onCircle$is
                if (length(indDblClickMZ$ind) > 0) {
                    if (input$order == "mz") {
                        replayPlot(PlotHighlightMZ)
                        highlight(mz_match, 
                                c(indDblClickMZ$ind), LinkDf_threshold()) 
                    }
                    
                    if (input$order == "retentionTime") {
                        replayPlot(PlotHighlightRT)
                        highlight(rt_match, 
                                c(indDblClickRT$ind), LinkDf_threshold()) 
                    }
                    
                    if (input$order == "clustering") {
                        replayPlot(PlotHighlightClust)
                        highlight(clust_match, 
                                c(indDblClickCluster$ind), LinkDf_threshold())  
                    }
                    
                } else {
                    if (input$order == "mz") {
                        replayPlot(PlotFilledMZ)
                        plotCircos(mz_match, LinkDf_threshold(), 
                                   initialize=FALSE, featureNames=FALSE, 
                                   groupSector=FALSE, groupName=FALSE, 
                                   links=TRUE, highlight=FALSE)
                    }
                    
                    
                    if (input$order == "retentionTime") {
                        replayPlot(PlotFilledRT)
                        plotCircos(rt_match, LinkDf_threshold(), 
                                   initialize=FALSE, featureNames = FALSE, 
                                   groupSector = FALSE, groupName = FALSE, 
                                   links = TRUE, highlight = FALSE)
                    }
                    
                    if (input$order == "clustering") {
                        replayPlot(PlotFilledClust)
                        plotCircos(clust_match, LinkDf_threshold(), 
                                   initialize=FALSE, featureNames=FALSE,
                                   groupSector=FALSE, groupName=FALSE, 
                                   links=TRUE, highlight=FALSE)
                    }
                }
            }
        })
        
        output$sized_plot <- renderUI({
            plotOutput("circos",
                       dblclick = "circosDblClick",
                       click = "circosClick",
                       width = ifelse(is.null(input$innerWidth), 0, input$innerWidth*0.5*input$plotSize), 
                       height = ifelse(is.null(input$innerWidth), 0, input$innerWidth*0.5*input$plotSize))
        })
        
        output$circosLegend <- renderPlot({
            if (!is.null(input$legend)) if(input$legend)
                circosLegend(rt_match, highlight=TRUE)
        })
        
        ## show when Clicking the feature which connects to it
        linkDfIndsClick <- reactive({
            getLinkDfIndices(GN()[indClick$ind], LinkDf_threshold())
        })
        
        output$clickConnectedFeature <- renderUI({ 
            if (!is.null(onCircle$is)) {
                if (onCircle$is)
                    HTML(printInformationSelect(select=GN()[indClick$ind], 
                                                spectra=spe$spectra, linkDfInd=linkDfIndsClick(), 
                                                linkDf=LinkDf_threshold(), 
                                                similarityMatrix=similarityMatrix, 
                                                roundDigits=input$precision))  
            }
        })
        
        output$dblClickFeature <- renderText({
            if (length(indDblClickMZ$ind) > 0) 
                c("(permanently) selected features: ", mz_match[indDblClickMZ$ind])
            else "no features (permanently) selected"
        })
        
        ## Tab Spectra
        output$querySpectra <- renderUI({
            selectInput("query", "Choose MS/MS spectra #2:", choices=GN())
        })
        
        output$subjectSpectra <- renderUI({
            selectInput("subject", "Choose MS/MS spectra #1:", choices=GN(), selected=GN()[indClick$ind])
        })
        
        output$plotSpec <- renderPlot({
            if (!is.null(input$subject) && !is.null(input$query)) {
                plotSpectra(spectra, input$subject, input$query)
            } 
        })
        
        ## on exit
        observe({
            if (input$stop == 0)
                return()
            else {
                circos.clear()
                selectedFeatures <- as.character(paste(
                    groupMZ[indDblClickMZ$ind], nameMZ[indDblClickMZ$ind], 
                    sep="_"))
                stopApp(
                    list(spectra=spe$spectra, selectedFeatures=selectedFeatures)
                )
            }
        })
        
        
    }
    
    app <- list(ui = ui, server = server)
    runApp(app)
}



#' @name printInformationSelect
#' @title Display information on connected features of selected features
#' @description Displays information on connected features of selected features.
#' @usage printInformationSelect(groupname, msp = NULL, ind, 
#'  lMatInd, linkMatrixThreshold, similarityMatrix, roundDigits = 2)
#' @param groupname \code{character} vector with groupname of selected feature,
#' vector containing "group" and "name" to display, that is 
#' a unique identifier of the features, "group" and "name" have to be separated
#' by \code{"_"} where "group" is the first and "name" is the last element
#' @param msp \code{MSP}, an S4 object of class \code{MSP} for information about 
#'  the selected feature
#' @param ind \code{numeric}
#' @param lMatInd \code{numeric} indices of selected features
#' @param linkMatrixThreshold \code{matrix} that contains information of linked 
#'  features for given thresholds
#' @param similarityMatrix \code{matrix} that is used to get information on the 
#' degree of similarity, \code{similarityMat} is an ordered version of a 
#' similarity matrix, see \code{?createOrderedSimMat}
#' @param roundDigits \code{numeric},  how many digits should be displayed?
#' @details \code{printInformationSelect} is for internal use. 
#' @return \code{character} that is in HTML format
#' @examples
#' data("idMSMStoMSP", package = "MetCirc")
#' data("binnedMSP", package = "MetCirc")
#' ## use only a selection
#' binnedMSP <- binnedMSP[c(1:20, 29:48, 113:132, 240:259),]
#' similarityMat <- createSimilarityMatrix(binnedMSP)
#' groupname <- rownames(similarityMat)
#' ## order similarityMat according to mz
#' simMat <- createOrderedSimMat(similarityMat, order = "mz") 
#' groupnameMZ <- rownames(simMat)
#' linkMat_thr <- createLinkMatrix(simMat, 0.8, 1) 
#' ind <- 2
#' indMZ <- which(groupname[ind] == truncateName(groupnameMZ, NULL, group = TRUE))
#' linkMatInds <- getLinkDfIndices(groupnameMZ[indMZ], linkMat_thr)
#' MetCirc:::printInformationSelect(groupname = groupname, 
#'  msp = NULL, ind = ind, lMatInd = linkMatInds, 
#'  linkMatrixThreshold = linkMat_thr, 
#'  similarityMatrix = similarityMat, roundDigits = 2)
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @return
printInformationSelect <- function(select=groupn[ind], spectra=NULL, 
                linkDfInd, linkDf, similarityMatrix, roundDigits = 2) {

    ## connected features: find
    connect <- as.character(unique(unlist(linkDf[linkDfInd, c("spectrum1", "spectrum2")])))
    
    ## remove duplicated hovFeat in connect
    if (select %in% connect) connect <- connect[-which(connect == select)]
    
    ## truncate group    
    connect_cut <- lapply(strsplit(connect, split="_"), "[", -1)
    connect_cut <- unlist(lapply(connect_cut, function(x) paste(x, collapse="_")))
    select_cut <- lapply(strsplit(select, split="_"), "[", -1)
    select_cut <- unlist(lapply(select_cut, function(x) paste(x, collapse="_")))
    
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
            
        connChar <- paste(connChar, collapse=" ")

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

