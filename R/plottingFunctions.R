#' @import shiny
#' @import circlize
#' @import ggplot2
#' @importFrom scales alpha

#' @name plotCircos
#'
#' @title Circular plot to visualise similarity
#'
#' @description
#' Circular plot to visualise similarity.
#'
#' @param
#' groupname `character` vector containing "group" and "name" to 
#' display, that is a unique identifier of the features, "group" and "name" have 
#' to be separated 
#' by `"_"` where "group" is the first and "name" is the last element
#'
#' @param 
#' linkDf `data.frame` containing linked features in each row, has 
#' five columns (group1, spectrum1, group2, spectrum2, similarity)
#' 
#' @param initialize `logical`, should plot be initialized?
#'
#' @param featureNames `logical`, should feature names be displayed?
#'
#' @param cexFeatureNames `numeric` size of feature names
#'
#' @param 
#' groupSector `logical`, should groups be displayed with background 
#' colours?
#'
#' @param 
#' groupName `logical`, should group names (e.g. compartment names or 
#' individual names) be displayed?
#'
#' @param links `logical`, should links be plotted?
#'
#' @param highlight `logical`, highlight is set to `TRUE`
#'
#' @param
#' colour `NULL` or `character`, colour defines the colours 
#' which are used for plotting, if `NULL` default colours are used
#'
#' @param transparency `numeric`, defines the transparency of the colours
#'
#' @details 
#' Internal use for `shinyCircos` or used outside of 
#' `shinyCircos` to reproduce figure
#'
#' @return The function will initialize a circlize plot and/or will plot 
#'  features of a circlize plot.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples 
#' data("spectra", package = "MetCirc")
#' similarityMat <- compare_Spectra(spectra_tissue[1:10],
#'     fun = normalizeddotproduct, binSize = 0.01)
#' ## order similarityMat according to retentionTime
#' simM <- orderSimilarityMatrix(similarityMat, spectra = spectra_tissue[1:10],
#'             type = "retentionTime", )
#' ## create link data.frame
#' linkDf <- createLinkDf(similarityMatrix = simM, spectra = spectra_tissue,
#'      condition = c("SPL", "LIM", "ANT", "STY"), lower = 0.5, upper = 1)
#' ## cut link data.frame (here: only display links between groups)
#' linkDf_cut <- cutLinkDf(linkDf, type = "inter")
#' ## set circlize paramters
#' circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0),
#'          track.margin = c(0.0, 0))
#' groupname <- c(as.character(linkDf_cut[, "spectrum1"]),
#'                 as.character(linkDf_cut[, "spectrum2"]))
#' groupname <- unique(groupname)
#' ## actual plotting
#' plotCircos(groupname, linkDf_cut, initialize = TRUE,
#'     featureNames = TRUE, cexFeatureNames = 0.3, groupSector = TRUE,
#'      groupName = FALSE, links = FALSE, highlight = FALSE, colour = NULL,
#'      transparency = 0.2)
#'
#' @export
plotCircos <- function(groupname, linkDf, initialize = c(TRUE, FALSE), 
    featureNames = c(TRUE, FALSE), cexFeatureNames = 0.3, 
    groupSector = c(TRUE, FALSE), groupName = c(TRUE, FALSE), 
    links = c(TRUE, FALSE), highlight = c(TRUE, FALSE), colour = NULL, 
    transparency = 0.2) {

    ## get group and name from groupname argument
    ## groupname is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    group <- lapply(strsplit(groupname, split = "_"), "[", 1)
    group <- unlist(group)
    name <- lapply(strsplit(groupname, split = "_"), function (x) x[length(x)])
    name <- unlist(name)
    
    ## get length of vector groupname
    groupname_l <- length(groupname)

    if (!is.numeric(cexFeatureNames)) stop("cexFeatureNames is not numeric")
    if (!is.logical(initialize)) stop("initialize is not logical")
    if (!is.logical(featureNames)) stop("featureNames is not logical")
    if (!is.logical(groupSector)) stop("groupSector is not logical")
    if (!is.logical(groupName)) stop("groupName is not logical")
    if (!is.logical(links)) stop("links is not logical")
    if (!is.logical(highlight)) stop("highlight is not logical")
    if (!is.null(transparency)) {
        if(!is.numeric(transparency)) stop("transparency is not numeric")
    }

    if (initialize) {
        circos.initialize(factor(groupname, levels = groupname),
            xlim = matrix(rep(c(0,1), groupname_l), ncol = 2, byrow = TRUE) )
        circos.trackPlotRegion(factor(groupname, levels = groupname), ylim = c(0,1))
    }

    ## display feature names
    if (featureNames) {
        ##groupnameFeatName <- paste(group, name, sep = "_")
        ##truncatedName <- truncateName(groupname)
        for (i in 1:groupname_l) {
            circos.text(x = 0.5, y = 0.5, labels = groupname[i], 
                sector.index = groupname[i], 
                facing = "clockwise", cex = as.numeric(cexFeatureNames),
                niceFacing = TRUE)
        }
    }

    ## create vector with unique groups
    uniqueGroup <- unique(group)

    if (!is.null(colour)) {
        if (length(colour) != length(uniqueGroup)) {
            if (length(colour) != 1) {
                stop("length of colour does not match with length of group")
            }
        }
    }

    ## group sector
    if (groupSector) {

        transparency <- if (highlight) {
            transparency - 0.1
        } else {
            transparency
        }
        
        if (is.null(colour)) {
            colour <- alpha(palette()[as.numeric(as.factor(uniqueGroup)) + 1],
                            transparency)
        } else {
            colour <- alpha(colour, transparency)
        }
        
        for( i in 1:length(uniqueGroup)) {
            ind <- which(uniqueGroup[i] == group)
            minInd <- min(ind)
            maxInd <- max(ind)
            circlize::highlight.sector(groupname[minInd:maxInd],
                col = colour[i])
        }
    }

    ## group name
    if (groupName) {
        for( i in 1:length(uniqueGroup)) {
            ind <- which(uniqueGroup[i] == group)
            minInd <- min(ind)
            maxInd <- max(ind)
            circlize::circos.text(x = 0.5, y = 1.5, labels = uniqueGroup[i],
                     sector.index = groupname[c(minInd:maxInd)[floor(length(minInd:maxInd) / 2)]],
                     facing = "downward")
        }
    }

    ## plot links
    if (links) {
        ##colourLink <- rep(alpha("black", 0.05), dim(linkDf)[1])

        if (dim(linkDf)[1] != 0) {
            for (i in 1:dim(linkDf)[1]) {
                circos.link(linkDf[i,]$"spectrum1", 0.5,
                    linkDf[i,]$"spectrum2", 0.5,
                    lwd = if (highlight) 0.3 else max(0.5, linkDf[i,][["similarity"]]),
                    ## transparency
                    col = rep(alpha("black", 0.05))) ##colourLink[i])
            }
        }
    }
}

#' @name highlight
#'
#' @title Add links and highlight sectors
#'
#' @description
#' A function to add links and highlight sectors to an initialised
#' and plotted `circlize` plot with one track.
#'
#' @param
#' groupname `character` vector containing "group" and "name" to 
#' display, that is a unique identifier of the features, "group" and "name" have 
#' to be separated by `"_"` where "group" is the first and "name" is the 
#' last element
#'
#' @param
#' ind `numeric`, indices which will be highlighted
#'
#' @param
#' linkDf `data.frame`, in each row there is information about 
#' features to be connected
#'
#' @param
#' colour `NULL` or `character`, colour defines the colours which 
#' are used for plotting, if `NULL` default colours are used
#'
#' @param
#' transparency `numeric`, defines the transparency of the colours
#'
#' @param
#' links `logical`, should links of unselected features be plotted
#'
#' @details
#' Internal use for `shinyCircos` or outside of `shinyCircos` to reproduce the 
#' figure.
#'
#' @return 
#' The function will update an existing plot by highlighting a 
#' specified sector and connected links.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples 
#' data("spectra", package = "MetCirc")
#' similarityMat <- compare_Spectra(spectra_tissue[1:10], 
#'     fun = normalizeddotproduct, binSize = 0.01)
#'  ## order similarityMat according to retentionTime and update rownames
#'  simM <- orderSimilarityMatrix(similarityMat, spectra = spectra_tissue[1:10], 
#'              type = "retentionTime")
#'  ## create link matrix
#'  linkDf <- createLinkDf(similarityMatrix = simM, spectra = spectra_tissue,
#'      condition = c("SPL", "LIM", "ANT", "STY"), lower = 0.5, upper = 1)
#'  ## cut link matrix (here: only display links between groups)
#'  linkDf_cut <- cutLinkDf(linkDf, type = "inter")
#'  ## set circlize parameters
#'  circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), 
#'          track.margin = c(0.0, 0))
#'  groupname <- c(as.character(linkDf_cut[, "spectrum1"]), 
#'                  as.character(linkDf_cut[, "spectrum2"]))
#'  groupname <- unique(groupname)
#'  ## here: set indSelected arbitrarily
#'  indSelected <- c(2,3)
#'  ## actual plotting
#'  plotCircos(groupname, linkDf_cut, initialize = TRUE, 
#'      featureNames = TRUE, cexFeatureNames = 0.2, groupSector = TRUE, 
#'      groupName = FALSE, links = FALSE, highlight = TRUE)
#'  ## highlight
#'  highlight(groupname = groupname, ind = indSelected, linkDf = linkDf_cut, 
#'      colour = NULL, transparency = 0.4, links = TRUE)
#'
#' @export
highlight <- function(groupname, ind, linkDf, colour = NULL, transparency = 0.4, 
                    links = TRUE) {

    ## get group and name from groupname argument
    ## groupname is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    group <- lapply(strsplit(groupname, split = "_"), "[", 1)
    group <- unlist(group)
    name <- lapply(strsplit(groupname, split = "_"), function (x) x[length(x)])
    name <- unlist(name)

    ##if (length(colour))
    ## get length of vector namegroup
    groupname_l <- length(groupname)

    ## create vector that contains selected (ind) groupname instances
    groupnameselected <- groupname[ind]
    nameselected <- name[ind]

    ## retrieve spectrum1 and spectrum2 from linkDf
    linkDfSpec1 <- linkDf[,"spectrum1"]
    linkDfSpec2 <- linkDf[,"spectrum2"]

    if (is.null(colour)) {
        colours <- alpha(palette()[as.numeric(as.factor(group))[ind]+1], transparency)
    } else {
        colours <- alpha(colour[as.numeric(as.factor(group))[ind]], transparency)
    }

    for (h in 1:length(ind)) {
        highlight.sector(sector.index = factor(groupnameselected[h],
            levels = groupnameselected[h]), col = colours[h])
    }

    ## get indices in linkDf of selected features
    if (nrow(linkDf) !=  0) {
        linkDfInd <- getLinkDfIndices(groupnameselected, linkDf)
    } else {linkDfInd <- NULL}

    ## plot all links
    if (links) {
        if (nrow(linkDf) != 0) {
            for (i in 1:dim(linkDf)[1]) {
                circos.link(linkDfSpec1[i], 0.5, linkDfSpec2[i], 0.5, lwd = 0.3,
                    col = alpha("black", 0.1))
            }
        }
    }

    ## plot highlighted links
    if (!is.null(linkDfInd)) {
        for (j in linkDfInd) {
            circos.link(linkDfSpec1[j], 0.5, linkDfSpec2[j], 0.5, lwd = 1,
                col = "black")
        }
    }
}

#' @name circosLegend
#' 
#' @title Plot a legend for circos plot
#' 
#' @description
#' `circosLegend` plots a legend for circos plot using group names.
#'
#' @param
#' groupname `character` vector containing "group" and "name" to 
#' display, that is  a unique identifier of the features, "group" and "name" have 
#' to be separated by `"_"` where "group" is the first and "name" is the 
#' last element
#'
#' @param
#' highlight `logical`, should colours be adjusted to highlight settings?
#'
#' @param
#' colour `NULL` or `character`, colour defines the colours which are 
#' used for plotting, if `NULL` default colours are used
#'
#' @param
#' cex `numeric`, parameter that controls size of the legend in the plot
#'
#' @details
#' Internal use in `shinyCircos` or outside of `shinyCircos` 
#' to reproduce figures.
#'
#' @return
#' The function will open a new plot and display colours together with labels.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
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
#' ## plot legend
#' circosLegend(groupname, highlight = TRUE, colour = NULL, cex = 1)
#'
#' @export
circosLegend <- function(groupname, highlight = TRUE, colour = NULL, cex = 1) {
    
    ## get group and name from groupname argument
    ## groupname is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    group <- lapply(strsplit(groupname, split = "_"), "[", 1)
    group <- unlist(group)
    group <- as.factor(group)
    
    uniqNumGroup <- unique(as.numeric(group))
    
    if (is.null(colour)) {
        colours <- palette()[uniqNumGroup + 1]
    } else {
        colours <- colour[uniqNumGroup + 1]
    }
    
    plot(x = c(0,0.5), y = c(0,0.5), type = "n", xlab = "", ylab = "",
         axes = FALSE, frame.plot = FALSE)
    if (highlight) {
        plot.new()
        leg <- legend(x = par("usr")[1:2], y =  par("usr")[3:4],
                legend = levels(group), bty = "n", fill = alpha(colours, 0.3),
                border = alpha(colours, 0.3), cex = 1)
    } else { ## if not highlight
        legend(x = par("usr")[1:2], y = par("usr")[3:4], legend = levels(group),
                bty = "n", fill = colours, border = colours, cex = cex)
    }
}

#' @name getLinkDfIndices
#'
#' @title Get indices in linkDf of feature
#'
#' @description Gets indices in linkDf of feature
#'
#' @param 
#' groupnameselected `character` vector with groupname of selected 
#' feature, vector containing "group" and "name" to display, that is 
#' a unique identifier of the features, "group" and "name" have to be separated
#' by `"_"` where "group" is the first and "name" is the last element
#'
#' @param 
#' linkDf `data.frame`, in each row there is information about 
#' features to be connected
#'
#' @details 
#' Internal use for function `highlight`
#'
#' @return 
#' `getLinkDfIndices` returns indices concerning `linkDf` to which 
#' `groupnameselected` connects
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples 
#' \dontrun{getLinkDfIndices(groupnameselected, linkMatrix)}
#'
#' @export
getLinkDfIndices <- function(groupnameselected, linkDf) {
    
    linkDfInd <- lapply(as.character(groupnameselected),
                            function(x) which(linkDf == x, arr.ind = TRUE))
    ## select only first column
    linkDfInd <- lapply(linkDfInd, function(x) x[,1])
    linkDfInd <- unlist(linkDfInd)
    
    linkDfInd <- as.numeric(linkDfInd)
    
    return(linkDfInd)
}

#' @name minFragCart2Polar
#'
#' @title Calculate the nearest feature in polar coordinates given cartesian
#' coordinates
#'
#' @description 
#' Calculates the nearest feature in polar coordinates given 
#' cartesian coordinates.
#'
#' @param x cartesian x coordinate
#'
#' @param y cartesian y coordinate
#'
#' @param degreeOfFeatures `list` of positions of features
#'
#' @details 
#' `minFragCart2Polar` is employed to find the feature with 
#' the smallest distance from given cartesian coordinates.
#'
#' @return 
#' `minFragCart2Polar` returns the index of the feature that has the
#' smallest distance to the given coordinates. As `minFragCart2Polar` is 
#' used in `shinyCircos` for the track 1 only polar r coordinates between
#' 0.8 and 1 will be used to find the feature with smallest distance.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
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
#'                 as.character(linkDf_cut[, "spectrum2"]))
#' groupname <- unique(groupname)
#' plotCircos(groupname, NULL, initialize = TRUE, featureNames = FALSE, 
#'     groupName = FALSE, groupSector = FALSE, links = FALSE, highlight = FALSE)
#' x <- 1
#' y <- 0
#' degreeFeatures <- lapply(groupname, 
#'  function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
#' minFragCart2Polar(x, y, degreeOfFeatures = degreeFeatures)
#'
#' @export
minFragCart2Polar <- function(x, y, degreeOfFeatures) {
    polar <- cart2Polar(x, y)
    minInd <- NA
    if (polar$r <= 1 & polar$r >= 0.8)
        minInd <- which.min(abs(polar$theta - unlist(degreeOfFeatures)))
    return(minInd)
}

#' @name cart2Polar
#'
#' @title Calculate polar coordinates from cartesian coordinates
#'
#' @description 
#' `cart2Polar` calculates polar coordinates from cartesian coordinates.
#'
#' @param x cartesian x coordinate
#'
#' @param y cartesian y coordinate
#'
#' @details
#' `cart2Polar` is employed to translate cartesian coordinates 
#'  into polar coordinates especially in interactive shiny applications when
#'  using hovering and clicking features.
#'
#' @return 
#' `cart2Polar` returns a list of colar coordinates r and theta
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples 
#' x <- 1; y <- 1
#' cart2Polar(x, y)
#'
#' @export
cart2Polar <- function(x, y) {
    r <- sqrt( x ^ 2 + y ^ 2)
    thetaP <- atan( y/x ) * 180 / pi
    if (x == 0 & y == 0) thetaP <- 0
    if (x >= 0 & y >= 0) theta <- thetaP ## 1st quadrant
    if (x < 0 & y >= 0) theta <- thetaP + 180 ## 2nd quadrant
    if (x < 0 & y < 0) theta <- thetaP + 180 ## 3rd quadrant
    if (x >= 0 & y < 0) theta <- thetaP + 360 ## 4th quadrant
    
    return(list(r = r, theta = theta))
}

#' @name plotSpectra
#'
#' @title Plot pair-wise spectra
#'
#' @description `plotSpectra` plots a spectra of a `subject` and `query`
#' spectra. `plotSpectra` uses `ggplot` plotting functionality.
#'
#' @param spectra `Spectra`` object
#'
#' @param
#' subject character, name of spectra that is aligned against, character
#' with preceding sample name
#'
#' @param
#' query character, name of spectra that is aligned to subject, character
#' with preceding sample name
#'
#' @details Internally, all intensities are normalized to 100\%.
#'
#' @return `ggplot2` plot
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples 
#' data("spectra", package = "MetCirc")
#' plotSpectra(spectra_tissue, subject = "SPL_1", query = "SPL_2")
#'
#' @export
plotSpectra <- function(spectra, subject, query) {

    ## strsplit subject and spectra to remove sample name
    subject <- strsplit(subject, split = "_")
    subject <- paste(subject[[1]][-1], collapse = "_")
    query <- strsplit(query, split = "_")
    query <- paste(query[[1]][-1], collapse = "_")

    mz_sub <- spectra[subject]@listData[[1]]@mz
    int_sub <- spectra[subject]@listData[[1]]@intensity
    mz_que <- spectra[query]@listData[[1]]@mz
    int_que <- spectra[query]@listData[[1]]@intensity

    ## normalize to 100% 
    int_sub <- int_sub / max(int_sub) * 100
    int_que <- int_que / max(int_que) * 100

    df_sub <- data.frame(mz = mz_sub, int = -int_sub, is = "MS/MS #1")
    df_que <- data.frame(mz = mz_que, int = int_que, is = "MS/MS #2")

    df <- rbind(df_que, df_sub)

    ggplot(df) + 
        geom_segment(aes(x = mz, xend = mz+0.001, y = int, yend = 0, col = is),
                        stat = "identity") +
        xlab("m/z") +
        scale_y_continuous("intensity (%)", breaks = c(-100, -50, 0, 50, 100),
                        labels = c("100", "50", "0", "50", "100")) +
        theme_light() + theme(legend.title = element_blank())
}