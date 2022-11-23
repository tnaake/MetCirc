#' @name plotCircos
#'
#' @title Circular plot to visualize similarity
#'
#' @description
#' Circular plot to visualize similarity.
#'
#' @param groupname \code{character} vector containing \code{"group"} and 
#' \code{"name"} to display that is a unique identifier of the features, 
#' \code{"group"} and \code{"name"} have to be separated
#' by \code{"_"} where \code{"group"} is the first and \code{"name"} is the 
#' last element
#'
#' @param linkDf \code{data.frame} containing linked features in each row, has
#' five columns (\code{group1}, \code{spectrum1}, \code{group2}, 
#' \code{spectrum2}, \code{similarity})
#' @param initialize \code{logical}, should plot be initialized?
#' @param featureNames \code{logical}, should feature names be displayed?
#' @param cexFeatureNames \code{numeric} size of feature names
#' @param groupSector \code{logical}, should groups be displayed with background
#' colours?
#' @param groupName \code{logical}, should group names (e.g. compartment names 
#' or individual names) be displayed?
#' @param links \code{logical}, should links be plotted?
#' @param highlight \code{logical}, \code{highlight} is set to \code{TRUE} 
#' by default
#' @param colour \code{NULL} or \code{character}, \code{colour} defines the 
#' colours which are used for plotting, if \code{NULL} default colours are used
#'
#' @param transparency \code{numeric}, defines the transparency of the colours
#'
#' @details
#' Internal use for \code{shinyCircos} or used outside of
#' \code{shinyCircos} to reproduce figure
#'
#' @return 
#' The function will initialize a circlize plot and/or will plot
#' features of a circlize plot.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' library("MsCoreUtils")
#' data("spectra", package = "MetCirc")
#' 
#' ## create similarity matrix
#' similarityMat <- Spectra::compareSpectra(sps_tissue[1:10],
#'     FUN = MsCoreUtils::ndotproduct, ppm = 20, m = 0.5, n = 2)
#' rownames(similarityMat) <- colnames(similarityMat) <- sps_tissue$name[1:10]
#' 
#' ## order similarityMat according to retentionTime
#' simM <- orderSimilarityMatrix(similarityMat, sps = sps_tissue[1:10],
#'             type = "retentionTime")
#'             
#' ## create link data.frame
#' linkDf <- createLinkDf(similarityMatrix = simM, sps = sps_tissue,
#'      condition = c("SPL", "LIM", "ANT", "STY"), lower = 0.01, upper = 1)
#' ## cut link data.frame (here: only display links between groups)
#' linkDf_cut <- cutLinkDf(linkDf, type = "inter")
#' 
#' ## set circlize paramters
#' circos.clear()
#' circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0),
#'          track.margin = c(0.0, 0))
#' groupname <- c(as.character(linkDf_cut[, "spectrum1"]),
#'                 as.character(linkDf_cut[, "spectrum2"]))
#' groupname <- unique(groupname)
#' 
#' ## actual plotting
#' plotCircos(groupname, linkDf_cut, initialize = TRUE,
#'     featureNames = TRUE, cexFeatureNames = 0.3, groupSector = TRUE,
#'     groupName = FALSE, links = FALSE, highlight = FALSE, colour = NULL,
#'     transparency = 0.2)
#'
#' @importFrom circlize circos.link circos.initialize circos.trackPlotRegion
#' @importFrom circlize circos.text circos.par circos.clear
#' @importFrom grDevices palette
#' @importFrom MsCoreUtils ndotproduct
#' @importFrom scales alpha
#' @importFrom Spectra compareSpectra
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
    name <- lapply(strsplit(groupname, split = "_"), function(x) x[length(x)])
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
        if (!is.numeric(transparency)) stop("transparency is not numeric")
    }

    if (initialize) {
        circlize::circos.initialize(factor(groupname, levels = groupname),
            xlim = matrix(rep(c(0, 1), groupname_l), ncol = 2, byrow = TRUE))
        circlize::circos.trackPlotRegion(
            factor(groupname, levels = groupname), ylim = c(0, 1))
    }

    ## display feature names
    if (featureNames) {
    
        for (i in 1:groupname_l) {
            circlize::circos.text(x = 0.5, y = 0.5, labels = groupname[i],
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
            colour <- scales::alpha(
                grDevices::palette()[as.numeric(as.factor(uniqueGroup)) + 1],
                transparency)
        } else {
            colour <- scales::alpha(colour, transparency)
        }

        for (i in seq_len(length(uniqueGroup))) {
            ind <- which(uniqueGroup[i] == group)
            minInd <- min(ind)
            maxInd <- max(ind)
            circlize::highlight.sector(groupname[minInd:maxInd],
                col = colour[i])
        }
    }

    ## group name
    if (groupName) {
        for( i in seq_len(length(uniqueGroup))) {
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

        if (dim(linkDf)[1] != 0) {
            for (i in 1:dim(linkDf)[1]) {
                circlize::circos.link(linkDf[i, ]$"spectrum1", 0.5,
                    linkDf[i, ]$"spectrum2", 0.5,
                    lwd = if (highlight) {
                        0.3
                    } else {
                        max(0.5, linkDf[i,][["similarity"]])},
                    ## transparency
                    col = rep(alpha("black", 0.05)))
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
#' and plotted \code{circlize} plot with one track.
#'
#' @param groupname \code{character} vector containing "group" and "name" to
#' display that is a unique identifier of the features, "group" and "name" have
#' to be separated by \code{"_"} where "group" is the first and "name" is the
#' last element
#' @param ind \code{numeric}, indices which will be highlighted
#' @param linkDf \code{data.frame}, in each row there is information about
#' features to be connected
#' @param colour \code{NULL} or \code{character}, \code{colour} defines the 
#' colours which are used for plotting, if `NULL` default colours are used
#' @param transparency \code{numeric}, defines the transparency of the colours
#' @param links \code{logical}, should links of unselected features be plotted
#'
#' @details
#' Internal use for \code{shinyCircos} or outside of \code{shinyCircos} to 
#' reproduce the figure.
#'
#' @return
#' The function will update an existing plot by highlighting a
#' specified sector and connected links.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("spectra", package = "MetCirc")
#' 
#' ## create similarity matrix
#' similarityMat <- Spectra::compareSpectra(sps_tissue[1:10],
#'     FUN = MsCoreUtils::ndotproduct, ppm = 20, m = 0.5, n = 2)
#' rownames(similarityMat) <- colnames(similarityMat) <- sps_tissue$name[1:10]
#' 
#' ## order similarityMat according to retentionTime and update rownames
#' simM <- orderSimilarityMatrix(similarityMat, sps = sps_tissue[1:10],
#'     type = "retentionTime")
#' 
#' ## create link matrix
#' linkDf <- createLinkDf(similarityMatrix = simM, sps = sps_tissue,
#'     condition = c("SPL", "LIM", "ANT", "STY"), lower = 0.01, upper = 1)
#' 
#' ## cut link matrix (here: only display links between groups)
#' linkDf_cut <- cutLinkDf(linkDf, type = "inter")
#' 
#' ## set circlize parameters
#' circos.clear()
#' circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0),
#'     track.margin = c(0.0, 0))
#' groupname <- c(as.character(linkDf_cut[, "spectrum1"]),
#'     as.character(linkDf_cut[, "spectrum2"]))
#' groupname <- unique(groupname)
#' 
#' ## here: set indSelected arbitrarily
#' indSelected <- c(2,3)
#' 
#' ## actual plotting
#' plotCircos(groupname, linkDf_cut, initialize = TRUE,
#'     featureNames = TRUE, cexFeatureNames = 0.2, groupSector = TRUE,
#'     groupName = FALSE, links = FALSE, highlight = TRUE)
#' 
#' ## highlight
#' highlight(groupname = groupname, ind = indSelected, linkDf = linkDf_cut,
#'     colour = NULL, transparency = 0.4, links = TRUE)
#' 
#' @importFrom circlize circos.link highlight.sector
#' @importFrom grDevices palette
#' @importFrom MsCoreUtils ndotproduct
#' @importFrom scales alpha
#' @importFrom Spectra compareSpectra
#' 
#' @export
highlight <- function(groupname, ind, linkDf, colour = NULL, transparency = 0.4,
    links = TRUE) {

    ## get group and name from groupname argument
    ## groupname is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    group <- lapply(strsplit(groupname, split = "_"), "[", 1) |>
        unlist()
    name <- lapply(strsplit(groupname, split = "_"), function (x) x[length(x)]) |>
        unlist()

    ## get length of vector namegroup
    groupname_l <- length(groupname)

    ## create vector that contains selected (ind) groupname instances
    groupnameselected <- groupname[ind]
    nameselected <- name[ind]

    ## retrieve spectrum1 and spectrum2 from linkDf
    linkDfSpec1 <- linkDf[, "spectrum1"]
    linkDfSpec2 <- linkDf[, "spectrum2"]

    if (is.null(colour)) {
        colours <- scales::alpha(
            grDevices::palette()[as.numeric(as.factor(group))[ind] + 1], transparency)
    } else {
        colours <- scales::alpha(
            colour[as.numeric(as.factor(group))[ind]], transparency)
    }

    for (h in seq_len(length(ind))) {
        circlize::highlight.sector(sector.index = factor(groupnameselected[h],
            levels = groupnameselected[h]), col = colours[h])
    }

    ## get indices in linkDf of selected features
    if (nrow(linkDf) !=  0) {
        linkDfInd <- getLinkDfIndices(groupnameselected, linkDf)
    } else {
        linkDfInd <- NULL
    }

    ## plot all links
    if (links) {
        if (nrow(linkDf) != 0) {
            for (i in 1:dim(linkDf)[1]) {
                circlize::circos.link(linkDfSpec1[i], 0.5, 
                    linkDfSpec2[i], 0.5, lwd = 0.3,
                    col = scales::alpha("black", 0.1))
            }
        }
    }

    ## plot highlighted links
    if (!is.null(linkDfInd)) {
        for (j in linkDfInd) {
            circlize::circos.link(linkDfSpec1[j], 0.5, linkDfSpec2[j], 0.5, lwd = 1,
                col = "black")
        }
    }
}

#' @name circosLegend
#'
#' @title Plot a legend for circos plot
#'
#' @description
#' \code{circosLegend} plots a legend for circos plot using group names.
#'
#' @param groupname \code{character} vector containing "group" and "name" to 
#' display that is  a unique identifier of the features, "group" and "name" have 
#' to be separated by \code{"_"} where "group" is the first and "name" is the
#' last element
#' @param highlight \code{logical}, should colours be adjusted to highlight 
#' settings?
#' @param colour \code{NULL} or \code{character}, \code{colour} defines the 
#' colours which are used for plotting, if \code{NULL} default colours are used
#' @param cex \code{numeric}, parameter that controls size of the legend in the 
#' plot
#'
#' @details
#' Internal use in \code{shinyCircos} or outside of \code{shinyCircos}
#' to reproduce figures.
#'
#' @return
#' The function will open a new plot and display colours together with labels.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("spectra", package = "MetCirc")
#' 
#' ## create similarity matrix
#' similarityMat <- Spectra::compareSpectra(sps_tissue[1:10],
#'     FUN = MsCoreUtils::ndotproduct, ppm = 20, m = 0.5, n = 2)
#' rownames(similarityMat) <- colnames(similarityMat) <- sps_tissue$name[1:10]
#' 
#' linkDf <- createLinkDf(similarityMatrix = similarityMat,
#'     sps = sps_tissue[1:10],
#'     condition = c("SPL", "LIM", "ANT", "STY"), lower = 0.01, upper = 1)
#' 
#' ## cut link data.frame (here: only display links between groups)
#' linkDf_cut <- cutLinkDf(linkDf, type = "inter")
#' groupname <- c(as.character(linkDf_cut[, "spectrum1"]),
#'             as.character(linkDf_cut[, "spectrum2"]))
#' groupname <- unique(groupname)
#' 
#' ## plot legend
#' circosLegend(groupname, highlight = TRUE, colour = NULL, cex = 1)
#' 
#' @importFrom graphics legend par
#' @importFrom grDevices palette
#' @importFrom MsCoreUtils ndotproduct
#' @importFrom scales alpha
#' @importFrom Spectra compareSpectra
#' 
#' @export
circosLegend <- function(groupname, highlight = TRUE, colour = NULL, cex = 1) {

    ## get group and name from groupname argument
    ## groupname is a vector containing information about group and name,
    ## where group is the first element and name the last element separated by _
    group <- lapply(strsplit(groupname, split = "_"), "[", 1) |>
        unlist() |>
        as.factor()

    uniqNumGroup <- as.numeric(group) |>
        unique()
    

    if (is.null(colour)) {
        colours <- grDevices::palette()[uniqNumGroup + 1]
    } else {
        colours <- colour[uniqNumGroup + 1]
    }

    plot(x = c(0, 0.5), y = c(0, 0.5), type = "n", xlab = "", ylab = "",
         axes = FALSE, frame.plot = FALSE)
    if (highlight) {
        plot.new()
        leg <- graphics::legend(
            x = graphics::par("usr")[1:2], 
            y = graphics::par("usr")[3:4],
            legend = levels(group), bty = "n", 
            fill = scales::alpha(colours, 0.3),
            border = scales::alpha(colours, 0.3), cex = 1)
    } else { ## if not highlight
        graphics::legend(
            x = graphics::par("usr")[1:2], 
            y = graphics::par("usr")[3:4], 
            legend = levels(group), bty = "n", fill = colours, 
            border = colours, cex = cex)
    }
}

#' @name getLinkDfIndices
#'
#' @title Get indices in \code{linkDf} of feature
#'
#' @description Gets indices in \code{linkDf} of feature
#'
#' @param groupnameselected \code{character} vector with \code{groupname} of 
#' selected feature, vector containing "group" and "name" to display, that is
#' a unique identifier of the features, "group" and "name" have to be separated
#' by \code{"_"} where "group" is the first and "name" is the last element
#' @param linkDf \code{data.frame}, in each row there is information about
#' features to be connected
#'
#' @details
#' Internal use for function \code{highlight}
#'
#' @return
#' \code{getLinkDfIndices} returns indices concerning \code{linkDf} to which 
#' \code{groupnameselected} connects
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
    
    ## select only first column, unlist and return
    lapply(linkDfInd, function(x) x[,1]) |>
        unlist() |>
        as.numeric()

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
#' @param x \code{numeric}, cartesian x coordinate
#' @param y \code{numeric}, cartesian y coordinate
#' @param degreeOfFeatures \code{list} of positions of features
#'
#' @details
#' \code{minFragCart2Polar} is employed to find the feature with
#' the smallest distance from given cartesian coordinates.
#'
#' @return
#' \code{minFragCart2Polar} returns the index of the feature that has the
#' smallest distance to the given coordinates. As \code{minFragCart2Polar} is
#' used in \code{shinyCircos} for the track 1 only polar \code{r} coordinates 
#' between 0.8 and 1 will be used to find the feature with smallest distance.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' library("MsCoreUtils")
#' data("spectra", package = "MetCirc")
#' 
#' ## create similarity matrix
#' similarityMat <- Spectra::compareSpectra(sps_tissue[1:10],
#'     FUN = MsCoreUtils::ndotproduct, ppm = 20, m = 0.5, n = 2)
#' rownames(similarityMat) <- colnames(similarityMat) <- sps_tissue$name[1:10]
#' 
#' linkDf <- createLinkDf(similarityMatrix = similarityMat,
#'     sps = sps_tissue[1:10],
#'     condition = c("SPL", "LIM", "ANT", "STY"), lower = 0.5, upper = 1)
#' 
#' ## cut link data.frame (here: only display links between groups)
#' linkDf_cut <- cutLinkDf(linkDf, type = "inter")
#' groupname <- c(as.character(linkDf_cut[, "spectrum1"]),
#'                 as.character(linkDf_cut[, "spectrum2"]))
#' groupname <- unique(groupname)
#'
#' ## set circlize parameters
#' circos.clear()
#' circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0),
#'     track.margin = c(0.0, 0))
#' plotCircos(groupname, NULL, initialize = TRUE, featureNames = FALSE,
#'     groupName = FALSE, groupSector = FALSE, links = FALSE, highlight = FALSE)
#' x <- 1
#' y <- 0
#' degreeFeatures <- lapply(groupname,
#'     function(x) 
#'         mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
#' minFragCart2Polar(x, y, degreeOfFeatures = degreeFeatures)
#'
#' @importFrom circlize circos.clear circos.par 
#' @importFrom MsCoreUtils ndotproduct
#' @importFrom Spectra compareSpectra
#'
#' @export
minFragCart2Polar <- function(x, y, degreeOfFeatures) {
    
    ## convert cartesian coordinates to polar
    polar <- cart2Polar(x, y)
    
    ## determine the closest index and return
    minInd <- NA
    if (polar$r <= 1 & polar$r >= 0.8)
        minInd <- which.min(abs(polar$theta - unlist(degreeOfFeatures)))
    minInd
}

#' @name cart2Polar
#'
#' @title Calculate polar coordinates from cartesian coordinates
#'
#' @description
#' \code{cart2Polar} calculates polar coordinates from cartesian coordinates.
#'
#' @param x \code{numeric} cartesian x coordinate
#' @param y \code{numeric} cartesian y coordinate
#'
#' @details
#' \code{cart2Polar} is employed to translate cartesian coordinates
#' into polar coordinates especially in interactive shiny applications when
#' using hovering and clicking features.
#'
#' @return
#' \code{cart2Polar} returns a \code{list} of colar coordinates r and theta
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
    thetaP <- atan(y / x) * 180 / pi
    if (x == 0 & y == 0) thetaP <- 0
    if (x >= 0 & y >= 0) theta <- thetaP ## 1st quadrant
    if (x < 0 & y >= 0) theta <- thetaP + 180 ## 2nd quadrant
    if (x < 0 & y < 0) theta <- thetaP + 180 ## 3rd quadrant
    if (x >= 0 & y < 0) theta <- thetaP + 360 ## 4th quadrant

    list(r = r, theta = theta)
}

#' @name plotSpectra
#'
#' @title Plot pair-wise spectra
#'
#' @description \code{plotSpectra} plots a spectra of a \code{subject} and 
#' \code{query} spectra. \code{plotSpectra} uses \code{ggplot} 
#' plotting functionality.
#'
#' @param sps \code{Spectra} object
#' @param subject \code{character}, name of spectra that is aligned against, 
#' \code{character} with preceding sample name
#' @param query \code{character}, name of spectra that is aligned to 
#' \code{subject}, \code{character} with preceding sample name
#'
#' @details Internally, all intensities are normalized to 100\%.
#'
#' @return \code{ggplot2} plot
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples 
#' data("spectra", package = "MetCirc")
#' plotSpectra(sps = sps_tissue, subject = "SPL_1", query = "SPL_2")
#' 
#' @importFrom ggplot2 ggplot geom_segment aes xlab scale_y_continuous 
#' @importFrom ggplot2 theme_light theme element_blank
#' 
#' @export
plotSpectra <- function(sps, subject, query) {

    ## strsplit subject and spectra to remove sample name
    subject <- strsplit(subject, split = "_")[[1]]
    subject <- paste(subject[-1], collapse = "_")
    query <- strsplit(query, split = "_")[[1]]
    query <- paste(query[-1], collapse = "_")

    mz_sub <- unlist(sps[sps$name %in% subject, ]$mz)
    int_sub <- unlist(sps[sps$name %in% subject, ]$intensity)
    mz_que <- unlist(sps[sps$name %in% query, ]$mz)
    int_que <- unlist(sps[sps$name %in% query, ]$intensity)

    ## normalize to 100%
    int_sub <- int_sub / max(int_sub, na.rm = TRUE) * 100
    int_que <- int_que / max(int_que, na.rm = TRUE) * 100

    df_sub <- data.frame(mz = mz_sub, mz_add = mz_sub + 0.001, 
        int = -int_sub, is = "MS/MS #1")
    df_que <- data.frame(mz = mz_que, mz_add = mz_que + 0.001, 
        int = int_que, is = "MS/MS #2")

    df <- rbind(df_que, df_sub)

    ggplot2::ggplot(df) +
        ggplot2::geom_segment(
            ggplot2::aes(x = !!ggplot2::sym("mz"), 
                xend = !!ggplot2::sym("mz_add"), y = !!ggplot2::sym("int"), 
                yend = 0, col = !!ggplot2::sym("is")),
            stat = "identity") +
        ggplot2::xlab("m/z") +
        ggplot2::scale_y_continuous("intensity (%)", 
            breaks = c(-100, -50, 0, 50, 100),
            labels = c("100", "50", "0", "50", "100")) +
        ggplot2::theme_light() + 
        ggplot2::theme(legend.title = ggplot2::element_blank())
}
