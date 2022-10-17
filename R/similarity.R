#' @name neutralloss
#'
#' @title Calculate similarity based on neutral losses
#'
#' @description Calculate similarity based on neutral losses (NLS)
#'
#' @param x \code{Spectra} object from \code{Spectra} containing
#' intensity and m/z values, first MS/MS spectrum
#' @param y \code{Spectra} object from \code{Spectra} containing 
#' intensity and m/z values, second MS/MS spectrum
#' @param m \code{numeric(1)}, exponent to calculate peak intensity-based 
#' weights
#' @param n \code{numeric(1)}, exponent to calculate m/z-based weights
#' @param na.rm \code{logical(1)}, if \code{NA} values should be removed
#' @param ... further arguments
#'
#' @details 
#' Similarities of spectra based on neutral losses are calculated according to 
#' the following formula: 
#' 
#' \deqn{NLS = \frac{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2}{ \sum(W_{S1, i} ^ 2) \cdot \sum(W_{S2, i} ^ 2) }}{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2 \sum(W_{S1, i} ^ 2) * \sum(W_{S2, i} ^ 2)},
#' 
#' with \eqn{W = [ peak intensity] ^{m} \cdot [ NL ]^n} and 
#' \eqn{NL = | m/z - precursor m/z |}. For further information 
#' see Li et al. (2015): Navigating natural variation in herbivory-induced
#' secondary metabolism in coyote tobacco populations using MS/MS structural 
#' analysis. PNAS, E4147--E4155. 
#' 
#' In here, the precursor m/z is taken by the m/z feature with the highest 
#' intensity. 
#' 
#' \code{neutralloss} returns a numeric value ranging between 0 and 1, where 0 
#' indicates no similarity between the two MS/MS features, while 1 indicates 
#' that the MS/MS features are identical. 
#' Prior to calculating \deqn{W_{S1}} or \deqn{W_{S2}}, all intensity values 
#' are divided by the maximum intensity value.
#'
#' @return
#' \code{neutralloss} returns a numeric similarity coefficient between 0 and 1
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples 
#' data("spectra", package = "MetCirc")
#' Spectra::compareSpectra(sps_tissue[1:10], FUN = neutralloss, m = 0.5, n = 2)
#'
#' @export
neutralloss <- function(x, y, m = 0.5, n = 2, na.rm = TRUE, ...) {

    ## calculate loss to mz value with highest intensity value
    x[, 1L] <- abs(x[, 1L] - x[which.max(x[, 2L]), 1L])
    y[, 1L] <- abs(y[, 1L] - y[which.max(y[, 2L]), 1L])
    
    ## normalize to % intensity
    x[, 2L] <- x[, 2L] / max(x[, 2L], na.rm = TRUE)*100
    y[, 2L] <- y[, 2L] / max(y[, 2L], na.rm = TRUE)*100

    wx <- MsCoreUtils:::.weightxy(x[, 1L], x[, 2L], m, n)
    wy <- MsCoreUtils:::.weightxy(y[, 1L], y[, 2L], m, n)
    # ws1 <- x[, 2L] ^ m * x[, 1L] ^ n
    # ws2 <- y[, 2L] ^ m * y[, 1L] ^ n

    sum(wx * wy, na.rm = na.rm)^2L/(sum(wx^2L, na.rm = na.rm) * 
        sum(wy^2L, na.rm = na.rm))
    # sum(ws1*ws2, na.rm = na.rm) ^ 2 / (sum(ws1^2, na.rm = na.rm) * 
    #     sum(ws2^2, na.rm = na.rm))
}

#' @name orderSimilarityMatrix
#'
#' @title
#' Order columns and rows of a similarity matrix according to 
#' m/z, retention time and clustering
#'
#' @description
#' Internal function for shiny application. May also be used 
#' outside of shiny to reconstruct figures.
#'
#' @param similarityMatrix \code{matrix}, \code{similarityMatrix} contains 
#' pair-wise similarity coefficients which give information about the similarity 
#' between precursors
#' @param sps \code{Spectra} object containing spectra that are compared
#' in \code{similarityMatrix}
#' @param type \code{character(1)}, one of \code{"retentionTime"}, \code{"mz"}, or 
#' \code{"clustering"}
#' @param group \code{logical(1)}, if \code{TRUE} \code{group} separated by 
#' \code{"_"} will be cleaved  from \code{rownames}/\code{colnames} of 
#' \code{similarityMatrix} and matched against names of \code{sps}
#' (\code{sps$name}), if \code{FALSE} \code{rownames}/\code{colnames} of 
#' \code{similarityMatrix} are taken as are and matched against names of 
#' \code{sps} (\code{sps$name})
#'
#' @details
#' \code{orderSimilarityMatrix} takes a similarity matrix,
#' \code{Spectra} object (\code{sps}, containing information on m/z and 
#' retention time), and a \code{character} vector as arguments. It will then 
#' reorder rows and columns of the \code{similarityMatrix} object such, 
#' that it orders rows and columns of \code{similarityMatrix} according to 
#' m/z, retention time or clustering in each group. 
#' 
#' \code{orderSimilarityMatrix} is employed in the \code{shinyCircos} 
#' function to create \code{similarityMatrix} objects which will allow to switch
#' between different types of ordering in between groups (sectors) in the 
#' circos plot. It may be used as well externally, to reproduce plots outside
#' of the reactive environment (see vignette for a workflow).
#'
#' @return \code{matrix}, \code{orderSimilarityMatrix} returns a similarity 
#' matrix with ordered \code{rownames} according to the \code{character} vector 
#' \code{type}
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("spectra", package = "MetCirc")
#' similarityMat <- Spectra::compareSpectra(sps_tissue[1:10], 
#'     FUN = MsCoreUtils::ndotproduct, ppm = 10)
#' rownames(similarityMat) <- colnames(similarityMat) <- sps_tissue$name[1:10]
#' 
#' ## order according to retention time 
#' orderSimilarityMatrix(similarityMatrix = similarityMat, 
#'     sps = sps_tissue, type = "retentionTime", group = FALSE)
#'     
#' @importFrom amap hcluster
#' @importFrom MsCoreUtils ndotproduct
#' 
#' @export
orderSimilarityMatrix <- function(similarityMatrix, sps, 
        type = c("retentionTime","mz", "clustering"), group = FALSE) {

    if (!(group %in% c(TRUE, FALSE))) stop("group has to be TRUE or FALSE")
    
    type <- match.arg(type)
    groupname <- rownames(similarityMatrix)
    
    ## set diagonal of similarityMatrix to 1
    diag(similarityMatrix) <- 1

    ## if a group is preciding the rownames of similarityMatrix
    if (group) {
        groupname <- strsplit(groupname, split = "_")
        groupname <- lapply(groupname, "[", -1)
        groupname <- lapply(groupname, function(x) paste(x, collapse = "_"))
        groupname <- unlist(groupname)
    }

    ## match colnames/rownames of similarityMatrix and names of sps and
    ## truncate sps
    sps <- sps[sps$name %in% groupname, ]

    ## retentionTime
    if (type == "retentionTime") {
        orderNew <- order(sps$rtime)
    }

    ## mz
    if (type == "mz") {
        orderNew <- order(sps$precursorMz)
    }

    ## clustering
    if (type == "clustering") {
        hClust <- amap::hcluster(similarityMatrix, method = "spearman")
        orderNew <- hClust$order
    }

    simM <- similarityMatrix[orderNew, orderNew]
    colnames(simM) <- rownames(simM) <- rownames(similarityMatrix)[orderNew]

    simM
}
