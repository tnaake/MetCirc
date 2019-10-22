#' @importFrom MSnbase Spectra compareSpectra intensity
#' @importMethodsFrom MSnbase Spectrum2 Spectra
#' @importFrom S4Vectors DataFrame
#' @importFrom utils combn
#'
#' @name compare_Spectra
#'
#' @title Create similarity matrix from `MSnbase::Spectra` object
#'
#' @description
#' `compare_Spectra`` creates a similarity matrix of all Spectrum objects in 
#' `object`
#'
#' @param object `Spectra`
#'
#' @param 
#' fun `function` or `character`, see `?MSnbase::compareSpectra` for further
#' information
#'
#' @param ... arguments passed to `compareSpectra`
#'
#' @details
#' Function inspired by `compareSpectra.OnDiskMSnExp`. Possibly
#' transfer to `MSnbase`.
#'``
#' @author Thomas Naake (inspired by `compareSpectra.OnDiskMSnExp`)
#'
#' @examples 
#' data("spectra", package = "MetCirc")
#' compare_Spectra(spectra_tissue[1:10], fun = "dotproduct")
#'
#' @export
compare_Spectra <- function(object, fun, ...) {

    nm <- names(object)
    cb <- combn(nm, 2)
    cb <- apply(cb, 2, function(x) compareSpectra(object[[x[1]]],
                    object[[x[[2]]]], fun = fun, ...)) ## "dotproduct"

    m <- matrix(NA, length(object), length(object),
                dimnames = list(nm, nm))

    ## fill lower triangle of the matrix
    m[lower.tri(m)] <- cb

    ## copy to upper triangle
    for (i in 1:nrow(m)) {
        m[i, ] <- m[, i]
    }

    diag(m) <- 1

    return(m)
}


#' @name normalizeddotproduct
#'
#' @title Calculate the normalized dot product
#'
#' @description Calculate the normalized dot product (NDP)
#'
#' @param 
#' x `Spectrum2` object from `MSnbase` containing intensity and m/z values, 
#' first MS/MS spectrum
#'
#' @param 
#' y `Spectrum2` object from `MSnbase` containing intensity and m/z values, 
#' second MS/MS spectrum
#'
#' @param m `numeric`, exponent to calculate peak intensity-based weights
#'
#' @param n `numeric`, exponent to calculate m/z-based weights
#'
#' @param ... further arguments passed to MSnbase:::bin_Spectra
#'
#' @details
#' The normalized dot product is calculated according to the 
#' following formula: 
#' \deqn{NDP = \frac{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2}{ \sum(W_{S1, i} ^ 2) * \sum(W_{S2, i} ^ 2) }}{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2 \sum(W_{S1, i} ^ 2) * \sum(W_{S2, i} ^ 2)},
#' with \eqn{W = [ peak intensity] ^{m} \cdot [m/z]^n}. For further information 
#' see Li et al. (2015): Navigating natural variation in herbivory-induced
#' secondary metabolism in coyote tobacco populations using MS/MS structural 
#' analysis. PNAS, E4147--E4155. `normalizeddotproduct} returns a numeric 
#' value ranging between 0 and 1, where 0 
#' indicates no similarity between the two MS/MS features, while 1 indicates 
#' that the MS/MS features are identical. Arguments can be passed to 
#' the function `MSnbase:::bin_Spectra}, e.g. to set the width of bins
#' (binSize). 
#' Prior to calculating \deqn{W_{S1}} or \deqn{W_{S2}}, all intensity values 
#' are divided by the maximum intensity value.
#'
#' @return
#' `normalizeddotproduct` returns a numeric similarity coefficient between
#' 0 and 1
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples 
#' data("spectra", package = "MetCirc")
#' x <- spectra_tissue[[1]]
#' y <- spectra_tissue[[2]]
#' normalizeddotproduct(x, y, m = 0.5, n = 2, binSize = 0.01)
#'
#' @export
normalizeddotproduct <- function(x, y, m = 0.5, n = 2, ...) {
    ## normalize to % intensity
    x@intensity <- x@intensity / max(x@intensity)*100
    y@intensity <- y@intensity / max(y@intensity)*100
    
    binnedSpectra <- MSnbase:::bin_Spectra(x, y, ...)
    inten <- lapply(binnedSpectra, intensity)
    mz <- lapply(binnedSpectra, mz)
    
    ws1 <- inten[[1]] ^ m * mz[[1]] ^ n
    ws2 <- inten[[2]] ^ m * mz[[2]] ^ n
    
    sum( ws1*ws2) ^ 2 / ( sum( ws1^2 ) * sum( ws2^2 ) )
}



#' @name neutralloss
#'
#' @title Calculate similarity based on neutral losses
#'
#' @description Calculate similarity based on neutral losses (NLS)
#'
#' @param x `Spectrum2` object from `MSnbase` containing
#' intensity and m/z values, first MS/MS spectrum
#'
#' @param y `Spectrum2` object from `MSnbase` containing 
#' intensity and m/z values, second MS/MS spectrum
#'
#' @param m `numeric`, exponent to calculate peak intensity-based weights
#'
#' @param n `numeric`, exponent to calculate m/z-based weights
#'
#' @param ... further arguments passed to `MSnbase:::bin_Spectra`
#'
#' @details Similarities of spectra based on neutral losses are calculated 
#' according to the following formula: 
#' \deqn{NLS = \frac{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2}{ \sum(W_{S1, i} ^ 2) * \sum(W_{S2, i} ^ 2) }}{\sum(W_{S1, i} \cdot W_{S2, i}) ^ 2 \sum(W_{S1, i} ^ 2) * \sum(W_{S2, i} ^ 2)},
#' with \eqn{W = [ peak intensity] ^{m} \cdot [ NL ]^n} and 
#' \eqn{NL = | m/z - precursor m/z |}. For further information 
#' see Li et al. (2015): Navigating natural variation in herbivory-induced
#' secondary metabolism in coyote tobacco populations using MS/MS structural 
#' analysis. PNAS, E4147--E4155. `neutralloss} returns a numeric 
#' value ranging between 0 and 1, where 0 
#' indicates no similarity between the two MS/MS features, while 1 indicates 
#' that the MS/MS features are identical. Arguments can be passed to 
#' the function `MSnbase:::bin_Spectra}, e.g. to set the width of bins
#' (binSize). 
#' Prior to calculating \deqn{W_{S1}} or \deqn{W_{S2}}, all intensity values 
#' are divided by the maximum intensity value.
#'
#' @return
#' `neutralloss` returns a numeric similarity coefficient between 0 and 1
#' 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples 
#' data("spectra", package = "MetCirc")
#' x <- spectra_tissue[[1]]
#' y <- spectra_tissue[[2]]
#' neutralloss(x, y, m = 0.5, n = 2, binSize = 0.01)
#'
#' @export
neutralloss <- function(x, y, m = 0.5, n = 2, ...) {

    ## normalize to % intensity
    x@intensity <- x@intensity / max(x@intensity)*100
    y@intensity <- y@intensity / max(y@intensity)*100

    ## calculate loss to precursorMz value
    x@mz <- abs(x@mz - x@precursorMz)
    y@mz <- abs(y@mz - y@precursorMz)

    binnedSpectra <- MSnbase:::bin_Spectra(x, y, ...)
    mz <- lapply(binnedSpectra, mz)
    inten <- lapply(binnedSpectra, intensity)

    ws1 <- inten[[1]] ^ m * mz[[1]] ^ n
    ws2 <- inten[[2]] ^ m * mz[[2]] ^ n

    sum( ws1*ws2) ^ 2 / ( sum( ws1^2 ) * sum( ws2^2 ) )
}

#' @import amap

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
#' @param
#' similarityMatrix `matrix`, `similarityMatrix` contains 
#' pair-wise similarity coefficients which give information about the similarity 
#' between precursors
#'
#' @param
#' spectra `Spectra` object containing spectra that are compared
#' in `similarityMatrix`
#'
#' @param
#' type `character`, one of "retentionTime", "mz" or "clustering"
#'
#' @param
#' group `logical`, if TRUE group separated by "_" will be cleaved
#' from rownames/colnames of similarityMatrix and matched against names of 
#' spectra, if FALSE rownames/colnames of similarityMatrix are taken as are 
#' and matched against names of spectra
#'
#' @details
#' `orderSimilarityMatrix` takes  a similarity matrix,
#' spectra (containing information on m/z and retentionTime and a 
#' `character` vector
#' as arguments. It will then reorder rows and columns of 
#' the similarityMatrix object such, that it orders rows and columns of 
#' similarityMatrix according to m/z, retention time or clustering in 
#' each group. `orderSimilarityMatrix` is employed in the `shinyCircos` 
#' function to create `similarityMatrix` objects which will allow to switch
#' between different types of ordering in between groups (sectors) in the 
#' circos plot. It may be used as well externally, to reproduce plots outside
#' of the reactive environment (see vignette for a workflow).
#'
#' @return `matrix`, `orderSimilarityMatrix` returns a similarity matrix 
#' with ordered rownames according to the `character` vector given to order
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("spectra", package = "MetCirc")
#' similarityMat <- compare_Spectra(spectra_tissue[1:10], 
#'     fun = normalizeddotproduct, binSize = 0.01)
#' ## order according to retention time 
#' orderSimilarityMatrix(similarityMatrix = similarityMat, 
#'     spectra_tissue, type = "retentionTime", group = FALSE)
#'
#' @export
orderSimilarityMatrix <- function(similarityMatrix, spectra, 
        type = c("retentionTime","mz", "clustering"), group = FALSE) {

    if (!(group %in% c(TRUE, FALSE))) stop("group has to be TRUE or FALSE")

    ## set diagonal of similarityMatrix to 1
    diag(similarityMatrix) <- 1
    type <- match.arg(type)
    groupname <- rownames(similarityMatrix)

    ## if a group is preciding the rownames of similarityMatrix
    if (group) {
        groupname <- strsplit(groupname, split = "_")
        groupname <- lapply(groupname, "[", -1)
        groupname <- lapply(groupname, function(x) paste(x, collapse = "_"))
        groupname <- unlist(groupname)
    }

    ## match colnames/rownames of similarityMatrix and names of spectra and
    ## truncate spectra
    spectra <- spectra[names(spectra) %in% groupname, ]

    ## retentionTime
    if (type == "retentionTime") {
        rt <- unlist(lapply(spectra, function(x) x@rt))
        orderNew <- order(rt)
    }

    ## mz
    if (type == "mz") {
        prec_mz <- unlist(lapply(spectra, function(x) x@precursorMz))
        orderNew <- order(prec_mz)
    }

    ## clustering
    if (type == "clustering") {
        hClust <- amap::hcluster(similarityMatrix, method = "spearman")
        orderNew <- hClust$order
    }

    simM <- similarityMatrix[orderNew, orderNew]
    colnames(simM) <- rownames(simM) <- rownames(similarityMatrix)[orderNew]

    return(simM)
}