#' @name spectraCondition
#'
#' @title Get MS/MS spectra that are present in condition
#'
#' @description
#' \code{spectraCondition} returns the names of \code{Spectra} that are
#' present in condition, corresponding to the slot
#' \code{metadata}.
#'
#' @param sps \code{Spectra} object of \code{Spectra} package
#' @param condition \code{character}, vector with conditions found as columns
#' in the metadata slot
#'
#' @details
#' Helper function in \code{createLink0df} and \code{shinyCircos}.
#'
#' @return
#' \code{list}, named \code{list} with \code{character} vector as entries that
#' contains the names of the MS/MS entries in \code{spectra} that are present 
#' in the \code{conditon} (tissues, stress conditions, time points, etc.)
#'
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("spectra", package = "MetCirc")
#' MetCirc:::spectraCondition(sps = sps_tissue,
#'     condition = c("SPL", "LIM", "ANT", "STY"))
spectraCondition <- function(sps, condition) {

    ## get the indices 
    inds <- lapply(condition, function(condition_i) {
        which(sps@metadata[[condition_i]] == 1)
    })
    
    ## return the names of the features that are present in condition
    inds <- lapply(inds, function(inds_i) {
        inds_i <- sps$name[inds_i]
        inds_i[!is.na(inds_i)]
    })
    names(inds) <- condition
    inds
}

#' @name createLink0df
#'
#' @title Create a link matrix
#'
#' @description
#' Create a link matrix which links every feature in similarity
#' matrix with another.
#'
#' @param similarityMatrix \code{matrix}, a similarity matrix that contains the
#' NDP similarity measure between all precursors in the data set
#' @param sps \code{Spectra} object
#' @param condition \code{character}, which conditions should be included?
#'
#' @details
#' \code{createLink0df} creates a \code{matrix} from a similarity matrix which 
#' includes all connections between features in the similarity matrix, but
#' exclude links which have a similarity of exactly 0.
#'
#' @return \code{createLink0df} returns a `matrix` that gives per each row
#' information on linked features
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("spectra", package = "MetCirc")
#' data("similarityMat", package = "MetCirc")
#' link0df <- createLink0df(similarityMatrix = similarityMat,
#'     sps = sps_tissue, condition = c("SPL", "LIM", "ANT", "STY"))
#'
#' @importFrom utils combn
#' 
#' @export
createLink0df <- function(similarityMatrix, sps, condition) {

    if (!all(colnames(similarityMatrix) == rownames(similarityMatrix))) {
        stop("colnames(similarityMatrix) != rownames(similarityMatrix)")
    }

    ## set diag to matrix (all similarities similar to 0 will be removed in
    ## a later step)
    diag(similarityMatrix) <- 0

    ## truncate sps that it has the same sps as in similarityMatrix
    sps <- sps[sps$name %in% colnames(similarityMatrix), ]

    inds <- spectraCondition(sps, condition)
    inds_cond <- lapply(seq_along(inds),
        function(x) {
            if (length(inds[[x]]) > 0) {
                paste(condition[x], inds[[x]], sep = "_")
            } else character()
    })

    inds_uniq <- unique(unlist(inds_cond))

    inds_uniq_combn <- utils::combn(inds_uniq, 2)

    ## get similarity values for all combinations
    sim <- lapply(seq_len(ncol(inds_uniq_combn)), function(x) {
        row_sim <- strsplit(inds_uniq_combn[1, x], split = "_")[[1]][2]
        col_sim <- strsplit(inds_uniq_combn[2, x], split = "_")[[1]][2]
        similarityMatrix[row_sim, col_sim]
    })

    sim <- unlist(sim)

    ## cbind spectrum1, spectrum2, similarity and remove rows where sim == 0
    mat0 <- data.frame(group1 = unlist(lapply(
            strsplit(inds_uniq_combn[1, ], split = "_"), "[", 1)),
        spectrum1 = inds_uniq_combn[1, ],
        group2 = unlist(lapply(
            strsplit(inds_uniq_combn[2, ], split = "_"), "[", 1)),
        spectrum2 = inds_uniq_combn[2, ],
        similarity = as.numeric(sim))

    mat0[!mat0[, "similarity"] == 0, ]
}

#' @name thresholdLinkDf
#'
#' @title Threshold a data frame containing information on links
#'
#' @description 
#' Threshold a link data frame based on lower and upper similarity values. 
#' The function will return the links that lie within the defined bounds.
#'
#' @param link0df \code{data.frame}, a link data frame that gives per each row
#' information on linked features
#'
#' @param lower \code{numeric}, threshold value for similarity values, below
#' this value linked features will not be returned
#'
#' @param
#' upper \code{numeric}, threshold value for similarity values, above
#' this value linked features will not be returned
#'
#' @details
#' \code{lower} and \code{upper} are numerical values
#' and truncate mass spectra based on their similarity values.
#'
#' @return
#' \code{thresholdLinkDf} returns a \code{data.frame} that gives per each row
#' information on linked features which are linked within certain thresholds.
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("spectra", package = "MetCirc")
#' data("similarityMat", package = "MetCirc")
#' link0df <- createLink0df(similarityMatrix = similarityMat,
#'     sps_tissue, condition = c("SPL", "LIM", "ANT", "STY"))
#' thresholdLinkDf(link0df = link0df, lower = 0.5, upper = 1)
#'
#' @export
thresholdLinkDf <- function(link0df, lower = 0.75, upper = 1) {

    if (!all(colnames(link0df) == c("group1", "spectrum1",  "group2", "spectrum2",  "similarity")))
        stop("linkDF does not have right colnames")

    sim <- link0df[, "similarity"]
    if (lower > upper) stop("lower greater than upper")
    if (upper > 1) stop("upper greater than 1")
    if (lower > max(sim))
        warning("lower greater than max similarity value in link0df")

    ## which rows have a coefficient >= threshold?
    indThr <- which(sim >= lower & sim <= upper)

    ## cut linkDf and return
    link0df[indThr, ]

}

#' @name createLinkDf
#'
#' @title Create a data frame which contains features to link (indices)
#'
#' @description
#' Create a data frame which contains features to link (indices).
#'
#' @param similarityMatrix \code{matrix}, a similarity matrix that contains the
#' similarity measure between all precursors in the data set
#'
#' @param sps \code{Spectra} object containing spectral data corresponding to 
#' features in \code{similarityMatrix}
#'
#' @param condition \code{character} \code{vector} containing the
#' conditions/samples for which a \code{linkDf} is created
#'
#' @param lower \code{numeric(1)}, threshold value for similarity values,
#' linked features below this value will not be included
#'
#' @param upper \code{numeric(1)}, threshold value for similarity values,
#' linked features above this value will not be included
#'
#' @details 
#' \code{lower} and \code{upper} are numerical values and truncate 
#' spectra based on their similarity. The function \code{createLinkDf} is 
#' a wrapper for the functions \code{createLink0df} and \code{thresholdLinkDf}.
#'
#' @return \code{createLinkDf} returns a \code{data.frame} that gives per each 
#' row information on linked features
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("spectra", package = "MetCirc")
#' data("similarityMat", package = "MetCirc")
#' link0df <- createLink0df(similarityMatrix = similarityMat,
#'     sps = sps_tissue, condition = c("SPL", "LIM", "ANT", "STY"))
#' createLinkDf(similarityMatrix = similarityMat, sps = sps_tissue,
#'     condition = c("SPL", "LIM", "ANT", "STY"), lower = 0.5, upper = 1)
#'
#' @export
createLinkDf <- function(similarityMatrix, sps, condition, lower, upper) {

    ## first create a link0Matrix
    link0df <- createLink0df(similarityMatrix = similarityMatrix,
        sps = sps, condition = condition)

    ## than threshold link0Matrix and return
    thresholdLinkDf(link0df = link0df, lower = lower, upper = upper)

}

#' @name cutLinkDf
#'
#' @title Create a cut data frame with information on links
#'
#' @description Create a cut link data frame
#'
#' @param linkDf \code{data.frame}, that gives per each row
#' information on linked features
#'
#' @param type \code{character}, one of "all", "inter" or "intra"
#'
#' @details
#' This function is used to truncate features from \code{linkDf}. If
#' \code{type = "all"}, \code{linkDf} will not be changed; if 
#' \code{type = "inter"} the returned \code{linkDf} will only contain entries 
#' of links which are between groups and  not inside groups; contrary to that,
#' if \code{type = "intra"} the returned \code{linkDf} will only contain entries 
#' of links which are inside groups and not between groups.
#'
#' @return
#' \code{cutLinkDf} returns a \code{data.frame} that gives per each row
#' information on linked features
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("spectra", package = "MetCirc")
#' data("similarityMat", package = "MetCirc")
#' linkDf <- createLinkDf(similarityMatrix = similarityMat,
#'     sps = sps_tissue, condition = c("SPL", "LIM", "ANT", "STY"),
#'     lower = 0.75, upper = 1)
#' cutLinkDf(linkDf = linkDf, type = "all")
#'
#' @export
cutLinkDf <- function(linkDf, type = c("all", "inter", "intra")) {

    ## match arguments of types
    type <- match.arg(type)

    if (type == "all")
        linkDf <- linkDf

    if (type == "inter")
        linkDf <- linkDf[which(linkDf[, "group1"] != linkDf[, "group2"]), ]

    if (type == "intra")
        linkDf <- linkDf[which(linkDf[, "group1"] == linkDf[, "group2"]), ]

    linkDf
}
