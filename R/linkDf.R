#' @name spectraCond
#' @title Get MS/MS spectra that are present in condition 
#' @description `spectraCond` returns the names of `spectra` that are 
#' present in condition, corresponding to the slot 
#' `elementMetadata@listData`.
#' @usage spectraCond(spectra, condition)
#' @param spectra `Spectra` object of `MSnbase` package
#' @param condition `character`, vector with conditions found as columns 
#' in the elementMetadata slot 
#' @details Helper function in `createLink0df` and `shinyCircos`.
#' @return `list`, named `list` with `character` vector as entries that contains 
#' the names of the MS/MS entries in `spectra` that are present in the 
#' `conditon` (tissues, stress conditions, time points, etc.)
#' @author Thomas Naake \email{thomasnaake@@googlemail.com}
#' @examples
#' data("spectra", package="MetCirc")
#' MetCirc:::spectraCond(spectra_tissue, condition=c("SPL", "LIM", "ANT", "STY"))
spectraCond <- function(spectra, condition) {
    inds <- lapply(condition, function(x) which(spectra@elementMetadata@listData[[x]] == 1))
    names(inds) <- condition
    inds <- lapply(inds, function(x) names(spectra)[x])
    
    return(inds)
}


#' @name createLink0df
#' @title Create a link matrix 
#' @description Create a link matrix which links every feature in similarity
#' matrix with another. 
#' @usage createLink0df(similarityMatrix, spectra, condition)
#' @param similarityMatrix `matrix`, a similarity matrix that contains the 
#' NDP similarity measure between all precursors in the data set
#' @param spectra `Spectra` object
#' @param condition `character`, which conditions should be included?
#' @details createLink0df creates a `matrix` from a similarity 
#' matrix which includes all connections between features in the 
#' similarity matrix, but 
#' exclude links which have a similarity of exactly 0.
#' @return createLink0df returns a `matrix` that gives per each row 
#' information on linked features
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("spectra", package="MetCirc")
#' data("similarityMat", package="MetCirc")
#' link0df <- createLink0df(similarityMatrix=similarityMat, 
#'     spectra_tissue, condition=c("SPL", "LIM", "ANT", "STY"))
#' @export
createLink0df <- function(similarityMatrix, spectra, condition) { 
    
    if (!all(colnames(similarityMatrix) == rownames(similarityMatrix))) {
        stop("colnames(similarityMatrix) != rownames(similarityMatrix)")
    } 
    
    ## set diag to matrix (all similarities similar to 0 will be removed in 
    ## a later step)
    diag(similarityMatrix) <- 0 
    
    ## truncate spectra that it has the same spectra as in similarityMatrix
    spectra <- spectra[names(spectra) %in% colnames(similarityMatrix),]
    
    inds <- MetCirc:::spectraCond(spectra, condition)
    inds_cond <- lapply(seq_along(inds), 
        function(x) {
            if (length(inds[[x]]) > 0) {
                paste(condition[x], inds[[x]], sep="_")
            } else character()
    })
    inds_uniq <- unique(unlist(inds_cond))
    
    inds_uniq_combn <- combn(inds_uniq, 2)
    
    ## get similarity values for all combinations
    sim <- lapply(1:ncol(inds_uniq_combn), function(x) {
        row_sim <- strsplit(inds_uniq_combn[1, x], split="_")[[1]][2]
        col_sim <- strsplit(inds_uniq_combn[2, x], split="_")[[1]][2]
        similarityMatrix[row_sim, col_sim]
    })
    sim <- unlist(sim)
    
    ## cbind spectrum1, spectrum2, similarity and remove rows where sim == 0
    mat0 <- data.frame(group1=unlist(lapply(strsplit(inds_uniq_combn[1,], split="_"), "[", 1)),
        spectrum1=inds_uniq_combn[1, ], 
        group2=unlist(lapply(strsplit(inds_uniq_combn[2,], split="_"), "[", 1)), 
        spectrum2=inds_uniq_combn[2, ], 
        similarity=as.numeric(sim))
    mat0 <- mat0[!mat0[, "similarity"] == 0,]
    return(mat0)
}

#' @name thresholdLinkDf
#' @title Threshold a data frame containing information on links
#' @description Threshold a link data frame
#' @usage thresholdLinkDf(link0df, lower=0.75, upper=1)
#' @param link0df `data.frame`, a link data frame that gives per each row 
#' information on linked features
#' @param lower `numeric`, threshold value for similarity values, below 
#' this value linked features will not be returned
#' @param upper `numeric`, threshold value for similarity values, above 
#' this value linked features will not be returned
#' @details `lower` and `upper` are numerical values 
#' and truncates mass spectra based on their similarity values. 
#' @return \code{thresholdLinkDf} returns a data frame that gives per each row 
#' information on linked features which are linked within certain thresholds
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("spectra", package="MetCirc")
#' data("similarityMat", package="MetCirc")
#' link0df <- createLink0df(similarityMatrix=similarityMat, 
#'     spectra_tissue, condition=c("SPL", "LIM", "ANT", "STY"))
#' thresholdLinkDf(link0df=link0df, lower=0.5, upper=1)
#' @export
thresholdLinkDf <- function(link0df, lower=0.75, upper=1) {
    
    if (!all(colnames(link0df) == c("group1", "spectrum1",  "group2", "spectrum2",  "similarity")))
        stop("linkDF does not have right colnames")
    
    sim <- link0df[, "similarity"]
    if (lower > upper) stop("lower greater than upper")
    if (upper > 1) stop("upper greater than 1")
    if (lower > max(sim)) 
        warning("lower greater than max similarity value in link0df")
    
    ## which rows have a coefficient >= threshold?
    indThr <- which(sim >= lower & sim <= upper)
    
    ## cut linkDf
    if (length(indThr) <= 1) {
        thrDf <- matrix(NA, ncol=ncol(link0df), nrow=length(indThr))
        thrDf[1:nrow(thrDf),1:ncol(thrDf)] <- link0df[indThr,]
    } else {
        thrDf <- link0df[indThr,]
    }
    
    colnames(thrDf) <- colnames(link0df)  
    
    return(thrDf)
}

#' @name createLinkDf
#' @title Create a data frame which contains features to link (indices)
#' @description Create a data frame which contains features to link (indices)
#' @usage createLinkDf(similarityMatrix, spectra, condition, lower, upper) 
#' @param similarityMatrix `matrix`, a similarity matrix that contains the 
#' similarity measure between all precursors in the data set
#' @param spectra Spectra object containing spectra of similarityMatrix
#' @param condition `character`, vector containing the 
#' conditions/samples for which a linkDf is created
#' @param lower `numeric`, threshold value for similarity values, 
#' below this value linked features will not be included
#' @param upper `numeric`, threshold value for similarity values, 
#' above this value linked features will not be included
#' @details `lower` and `upper` are numerical values 
#' and truncate similar spectra. The function createLinkDf is a wrapper 
#' for the functions `createLink0df` and `thresholdLinkDf`.
#' @return `createLinkDf` returns a data frame that gives per each row 
#' information on linked features
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("spectra", package="MetCirc")
#' data("similarityMat", package="MetCirc")
#' link0df <- createLink0df(similarityMatrix=similarityMat, 
#'     spectra_tissue, condition=c("SPL", "LIM", "ANT", "STY"))
#' createLinkDf(similarityMatrix=similarityMat, spectra=spectra_tissue,
#'     condition=c("SPL", "LIM", "ANT", "STY"), lower=0.5, upper=1)
#' @export
createLinkDf <- function(similarityMatrix, spectra, condition, lower, upper) {
    ## first create a link0Matrix
    link0df <- createLink0df(similarityMatrix=similarityMatrix, 
        spectra=spectra, condition=condition)
    ## than threshold link0Matrix
    thrLinkDf <- thresholdLinkDf(link0df=link0df, lower=lower, upper=upper)
     return(thrLinkDf)
}

#' @name cutLinkDf
#' @title Create a cut data frame with information on links
#' @description Create a cut link data frame
#' @usage cutLinkDf(linkDf, type=c("all", "inter", "intra"))
#' @param linkDf `data.frame`, that gives per each row 
#' information on linked features
#' @param type `character`, one of "all", "inter" or "intra"
#' @details This function is used to truncate features from linkDf. If 
#' type="all", linkDf will not be changed; if type="inter" the returned
#' linkDf will only contain entries of links which are between groups and 
#' not inside groups; contrary to that, if type="intra" the returned linkDf 
#' will only contain entries of links which are inside groups and not between 
#' groups.
#' @return cutLinkDf returns a data.frame that gives per each row 
#' information on linked features
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("spectra", package="MetCirc")
#' data("similarityMat", package="MetCirc")
#' linkDf <- createLinkDf(similarityMatrix=similarityMat, spectra=spectra_tissue, 
#'     condition=c("SPL", "LIM", "ANT", "STY"), lower=0.75, upper=1)
#' cutLinkDf(linkDf=linkDf, type="all")
#' @export
cutLinkDf <- function(linkDf, type=c("all", "inter", "intra")) {
    
    ## match arguments of types
    type <- match.arg(type)
    
    if (type == "all") 
        linkDf <- linkDf
    if (type == "inter")
        linkDf <- linkDf[which(linkDf[,"group1"] != linkDf[,"group2"]), ]
    if (type == "intra") 
        linkDf <- linkDf[which(linkDf[,"group1"] == linkDf[,"group2"]), ]
    return(linkDf)
}

