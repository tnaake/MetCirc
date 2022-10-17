#' @name similarityMat
#'
#' @title Example data for `MetCirc`: `similarityMat`
#'
#' @description
#' `similarityMat` is a `matrix` containing the pair-wise
#' similarity scores derived from the `idMSMStissueproject` data set.
#' See the vignette for a workflow to reproduce the object `similarityMat`.
#'
#' @docType data
#'
#' @return `matrix`
#'
#' @format `matrix`
#'
#' @source
#' data("spectra", package = "MetCirc")
#' similarityMat <- Spectra::compareSpectra(sps_tissue,
#'     fun = ndotproduct, ppm = 10)
#' rownames(similarityMat) <- colnames(similarityMat) <- sps_tissue$name
#' save(similarityMat, file = "similarityMat.RData", compress = "xz")
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
NULL  
