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
#' similarityMat <- compare_Spectra(spectra_tissue,
#'     fun = normalizeddotproduct, binSize = 0.01)
#' save(similarityMat, file = "similarityMat.RData", compress = "xz")
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
NULL  
