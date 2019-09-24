#' @name spectra_tissue
#' @title Example data for \code{MetCirc}: \code{spectra_tissue}
#' @description \code{spectra_tissue} is a \code{Spectra} object containing
#' \code{Spectrum2} objects derived from the \code{idMSMStissueproject} data set.
#' See the vignette for a workflow to reproduce the object \code{spectra}.
#' @docType data
#' @usage similarityMat
#' @return \code{matrix}
#' @format \code{matrix}
#' @source
#' data("idMSMStissueproject", package = "MetCirc")
#' ## get all MS/MS spectra
#' tissue <- tissue[tissue[, "id"] %in% compartmentTissue[, "mz_rt_pcgroup"], ]
#' id_uniq <- unique(tissue[, "id"])
#'
#' ## obtain precursor m/z from id_uniq
#' prec_mz_l <- lapply(strsplit(as.character(id_uniq), split = "_"), "[", 1)
#' prec_mz_l <- lapply(prec_mz_l, as.numeric)
#' 
#' ## obtain m/z from fragments per precursor m/z
#' mz_l <- lapply(id_uniq, function(x) tissue[tissue[, "id"] == x, "mz"])
#' ## obtain corresponding intensity values
#' int_l <- lapply(id_uniq, function(x) tissue[tissue[, "id"] == x, "intensity"])
#' ## obtain retention time by averaging all retention time values
#' rt_l <- lapply(id_uniq, function(x) tissue[tissue[, "id"] == x, "rt"])
#' rt_l <- lapply(rt_l, mean)
#' 
#' ## create list of spectrum2 objects
#' spectrum2_tissue <- lapply(1:length(mz_l), function(x) {
#'         new("Spectrum2", rt = rt_l[[x]], precursorMz = prec_mz_l[[x]], 
#'         mz = mz_l[[x]], intensity = int_l[[x]])})
#' 
#' ## combine list of spectrum2 objects to Spectra object, 
#' ## use SPL, LIM, ANT, STY for further analysis
#' spectra_tissue <- Spectra(spectrum2_tissue, 
#'     elementMetadata = DataFrame(compartmentTissue[, c("SPL", "LIM", "ANT", "STY")])) 
#' 
#' save(spectra_tissue, file = "spectra.RData", compress = "xz")
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com} 
NULL  
