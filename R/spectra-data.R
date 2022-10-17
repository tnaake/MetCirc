#' @name sps_tissue
#'
#' @title Example data for `MetCirc`: `sps_tissue`
#'
#' @description `sps_tissue` is a `Spectra` object
#' derived from the `idMSMStissueproject` data set.
#' See the vignette for a workflow to reproduce the object `spectra`.
#'
#' @docType data
#'
#' @return `matrix`
#'
#' @format `matrix`
#'
#' @source
#' data("idMSMStissueproject", package = "MetCirc")
#' ## get all MS/MS spectra
#' tissue <- tissue[tissue[, "id"] %in% compartmentTissue[, "mz_rt_pcgroup"], ]
#' id_uniq <- unique(tissue[, "id"])
#'
#' ## obtain precursor m/z from id_uniq
#' prec_mz <- lapply(strsplit(as.character(id_uniq), split = "_"), "[", 1) |>
#'     unlist() |>
#'     as.numeric()
#'
#' ## obtain m/z from fragments per precursor m/z
#' mz_l <- lapply(id_uniq, function(id_i) tissue[tissue[, "id"] == id_i, "mz"])
#' 
#' ## obtain corresponding intensity values
#' int_l <- lapply(id_uniq, function(id_i) tissue[tissue[, "id"] == id_i, "intensity"])
#'     
#' ## order mz and intensity values
#' int_l <- lapply(seq_along(int_l), function(i) int_l[[i]][order(mz_l[[i]])])
#' mz_l <- lapply(seq_along(mz_l), function(i) sort(mz_l[[i]]))
#' 
#' ## obtain retention time by averaging all retention time values
#' rt <- lapply(id_uniq, function(id_i) tissue[tissue[, "id"] == id_i, "rt"]) |>
#'     lapply(mean) |>
#'     unlist()
#'
#' ## create list of Spectra objects and concatenate
#' sps_l <- lapply(seq_len(length(mz_l)), function(i) {
#'     spd <- S4Vectors::DataFrame(
#'         name = as.character(i),
#'         rtime = rt[i], 
#'         msLevel = 2L,
#'         precursorMz = prec_mz[i])
#'     spd$mz <- list(mz_l[[i]])
#'     spd$intensity <- list(int_l[[i]])
#'     Spectra::Spectra(spd)})
#' sps_tissue <- Reduce(c, sps_l)
#'
#' ## combine list of spectrum2 objects to MSpectra object, 
#' ## use SPL, LIM, ANT, STY for further analysis
#' sps_tissue@metadata <- data.frame(
#'     compartmentTissue[, c("SPL", "LIM", "ANT", "STY")])
#'
#' save(sps_tissue, file = "spectra.RData", compress = "xz")
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
NULL  
