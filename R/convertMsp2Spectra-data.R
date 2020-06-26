#' @name msp2spectra
#'
#' @title Example data for `MetCirc`: `msp2spectra`
#'
#' @description
#' `convertMsp2Spectra` contains the object `msp2spectra`
#' that is a data frame in .MSP format, a typical format for MS/MS library
#' building. Each entry consists of the metabolite name (NAME), the precursor
#' mz (PRECURSORMZ), the retention time (RETENTIONTIME), number of peaks
#' (Num Peaks), together with fragments and their
#' intensity values. In the example used in the function
#' `convertMsp2Spectra` the `data.frame` `msp2spectra`
#' is used to construct an object of class `MSpectra`.
#'
#' @docType data
#'
#' @return `data.frame`
#'
#' @format `data.frame`
#'
#' @source http://prime.psc.riken.jp/Metabolomics_Software/MS-DIAL/, truncated
#' .MSP file of GNPS MS/MS Negative (contains 22 entries):
#' http://prime.psc.riken.jp/Metabolomics_Software/MS-DIAL/MSMS-GNPS-Curated-Neg.msp
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
NULL
