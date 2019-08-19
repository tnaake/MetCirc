#' @name msp2msp
#' @title Example data for \code{MetCirc}: \code{msp2spectra}
#' @description \code{convertMsp2Spectra} contains the object \code{msp2spectra} 
#' that is a data frame in .MSP format, a typical format for MS/MS library 
#' building. Each entry consists of the metabolite name (NAME), the precursor 
#' mz (PRECURSORMZ), the retention 
#' time (RETENTIONTIME), number of peaks (Num Peaks), together with fragments and their 
#' intensity values. In the example used in the function
#' \code{convertMsp2Spectra} the \code{data.frame} \code{msp2spectra} 
#' is used to construct an object of class \code{Spectra}. 
#' @docType data
#' @usage msp2spectra
#' @return \code{data.frame}
#' @format \code{data.frame}
#' @source http://prime.psc.riken.jp/Metabolomics_Software/MS-DIAL/, truncated
#' .MSP file of GNPS MS/MS Negative (contains 22 entries): 
#' http://prime.psc.riken.jp/Metabolomics_Software/MS-DIAL/MSMS-GNPS-Curated-Neg.msp
#' @author Thomas Naake, \email{thomasnaake@googlemail.com}
NULL  
