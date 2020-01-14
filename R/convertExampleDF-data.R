#' @name convertExampleDF
#'
#' @title Example data for `MetCirc`: convertExampleDF
#'
#' @description
#' `convertExampleDF` is a `data.frame` which comprises
#' information on a specific metabolite per row stating the average retention
#' time, average m/z, the name of the metabolite, the adduct ion name and the
#' spectrum reference file name. The function `allocatePrecursor2mz` uses
#' `data.frame`s of the kind of `sd01\_outputXCMS` and
#' `sd02\_deconvoluted` to create a `data.frame` of the kind of
#' `convertExampleDF`. Allocation of precursor ions to candidate m/z values is
#' based on minimal distance of m/z and deviance of retention time
#' based on an objective function. See `?allocatePrecursor2mz` for further
#' information.
#'
#' @docType data
#'
#' @return `data.frame`
#'
#' @format `data.frame`
#'
#' @source internal
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
NULL
