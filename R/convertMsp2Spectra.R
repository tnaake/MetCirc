#' @name convertMsp2Spectra
#'
#' @title Convert MSP data frame into object of class \code{Spectra}
#'
#' @description Convert msp data frame into object of class [Spectra()]
#'
#' @param
#' msp \code{data.frame} that mimicks the .msp file format, see Details for 
#' further information.
#'
#' @details
#' msp is a data frame of a .msp file, a typical data file for
#' MS/MS libraries. The data frame has two columns and contains in the first
#' column the entries "NAME:",
#' "PRECURSORMZ:" (or "EXACTMASS:"), "RETENTIONTIME:", Num Peaks:"
#' and information on fragments and peak areas/intensities and will
#' extract the respective information in the second column.
#'
#' @return \code{convertMsp2Spectra} returns an object of class `Spectra`
#'
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#'
#' @examples
#' data("convertMsp2Spectra", package = "MetCirc")
#' convertMsp2Spectra(msp = msp2spectra)
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom Spectra Spectra
#' 
#' @export
convertMsp2Spectra <- function(msp) {

    nameInd <- grep("NAME:", msp[, 1]) ## indices of Name

    ## get indices of precursor mz or exactmass
    mzInd_1 <- grep("PRECURSORMZ:", msp[, 1])
    mzInd_2 <- grep("EXACTMASS:", msp[, 1])

    mzInd <- if (length(mzInd_1) >= length(mzInd_2)) {
        mzInd_1
    } else {
        mzInd_2
    }

    numpeaksInd <- grep("Num Peaks:", msp[, 1])

    ## get indices of "RETENTIONTIME:"
    rtInd <- grep("RETENTIONTIME:", msp[, 1])

    NAMES <- as.character(msp[nameInd,  2])
    MZ <- as.numeric(as.character(msp[mzInd, 2]))
    NUMPEAKS <- as.numeric(as.character(msp[numpeaksInd, 2]))

    ## how many entries are in msp
    numEntries <- length(NAMES)

    if (numEntries != length(MZ)) 
        stop("length of precursor mz != length of names")

    ## create annotation vectors
    RT <- if (length(rtInd) == numEntries) {
        as.numeric(as.character(msp[rtInd, 2]))
    } else {
        rep(NaN, numEntries)
    }

    sps_l <- vector("list", numEntries)

    for (i in seq_len(numEntries)) {
        
        beg <- numpeaksInd[i] + 1
        end <- numpeaksInd[i] + NUMPEAKS[i]
        mzs <- as.character(msp[beg:end, 1]) |>
            as.numeric()
        intensity <- as.character(msp[beg:end, 2]) |>
            as.numeric()

        ## sort fragments according to increasing fragment values
        sorted_intensity <- intensity[order(mzs)]
        sorted_mzs <- sort(mzs)

        ## delete double entries
        sorted_intensity <- sorted_intensity[!duplicated(sorted_mzs)]
        sorted_mzs <- sorted_mzs[!duplicated(sorted_mzs)]

        ## calculate percentages
        ##intensity <- intensity / max(intensity) * 100
        spd <- S4Vectors::DataFrame(
            msLevel = 2L,
            rtime = RT[i],
            precursorMz = MZ[i],
            name = NAMES[i]
        )
        spd$intensity <- list(sorted_intensity)
        spd$mz <- list(sorted_mzs)
        
        ## create Spectra object and write to list entry
        sps_l[[i]] <- Spectra::Spectra(spd)
    }
    
    ## combine the different entries and return 
    sps <- Reduce(c, sps_l)
    sps
}

