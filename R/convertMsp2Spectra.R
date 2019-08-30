#' @name convertMsp2Spectra
#' @title Convert MSP data frame into object of class \code{Spectra}
#' @description Convert msp data frame into object of class \code{Spectra}
#' @usage convertMsp2Spectra(msp)
#' @param msp \code{data.frame} that mimicks the .msp file 
#' format, see Details for further information
#' @details msp is a data frame of a .msp file, a typical data file for 
#' MS/MS libraries. The data frame has two columns and contains in the first
#' column the entries "NAME:", 
#' "PRECURSORMZ:" (or "EXACTMASS:"), "RETENTIONTIME:", Num Peaks:"  and information on fragments and 
#' peak areas/intensities and will
#' extract the respective information in the second column.
#' @return \code{convertMsp2Spectra} returns an object of class \code{Spectra}. 
#' @author Thomas Naake, \email{thomasnaake@@googlemail.com}
#' @examples 
#' data("convertMsp2Spectra", package="MetCirc")
#' convertMsp2Spectra(msp=msp2spectra)
#' @export
convertMsp2Spectra <- function(msp) {
    
    nameInd <- grep("NAME:", msp[,1]) ## indices of Name
    
    ## get indices of precursor mz or exactmass
    mzInd_1 <- grep("PRECURSORMZ:", msp[,1])
    mzInd_2 <- grep("EXACTMASS:", msp[,1])
    
    mzInd <- if (length(mzInd_1) >= length(mzInd_2)) {
        mzInd_1 
    } else {
        mzInd_2
    }
    
    numpeaksInd <- grep("Num Peaks:", msp[,1])
    
    ## get indices of "RETENTIONTIME:"
    rtInd <- grep("RETENTIONTIME:", msp[,1])
    
    NAMES <- as.character(msp[nameInd, 2])
    MZ <- as.numeric(as.character(msp[mzInd, 2]))
    NUMPEAKS <- as.numeric(as.character(msp[numpeaksInd, 2]))
    
    ## how many entries are in msp
    numEntries <- length(NAMES)
    
    if (numEntries != length(MZ)) 
        stop("length of precursor mz != length of names")
    
    ## create annotation vectors
    RT <- if(length(rtInd) == numEntries) {
        as.numeric(as.character(msp[rtInd, 2]))
    } else {
        rep(NaN, numEntries)   
    }
       
    spN_l <- vector("list", numEntries)
    
    for (i in 1:numEntries) {
        beg <- numpeaksInd[i] + 1
        end <- numpeaksInd[i] + NUMPEAKS[i]
        fragment <- as.numeric(as.character(msp[beg:end, 1]))
        intensity <- as.numeric(as.character(msp[beg:end, 2]))
        
        ## sort fragments according to increasing fragment values
        fragment <- sort(fragment)
        intensity <- intensity[order(fragment)]
        
        ## delete double entries
        intensity <- intensity[!duplicated(fragment)]
        fragment <- fragment[!duplicated(fragment)]
        
        ## calculate percentages
        ##intensity <- intensity / max(intensity) * 100
        
        spN_l[[i]] <- new("Spectrum2", rt=RT[i], precursorMz=MZ[i], 
                         mz=fragment, intensity=intensity)
    }
    
    spl <- Spectra(spN_l, elementMetadata=DataFrame(names=NAMES)) 
    return(spl)
}
