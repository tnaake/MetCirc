data("convertMsp2Spectra", package="MetCirc")
spl <- convertMsp2Spectra(msp=msp2spectra)
## START unit test convertMSP2MSP
test_convertMsp2Spectra <- function() {
    checkEquals(length(spl), 22)
    checkEquals(is(spl), c("Spectra", "SimpleList", "List", "Vector", "list_OR_List", "Annotated", "vector_OR_Vector"))
    checkEquals(dim(spl@elementMetadata), c(22, 1))
    checkTrue(is.character(spl@elementMetadata$names))
    checkEquals(spl@elementType, "Spectrum")
    checkEquals(spl@metadata, list())
    checkTrue(is.list(spl@listData))
    checkEquals(is(spl@listData[[1]]), c("Spectrum2", "Spectrum", "Versioned"))
    checkTrue(is.numeric(spl@listData[[1]]@intensity))
    checkTrue(is.numeric(spl@listData[[1]]@mz))
    checkTrue(is.numeric(spl@listData[[1]]@precursorMz))
}
## END unit test convertMSP2MSP

