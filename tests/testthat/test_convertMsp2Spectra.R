data("convertMsp2Spectra", package = "MetCirc")

spl <- convertMsp2Spectra(msp = msp2spectra)
## START unit test convertMSP2MSP
test_that("convertMsp2Spectra", {
    expect_equal(length(spl), 22)
    expect_equal(is(spl), c("Spectra", "SimpleList", "List", "Vector", 
        "list_OR_List", "Annotated", "vector_OR_Vector"))
    expect_equal(dim(spl@elementMetadata), c(22, 1))
    expect_true(is.character(spl@elementMetadata$names))
    expect_equal(spl@elementType, "Spectrum")
    expect_equal(spl@metadata, list())
    expect_true(is.list(spl@listData))
    expect_equal(is(spl@listData[[1]]), c("Spectrum2", "Spectrum", "Versioned"))
    expect_true(is.numeric(spl@listData[[1]]@intensity))
    expect_true(is.numeric(spl@listData[[1]]@mz))
    expect_true(is.numeric(spl@listData[[1]]@precursorMz))
})
## END unit test convertMSP2MSP