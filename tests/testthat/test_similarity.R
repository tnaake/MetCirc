data("convertMsp2Spectra", package = "MetCirc")
spl <- convertMsp2Spectra(msp2spectra)

simMat <- compare_Spectra(spl[1:2], fun = "dotproduct")

## START unit test compare_Spectra
test_that("compare_Spectra", {
    expect_equal(simMat, 
        matrix(c(1, 1.004575e-05, 1.004575e-05, 1),
            ncol = 2, dimnames = list(1:2, 1:2)), tolerance = 1e-6)
    expect_error(compare_Spectra(spl[1], "dotproduct"), "n < m")
    expect_error(compare_Spectra(spl[1:2], "dotproduct2"), 
        "'arg' should be one of ")
})
## END unit test compare_Spectra

## START unit test normalizeddotproduct
ndp <- normalizeddotproduct(spl[[1]], spl[[2]])
test_that("normalizeddotproduct", {
    expect_equal(ndp, 1.591002e-07, tolerance = 1e-6)
    expect_error(normalizeddotproduct(spl[[1]], 1), 
        "trying to get slot \"intensity\"")
    expect_error(normalizeddotproduct(1, spl[[1]]),
        "trying to get slot \"intensity\"")
    expect_error(normalizeddotproduct(spl[[1]], spl[[2]], "a"),
        "non-numeric argument to binary operator")
    expect_error(normalizeddotproduct(spl[[1]], spl[[2]], 1, "a"),
        "non-numeric argument to binary operator")
})
## END unit test normalizeddotproduct

## START unit test neutralloss
nl <- neutralloss(spl[[1]], spl[[2]])
test_that("neutralloss", {
    expect_equal(nl, 5.899399e-05, tolerance = 1e-6)
    expect_error(neutralloss(spl[[1]], 1),
        "trying to get slot \"intensity\"")
    expect_error(neutralloss(1, spl[[1]]),
        "trying to get slot \"intensity\"")
    expect_error(neutralloss(spl[[1]], spl[[2]], "a"),
        "non-numeric argument to binary operator")
    expect_error(neutralloss(spl[[1]], spl[[2]], 1, "a"),
        "non-numeric argument to binary operator")
})
## END unit test neutralloss


## START unit test orderSimilarityMatrix
simMat_o_mz <- orderSimilarityMatrix(simMat, spl, type = "mz", group = FALSE)
simMat_o_rt <- orderSimilarityMatrix(simMat, spl, type = "retentionTime",
    group = FALSE)
simMat_o_cl <- orderSimilarityMatrix(simMat, spl, type = "clustering",
    group = FALSE)

## create a matrix with groups
simMat_gr <- simMat
rownames(simMat_gr) <- colnames(simMat_gr) <- paste("A", rownames(simMat),
    sep = "_")

test_that("orderSimilarityMatrix", {
    expect_error(orderSimilarityMatrix(simMat, spl, type = "foo"),
        "'arg' should be one of ")
    expect_error(orderSimilarityMatrix(simMat, spl, type = "mz", group = "a"),
        "group has to be TRUE or FALSE")
    #expect_error(orderSimilarityMatrix(type = "mz"), 
    #    "is missing, with no default")
    #expect_error(orderSimilarityMatrix(type = "retentionTime"), 
    #    "is missing, with no default")
    expect_error(orderSimilarityMatrix(type = "clustering"),
        "is missing, with no default")
    expect_equal(colnames(simMat), colnames(simMat_o_mz))
    expect_equal(rownames(simMat), rownames(simMat_o_mz))
    expect_equal(colnames(simMat_o_mz), rownames(simMat_o_mz))
    expect_equal(dim(simMat), dim(simMat_o_mz))
    expect_true(is.matrix(simMat_o_mz))
    expect_true(is.numeric(simMat_o_mz))
    expect_equal(colnames(simMat_o_rt), rownames(simMat_o_rt))
    expect_equal(dim(simMat), dim(simMat_o_rt))
    expect_true(is.matrix(simMat_o_rt))
    expect_true(is.numeric(simMat_o_rt))
    expect_equal(colnames(simMat_o_cl), rownames(simMat_o_cl))
    expect_equal(dim(simMat), dim(simMat_o_cl))
    expect_true(is.matrix(simMat_o_cl))
    expect_true(is.numeric(simMat_o_cl))
    expect_true(is.matrix(
        orderSimilarityMatrix(simMat_gr, spl, type = "mz", group = TRUE)))
    expect_true(is.numeric(
        orderSimilarityMatrix(simMat_gr, spl, type = "mz", group = TRUE)))
    expect_equal(rownames(orderSimilarityMatrix(simMat_gr, spl, type = "mz",
            group = TRUE)), c("A_1", "A_2"))
    expect_equal(colnames(orderSimilarityMatrix(simMat_gr, spl, type = "mz",
            group = TRUE)), c("A_1", "A_2"))
})
## END unit test orderSimilarityMatrix