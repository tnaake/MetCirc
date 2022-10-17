data("convertMsp2Spectra", package = "MetCirc")
spl <- convertMsp2Spectra(msp2spectra)

simMat <- Spectra::compareSpectra(spl, FUN = ndotproduct, m = 0.5, n = 2)
rownames(simMat) <- colnames(simMat) <- spl$name

## START unit test compareSpectra
test_that("compareSpectra", {
    expect_equal(sum(simMat), 38.22854, tolerance = 1e-06)
    expect_equal(simMat["Isoquercitrin M-H     ", "Rutin M-H     "], 0.9887026,
        tolerance = 1e-06)
})
## END unit test compareSpectra

## START unit test neutralloss
nl <- neutralloss(
    matrix(c(unlist(spl[1]$mz), unlist(spl[1]$intensity)), ncol = 2), 
    matrix(c(unlist(spl[1]$mz), unlist(spl[1]$intensity)), ncol = 2))
test_that("neutralloss", {
    expect_equal(nl, 1, tolerance = 1e-6)
    simMat_nl <- Spectra::compareSpectra(spl, FUN = neutralloss, m = 0.5, n = 2)
    rownames(simMat_nl) <- colnames(simMat_nl) <- spl$name
    expect_equal(sum(simMat_nl), 25.9532, tolerance = 1e-04)
    expect_equal(simMat_nl["Isoquercitrin M-H     ", "Rutin M-H     "], 0.182,
                 tolerance = 1e-03)
})
## END unit test neutralloss


## START unit test orderSimilarityMatrix
simMat_o_mz <- orderSimilarityMatrix(simMat, sps = spl, type = "mz", group = FALSE)
simMat_o_rt <- orderSimilarityMatrix(simMat, sps = spl, type = "retentionTime",
    group = FALSE)
simMat_o_cl <- orderSimilarityMatrix(simMat, sps = spl, type = "clustering",
    group = FALSE)

## create a matrix with groups
simMat_gr <- simMat
rownames(simMat_gr) <- colnames(simMat_gr) <- paste("A", rownames(simMat),
    sep = "_")

test_that("orderSimilarityMatrix", {
    expect_error(orderSimilarityMatrix(simMat, sps = spl, type = "foo"),
        "'arg' should be one of ")
    expect_error(orderSimilarityMatrix(simMat, sps = spl, type = "mz", group = "a"),
        "group has to be TRUE or FALSE")
    expect_error(orderSimilarityMatrix(type = "mz"), 
        "is missing, with no default")
    expect_error(orderSimilarityMatrix(type = "retentionTime"), 
        "is missing, with no default")
    expect_error(orderSimilarityMatrix(type = "clustering"),
        "is missing, with no default")
    expect_equal(colnames(simMat_o_mz), spl$name[order(spl$precursorMz)])
    expect_equal(rownames(simMat_o_mz), spl$name[order(spl$precursorMz)])
    expect_equal(colnames(simMat_o_mz), rownames(simMat_o_mz))
    expect_equal(dim(simMat), dim(simMat_o_mz))
    expect_true(is.matrix(simMat_o_mz))
    expect_true(is.numeric(simMat_o_mz))
    expect_equal(colnames(simMat_o_rt), spl$name[order(spl$rtime)])
    expect_equal(rownames(simMat_o_rt), spl$name[order(spl$rtime)])
    expect_equal(colnames(simMat_o_rt), rownames(simMat_o_rt))
    expect_equal(dim(simMat), dim(simMat_o_rt))
    expect_true(is.matrix(simMat_o_rt))
    expect_true(is.numeric(simMat_o_rt))
    expect_equal(colnames(simMat_o_cl), rownames(simMat_o_cl))
    expect_equal(dim(simMat), dim(simMat_o_cl))
    expect_true(is.matrix(simMat_o_cl))
    expect_true(is.numeric(simMat_o_cl))
    expect_true(is.matrix(
        orderSimilarityMatrix(simMat_gr, sps = spl, type = "mz", group = TRUE)))
    expect_true(is.numeric(
        orderSimilarityMatrix(simMat_gr, sps = spl, type = "mz", group = TRUE)))
    expect_equal(rownames(orderSimilarityMatrix(simMat_gr, sps = spl, 
        type = "mz", group = TRUE)), 
        c("A_Chrysin M-H     ",  "A_Apigenin M-H     ", "A_Baicalein M-H     ",
            "A_Galangin M-H     ", "A_Pinostrobin M-H     ",
            "A_Naringenin M-H     ", "A_Acacetin M-H     ", 
            "A_Kaempferol M-H     ", "A_Luteolin M-H     ",
            "A_Quercetin M-H     ", "A_Homoeriodictyol M-H     ",
            "A_Hesperetin M-H     ", "A_3-O-methylquercetin M-H     ",
            "A_Apigetrin M-H     ", "A_Avicularin M-H     ",
            "A_Swertisin M-H     ", "A_Isoquercitrin M-H     ",
            "A_2O-rhamnosyl-swertisin M-H     ", "A_Tiliroside M-H     ",
            "A_Vicenin-2 M-H     ", "A_Rutin M-H     ", "A_Hesperedin M-H     "))
    expect_equal(colnames(orderSimilarityMatrix(simMat_gr, sps = spl, 
        type = "mz", group = TRUE)), 
        c("A_Chrysin M-H     ",  "A_Apigenin M-H     ", "A_Baicalein M-H     ",
          "A_Galangin M-H     ", "A_Pinostrobin M-H     ",
          "A_Naringenin M-H     ", "A_Acacetin M-H     ", 
          "A_Kaempferol M-H     ", "A_Luteolin M-H     ",
          "A_Quercetin M-H     ", "A_Homoeriodictyol M-H     ",
          "A_Hesperetin M-H     ", "A_3-O-methylquercetin M-H     ",
          "A_Apigetrin M-H     ", "A_Avicularin M-H     ",
          "A_Swertisin M-H     ", "A_Isoquercitrin M-H     ",
          "A_2O-rhamnosyl-swertisin M-H     ", "A_Tiliroside M-H     ",
          "A_Vicenin-2 M-H     ", "A_Rutin M-H     ", "A_Hesperedin M-H     "))
})
## END unit test orderSimilarityMatrix