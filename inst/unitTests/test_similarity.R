data("convertMsp2Spectra", package="MetCirc")
spl <- convertMsp2Spectra(msp2spectra)

simMat <- compare_Spectra(spl[1:2], fun="dotproduct")

## START unit test compare_Spectra
test_compare_Spectra <- function() {
    checkEquals(simMat, matrix(c(1, 1.004575e-05, 1.004575e-05, 1), ncol=2, dimnames=list(1:2, 1:2)), tolerance=1e-6)
    checkException(compare_Spectra(spl[1], "dotproduct"))
    checkException(compare_Spectra(spl[1:2], "dotproduct2"))
    
}
## END unit test compare_Spectra

## START unit test normalizeddotproduct
ndp <- normalizeddotproduct(spl[[1]], spl[[2]])
test_normalizeddotproduct <- function() {
    checkEquals(ndp, 1.591002e-07, tolerance=1e-6)
    checkException(normalizeddotproduct(spl[[1]], 1))
    checkException(normalizeddotproduct(1, spl[[1]]))
    checkException(normalizeddotproduct(spl[[1]], spl[[2]], "a"))
    checkException(normalizeddotproduct(spl[[1]], spl[[2]], 1, "a"))
}
## END unit test normalizeddotproduct

## START unit test neutralloss
nl <- neutralloss(spl[[1]], spl[[2]])
test_neutralloss <- function() {
    checkEquals(nl, 5.899399e-05, tolerance=1e-6)
    checkException(neutralloss(spl[[1]], 1))
    checkException(neutralloss(1, spl[[1]]))
    checkException(neutralloss(spl[[1]], spl[[2]], "a"))
    checkException(neutralloss(spl[[1]], spl[[2]], 1, "a"))
}
## END unit test neutralloss


## START unit test createOrderedSimMat
simMat_o_mz <- orderSimilarityMatrix(simMat, spl, type="mz", group=FALSE)
simMat_o_rt <- orderSimilarityMatrix(simMat, spl, type="retentionTime", group=FALSE)
simMat_o_cl <- orderSimilarityMatrix(simMat, spl, type="clustering", group=FALSE)

## create a matrix with groups
simMat_gr <- simMat
rownames(simMat_gr) <- colnames(simMat_gr) <- paste("A", rownames(simMat), sep="_")

test_createOrderedSimMat <- function() {
    checkException(orderSimilarityMatrix(simMat, spl, type="foo"))
    checkException(orderSimilarityMatrix(simMat, spl, type="mz", group="a"))
    checkException(orderSimilarityMatrix(type="mz"))
    checkException(orderSimilarityMatrix(type="retentionTime"))
    checkException(orderSimilarityMatrix(type="clustering"))
    checkEquals(colnames(simMat), colnames(simMat_o_mz))
    checkEquals(rownames(simMat), rownames(simMat_o_mz))
    checkEquals(colnames(simMat_o_mz), rownames(simMat_o_mz))
    checkEquals(dim(simMat), dim(simMat_o_mz))
    checkTrue(is.matrix(simMat_o_mz))
    checkTrue(is.numeric(simMat_o_mz))
    checkEquals(colnames(simMat_o_rt), rownames(simMat_o_rt))
    checkEquals(dim(simMat), dim(simMat_o_rt))
    checkTrue(is.matrix(simMat_o_rt))
    checkTrue(is.numeric(simMat_o_rt))
    checkEquals(colnames(simMat_o_cl), rownames(simMat_o_cl))
    checkEquals(dim(simMat), dim(simMat_o_cl))
    checkTrue(is.matrix(simMat_o_cl))
    checkTrue(is.numeric(simMat_o_cl))
    checkTrue(is.matrix(
        orderSimilarityMatrix(simMat_gr, spl, type="mz", group=TRUE)))
    checkTrue(is.numeric(
        orderSimilarityMatrix(simMat_gr, spl, type="mz", group=TRUE)))
    checkEquals(
        rownames(orderSimilarityMatrix(simMat_gr, spl, type="mz", group=TRUE)),
        c("A_1", "A_2"))
    checkEquals(
        colnames(orderSimilarityMatrix(simMat_gr, spl, type="mz", group=TRUE)),
        c("A_1", "A_2"))
}
## END unit test createOrderedSimMat