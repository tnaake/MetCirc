## START unit test createLink0Matrix
## create objects which will be used in unit tests
data("spectra", package = "MetCirc")

## use only a selection
condition <- c("SPL", "LIM", "ANT", "STY")
sps_tissue <- sps_tissue[c(1:20, 29:48, 113:132, 240:259), ]
similarityMat <- Spectra::compareSpectra(sps_tissue, 
    FUN = ndotproduct, m = 0.5, n = 2)
colnames(similarityMat) <- rownames(similarityMat) <- sps_tissue$name
groupname <- rownames(similarityMat)
inds <- spectraCondition(sps = sps_tissue, condition = condition)
inds_match <- lapply(inds, function(x) {
    inds_match <- match(groupname, x)
    inds_match <- inds_match[!is.na(inds_match)]
    x[inds_match]
})
inds_cond <- lapply(seq_along(inds_match), function(x) {
    if (length(inds_match[[x]]) > 0) {
        paste(condition[x], inds_match[[x]], sep = "_")
    } else character()
})
inds_cond <- unique(unlist(inds_cond))
group <- unlist(lapply(strsplit(inds_cond, "_"), "[", 1))

## create link0df
link0df <- createLink0df(similarityMat, sps_tissue, condition)
ndps <- as.numeric(link0df[, "similarity"])

test_that("link0df",  {
    expect_equal(dim(link0df)[2], 5)
    expect_true(is.data.frame(link0df))
    expect_true(all(
        colnames(link0df) == c("group1", "spectrum1", "group2", "spectrum2", 
                "similarity")))
    expect_true(all(unique(c(as.character(link0df$group1), 
        as.character(link0df$group2))) %in% unique(group)))
    expect_true(all(unique(c(as.character(link0df$name1), 
        as.character(link0df$name2))) %in% unique(inds_cond)))
    expect_true(all(0 < ndps & ndps <= 1))
    simMat_mock <- similarityMat
    rownames(simMat_mock)[1] <- "a"
    expect_error(createLink0df(simMat_mock, sps_tissue, condition), 
        "colnames[(]similarityMatrix[)] != rownames[(]similarityMatrix[)]")
    simMat_mock <- similarityMat
    rownames(simMat_mock) <- NULL
    expect_error(createLink0df(simMat_mock, sps_tissue, condition),
        "subscript out of bounds")
    simMat_mock <- similarityMat
    colnames(simMat_mock) <- NULL
    expect_error(createLink0df(simMat_mock, sps_tissue, condition), "n < m")
    similarityMat <- Spectra::compareSpectra(sps_tissue[1:2], 
        FUN = ndotproduct, ppm = 10, m = 0.5, n = 2)
    colnames(similarityMat) <- rownames(similarityMat) <- sps_tissue$name[1:2]
    expect_true(
        is.data.frame(createLink0df(similarityMat, sps_tissue, condition)))
})
## END unit test link0df

## START unit test thresholdLinkDf
test_that("thresholdLinkDf", {
    expect_equal(dim(thresholdLinkDf(link0df, 0, 1)), dim(link0df))
    expect_error(thresholdLinkDf(similarityMat, 0, 1), 
        "linkDF does not have right colnames")
    expect_error(thresholdLinkDf(similarityMat, 0.6, 0.5),
        "linkDF does not have right colnames")
    expect_error(thresholdLinkDf(link0df, 1.05, 1.1), "upper greater than 1")
    expect_true(
        dim(thresholdLinkDf(link0df, 0.2, 1))[1] >=
            dim(thresholdLinkDf(link0df, 0.3, 1))[1])
    expect_error(thresholdLinkDf(link0df, 0.9, 0.1), "lower greater than upper")
    expect_true(nrow(thresholdLinkDf(link0df, 0.9, 1)) == 12)
})
## END unit test thresholdLinkDf

## START unit test createLinkDf
tLinkDf1 <- thresholdLinkDf(link0df, 0.9, 1)
tLinkDf2 <- createLinkDf(similarityMat, sps_tissue, condition, 0.9, 1)

test_that("createLinkDf",  {
    expect_true(identical(tLinkDf1, tLinkDf2))
})
## END unit test createLinkDf

## START unit test cutLinkMatrix
cutLDFInter <- cutLinkDf(tLinkDf1, type = "inter")
cutLDFIntra <- cutLinkDf(tLinkDf1, type = "intra")

test_that("cutLinkDf", {
    expect_true(all(dim(cutLinkDf(tLinkDf1, type = "all")) == dim(tLinkDf2)))
    expect_error(cutLinkDf(tLinkDf1, type = "foo"), "'arg' should be one of ")
    expect_true(all(unlist(lapply(1:dim(cutLDFInter)[1],
        function(x) cutLDFInter[x, 1] != cutLDFInter[x, 3]))))
    expect_true(all(unlist(lapply(1:dim(cutLDFIntra)[1],
        function(x) cutLDFIntra[x, 1] == cutLDFIntra[x, 3]))))
})
## END unit test cutLinkDf
