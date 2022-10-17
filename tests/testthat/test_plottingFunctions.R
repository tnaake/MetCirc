## create objects which will be used in unit tests
data("spectra", package = "MetCirc")

## use only a selection
condition <- c("SPL", "LIM", "ANT", "STY")
sps_tissue <- sps_tissue[c(1:20, 29:48, 113:132, 240:259), ]
similarityMat <- Spectra::compareSpectra(sps_tissue, 
    FUN = MsCoreUtils::ndotproduct, ppm = 10, m = 0.5, n = 2)
rownames(similarityMat) <- colnames(similarityMat) <- sps_tissue$name
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
linkDf <- createLinkDf(similarityMat, sps = sps_tissue, condition, 
    lower = 0.95, upper = 1)

## START unit test for plotCircos
circos.clear()
circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0),
           track.margin = c(0.0, 0))
test_that("plotCircos",  {
    expect_error(plotCircos(inds_cond, NULL, initialize = TRUE,
        featureNames = FALSE, groupSector = FALSE, groupName = FALSE,
        links = TRUE, highlight = FALSE, colour = NULL, transparency = 0.2),
        "argument is of length zero")
    expect_error(plotCircos(groupname, linkDf, initialize = TRUE,
        featureNames = FALSE, groupSector = FALSE, groupName = FALSE,
        links = TRUE, highlight = FALSE, colour = NULL, transparency = 0.2),
        "must be the same length as the vector")
    expect_error(plotCircos(featureNames = TRUE), 
        "argument \"groupname\" is missing, with no default")
    expect_error(plotCircos(inds_cond, linkDf, cexFeatureNames = 0.3,
        initialize = TRUE, featureNames = FALSE, groupSector = FALSE,
        groupName = FALSE, links = TRUE, highlight = FALSE, transparency = 0.2,
        colour = c("red", "blue")), 
        "length of colour does not match with length of group")
    expect_error(plotCircos(inds_cond, linkDf, cexFeatureNames = 0.3,
        initialize = TRUE, featureNames = FALSE, groupSector = TRUE,
        groupName = FALSE, links = TRUE, highlight = FALSE,
        colour = c("red", "blue"), transparency = 0.2), 
        "length of colour does not match with length of group")
    expect_error(plotCircos(inds_cond, linkDf, cexFeatureNames = FALSE,
        initialize = TRUE, featureNames = FALSE, groupSector = TRUE,
        groupName = FALSE, links = TRUE, highlight = FALSE, colour = NULL,
        transparency = 0.2), 
        "cexFeatureNames is not numeric")
    expect_error(plotCircos(inds_cond, linkDf, cexFeatureNames = 0.3,
        initialize = "a", featureNames = FALSE, groupSector = TRUE,
        groupName = FALSE, links = TRUE, highlight = FALSE, colour = NULL,
        transparency = 0.2),
        "initialize is not logical")
    expect_error(plotCircos(inds_cond, linkDf, cexFeatureNames = 0.3,
        initialize = TRUE, featureNames = "a", groupSector = TRUE,
        groupName = FALSE, links = TRUE, highlight = FALSE, colour = NULL,
        transparency = 0.2), 
        "featureNames is not logical")
    expect_error(plotCircos(inds_cond, linkDf, cexFeatureNames = 0.3,
        initialize = TRUE, featureNames = FALSE, groupSector = "a",
        groupName = FALSE, links = TRUE, highlight = FALSE, colour = NULL,
        transparency = 0.2),
        "groupSector is not logical")
    expect_error(plotCircos(inds_cond, linkDf, cexFeatureNames = 0.3,
        initialize = TRUE, featureNames = FALSE, groupSector = TRUE,
        groupName = "a", links = TRUE, highlight = FALSE, colour = NULL,
        transparency = 0.2),
        "groupName is not logical")
    expect_error(plotCircos(inds_cond, linkDf, cexFeatureNames = 0.3,
        initialize = TRUE, featureNames = FALSE, groupSector = TRUE,
        groupName = FALSE, links = "a", highlight = FALSE, colour = NULL,
        transparency = 0.2),
        "links is not logical")
    expect_error(plotCircos(inds_cond, linkDf, cexFeatureNames = 0.3,
        initialize = TRUE, featureNames = FALSE, groupSector = TRUE,
        groupName = FALSE, links = TRUE, highlight = "a", colour = NULL,
        transparency = 0.2),
        "highlight is not logical")
    expect_error(plotCircos(inds_cond, linkDf, cexFeatureNames = 0.3,
        initialize = TRUE, featureNames = FALSE, groupSector = TRUE,
        groupName = FALSE, links = TRUE, highlight = FALSE, colour = NULL,
        transparency = "a"),
        "transparency is not numeric")
})

## plot with different arguments are TRUE
plotCircos(inds_cond, linkDf, cexFeatureNames = 0.3, initialize = TRUE,
    featureNames = TRUE, groupSector = TRUE, groupName = FALSE, links = TRUE,
    highlight = FALSE, colour = NULL, transparency = 0.3)
plotCircos(inds_cond, linkDf, cexFeatureNames = 0.3, initialize = TRUE,
    featureNames = FALSE, groupSector = TRUE, groupName = FALSE, links = TRUE,
    highlight = FALSE, colour = c(2, 3, 4), transparency = 0.3)
plotCircos(inds_cond, linkDf, cexFeatureNames = 0.3, initialize = TRUE,
    featureNames = TRUE, groupSector = FALSE, groupName = TRUE, links = TRUE,
    highlight = FALSE, colour = NULL, transparency = 0.3)
## END unit test for plotCircos


## START unit test for highlight
test_that("highlight", {
    expect_error(highlight(inds_cond, 1, NULL, NULL, 0.4), 
        "argument is of length zero")
    expect_error(highlight(inds_cond, length(inds_cond) + 1, NULL, NULL, 0.4),
        "contains index that does not beling to available sectors")
    expect_error(highlight(inds_cond, length(inds_cond) + 1, NULL, NULL, 0.4,
        colour = "red"), 
        "contains index that does not beling to available sectors")
    expect_error(highlight(inds_cond,
        linkDf = cbind("spectrum1" = 0, "spectrum2" = 0), ind = 1),
        "must be the same length as the vector")
    expect_error(highlight(inds_cond,
        linkDf = cbind("spectrum1" = c(0, 2), "spectrum2" = c(2, 0)), ind = 1),
        "must be the same length as the vector")
    expect_error(highlight(groupname, linkDf = linkDf, links = FALSE,
        ind = 100), "contains index that does not beling to available sectors")
    expect_error(highlight(groupname, linkDf = linkDf, links = TRUE,
        ind = 100), "contains index that does not beling to available sectors")
    ## names in linkMat do not match names in groupname
    expect_error(highlight(groupname, 1, linkDf, NULL, 0.4),
        "contains index that does not beling to available sectors")
})

## plot which does not fail
highlight(inds_cond, ind = 10, linkDf, links = TRUE)
## END unit test for highlight

## START unit test for circosLegend
test_that("circosLegend", {
    expect_error(circosLegend(1, highlight = TRUE, colour = NULL), 
        "non-character argument")
    expect_error(circosLegend(1, highlight = TRUE, colour = "red"),
        "non-character argument")
    expect_error(circosLegend(1, highlight = FALSE, colour = NULL),
        "non-character argument")
    expect_error(circosLegend(1, highlight = FALSE, colour = "red"),
        "non-character argument")
})

## plot which does not fail
colours <- c("red", "blue", "green", "yellow")
circosLegend(inds_cond)
circosLegend(inds_cond, colour = colours) 
circosLegend(inds_cond, colour = colours, highlight = TRUE)
circosLegend(inds_cond, colour = colours, highlight = FALSE)
## END unit test for circosLegend


## START unit test for getLinkDfIndices
circos.clear()
## set circlize paramters
circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0),
           track.margin = c(0.0, 0))
plotCircos(inds_cond, NULL, initialize = TRUE, featureNames = FALSE,
    groupSector = FALSE, groupName = FALSE, links = FALSE, highlight = FALSE)
test_that("getLinkDfIndices",  {
    expect_equal(getLinkDfIndices(inds_cond[4], linkDf), numeric())
    expect_equal(getLinkDfIndices(inds_cond[10], linkDf), c(1, 2))
    expect_equal(getLinkDfIndices(inds_cond[17], linkDf), c(3, 1))
    expect_equal(getLinkDfIndices(inds_cond[4:12], linkDf), c(1:2))
    expect_error(getLinkDfIndices(inds_cond[1], NULL), 
        "incorrect number of dimensions")
})
## END unit test for getLinkMatrixIndices

## START unit test minFragCart2Polar
degreeFeatures <- lapply(inds_cond,
    function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
test_that("minFragCart2Polar",  {
    expect_equal(minFragCart2Polar(1, 0, degreeFeatures), 107)
    expect_equal(minFragCart2Polar(0.1, 0.9, degreeFeatures), 83)
    expect_equal(minFragCart2Polar(1, 1, degreeFeatures), NA)
    expect_equal(minFragCart2Polar(1, 0, NULL), integer())
    expect_error(minFragCart2Polar(NA, NA, degreeFeatures), 
        "missing value where TRUE/FALSE needed")
    expect_error(minFragCart2Polar(1, NA, degreeFeatures),
        "missing value where TRUE/FALSE needed")
    expect_error(minFragCart2Polar(NA, 1, degreeFeatures),
        "missing value where TRUE/FALSE needed")
})
## END unit test minFragCart2Polar


## START unit test cart2Polar
test_that("cart2Polar",  {
    expect_equal(cart2Polar(0, 0), list(r = 0, theta = 0))
    expect_equal(cart2Polar(1, 1), list(r = 1.414214, theta = 45),
        tolerance = 0.00001)
    expect_equal(cart2Polar(0, 1), list(r = 1, theta = 90))
    expect_equal(cart2Polar(-1, 1), list(r = 1.414214, theta = 135),
        tolerance = 0.00001)
    expect_equal(cart2Polar(-1, -1), list(r = 1.414214, theta = 225),
        tolerance = 0.00001)
    expect_equal(cart2Polar(1, -1), list(r = 1.414214, theta = 315),
        tolerance = 0.00001)
    expect_error(cart2Polar(NA, NA), "missing value where TRUE/FALSE needed")
    expect_error(cart2Polar(1, NA), "missing value where TRUE/FALSE needed")
    expect_error(cart2Polar(NA, 1), "missing value where TRUE/FALSE needed")
})
## END unit test cart2Polar

## START unit test plotSpectra
## plot which does not fail
gg <- plotSpectra(sps = sps_tissue, subject = "SPL_1", query = "SPL_2")

test_that("plotSpectra", {
    suppressWarnings(expect_error(
        plotSpectra(sps = NULL, subject = "LIM_1", query = "LIM_1"),
        "arguments imply differing number of rows"))
    expect_error(
        plotSpectra(sps = sps_tissue, subject = NULL, query = "LIM_1"),
        "non-character argument")
    expect_error(
        plotSpectra(sps = sps_tissue, subject = "LIM_1", query = NULL),
        "non-character argument")
    suppressWarnings(expect_error(
        plotSpectra(sps = sps_tissue, subject = "LIM_1", query = "LIM_0"),
        "arguments imply differing number of rows"))
    expect_is(gg, "gg")
})
## END unit test plotSpectra