## load ggplot2
library(ggplot2)

## create objects which will be used in unit tests
data("spectra", package="MetCirc")
## use only a selection 
condition <- c("SPL", "LIM", "ANT", "STY")
spectra_tissue <- spectra_tissue[c(1:20, 29:48, 113:132, 240:259),]
similarityMat <- compare_Spectra(spectra_tissue, fun=normalizeddotproduct)  
groupname <- rownames(similarityMat)
inds <- MetCirc:::spectraCond(spectra_tissue, condition=condition)
inds_match <- lapply(inds, function(x) {inds_match <- match(groupname, x)
inds_match <- inds_match[!is.na(inds_match)]; x[inds_match]})
inds_cond <- lapply(seq_along(inds_match), 
    function(x) {
        if (length(inds_match[[x]]) > 0) {
            paste(condition[x], inds_match[[x]], sep="_")
        } else character()
})
inds_cond <- unique(unlist(inds_cond))
group <- unlist(lapply(strsplit(inds_cond, "_"), "[", 1))

## create link0df
linkDf <- createLinkDf(similarityMat, spectra_tissue, condition, lower=0.95, upper=1)

## START unit test for plotCircos
circos.clear()
circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), 
           track.margin = c(0.0, 0))
test_plotCircos <- function() {
    checkException(plotCircos(inds_cond, NULL, initialize=TRUE, 
        featureNames=FALSE, groupSector=FALSE, groupName=FALSE, links=TRUE, 
        highlight=FALSE, colour=NULL, transparency=0.2))
    checkException(plotCircos(groupname, linkDf, initialize=TRUE, 
        featureNames=FALSE, groupSector=FALSE, groupName=FALSE, links=TRUE, 
        highlight=FALSE, colour=NULL, transparency=0.2)) ## names are different
    checkException(plotCircos(featureNames=TRUE))
    checkException(plotCircos(inds_cond, linkDf, cexFeatureNames=0.3,
        initialize=TRUE, featureNames=FALSE, groupSector=FALSE, groupName=FALSE, 
        links=TRUE, highlight=FALSE, transparency=0.2, colour=c("red", "blue")))
    checkException(plotCircos(inds_cond, linkDf, cexFeatureNames=0.3,
        initialize=TRUE, featureNames=FALSE, groupSector=TRUE, groupName=FALSE, 
        links=TRUE, highlight=FALSE, colour=c("red", "blue"), transparency=0.2))
    checkException(plotCircos(inds_cond, linkDf, cexFeatureNames=FALSE,
        initialize=TRUE, featureNames=FALSE, groupSector=TRUE, groupName=FALSE, 
        links=TRUE, highlight=FALSE, colour=NULL, transparency=0.2))
    checkException(plotCircos(inds_cond, linkDf, cexFeatureNames=0.3,
        initialize="a", featureNames=FALSE, groupSector=TRUE, groupName=FALSE, 
        links=TRUE, highlight=FALSE, colour=NULL, transparency=0.2))
    checkException(plotCircos(inds_cond, linkDf, cexFeatureNames=0.3,
        initialize=TRUE, featureNames="a", groupSector=TRUE, groupName=FALSE, 
        links=TRUE, highlight=FALSE, colour=NULL, transparency=0.2))
    checkException(plotCircos(inds_cond, linkDf, cexFeatureNames=0.3,
        initialize=TRUE, featureNames=FALSE, groupSector="a", groupName=FALSE, 
        links=TRUE, highlight=FALSE, colour=NULL, transparency=0.2))
    checkException(plotCircos(inds_cond, linkDf, cexFeatureNames=0.3,
        initialize=TRUE, featureNames=FALSE, groupSector=TRUE, groupName="a", 
        links=TRUE, highlight=FALSE, colour=NULL, transparency=0.2))
    checkException(plotCircos(inds_cond, linkDf, cexFeatureNames=0.3,
        initialize=TRUE, featureNames=FALSE, groupSector=TRUE, groupName=FALSE, 
        links="a", highlight=FALSE, colour=NULL, transparency=0.2))
    checkException(plotCircos(inds_cond, linkDf, cexFeatureNames=0.3,
        initialize=TRUE, featureNames=FALSE, groupSector=TRUE, groupName=FALSE, 
        links=TRUE, highlight="a", colour=NULL, transparency=0.2))
    checkException(plotCircos(inds_cond, linkDf, cexFeatureNames=0.3,
        initialize=TRUE, featureNames=FALSE, groupSector=TRUE, groupName=FALSE, 
        links=TRUE, highlight=FALSE, colour=NULL, transparency="a"))
}

## plot with different arguments are TRUE
plotCircos(inds_cond, linkDf, cexFeatureNames=0.3,
           initialize=FALSE, featureNames=TRUE, groupSector=TRUE, groupName=FALSE, 
           links=TRUE, highlight=FALSE, colour=NULL, transparency=0.3)
plotCircos(inds_cond, linkDf, cexFeatureNames=0.3,
           initialize=FALSE, featureNames=FALSE, groupSector=TRUE, groupName=FALSE, 
           links=TRUE, highlight=FALSE, colour=c(1,2,3,4), transparency=0.3)
plotCircos(inds_cond, linkDf, cexFeatureNames=0.3,
           initialize=FALSE, featureNames=TRUE, groupSector=FALSE, groupName=TRUE, 
           links=TRUE, highlight=FALSE, colour=NULL, transparency=0.3)
## END unit test for plotCircos


## START unit test for highlight
test_highlight <- function() {
    checkException(highlight(inds_cond, 1, NULL, NULL, 0.4))
    checkException(highlight(inds_cond, length(inds_cond)+1,NULL, NULL, 0.4))
    checkException(highlight(inds_cond, length(inds_cond)+1,NULL, NULL, 0.4,
                             colour = "red"))
    checkException(highlight(inds_cond, 
        linkDf=cbind("spectrum1"=0, "spectrum2"=0), ind=1))
    checkException(highlight(inds_cond, 
        linkDf=cbind("spectrum1"=c(0, 2), "spectrum2"=c(2, 0)), ind=1))
    checkException(highlight(groupname, linkDf=linkDf, links=FALSE, ind=100))
    checkException(highlight(groupname, linkDf=linkDf, links=TRUE, ind=100))
    ## names in linkMat do not match names in groupnameO
    checkException(highlight(groupname, 1, linkDf, NULL, 0.4))
}

## plot which does not fail
highlight(inds_cond, ind=10, linkDf, links=TRUE)
## END unit test for highlight

## START unit test for circosLegend
test_circosLegend <- function() {
    checkException(circosLegend(1, highlight=TRUE, colour=NULL))
    checkException(circosLegend(1, highlight=TRUE, colour="red"))
    checkException(circosLegend(1, highlight=FALSE, colour=NULL))
    checkException(circosLegend(1, highlight=FALSE, colour="red"))
}

## plot which does not fail
circosLegend(inds_cond)
circosLegend(inds_cond, colour=1:4)
circosLegend(inds_cond, colour=1:4, highlight=TRUE)
circosLegend(inds_cond, colour=1:4, highlight=FALSE)
## END unit test for circosLegend


## START unit test for getLinkDfIndices
circos.clear()
## set circlize paramters
circos.par(gap.degree = 0, cell.padding = c(0.0, 0, 0.0, 0), 
           track.margin = c(0.0, 0))
plotCircos(inds_cond, NULL, initialize = TRUE, 
    featureNames = FALSE, groupSector = FALSE, groupName = FALSE, links = FALSE, highlight = FALSE)
test_getLinkDfIndices <- function() {
    checkEquals(getLinkDfIndices(inds_cond[4], linkDf), numeric())
    checkEquals(getLinkDfIndices(inds_cond[9], linkDf), 1)
    checkEquals(getLinkDfIndices(inds_cond[12], linkDf), c(2, 3))
    checkEquals(getLinkDfIndices(inds_cond[4:12], linkDf), c(1:3))
    checkException(getLinkDfIndices(inds_cond[1], NULL))
}
## END unit test for getLinkMatrixIndices

## START unit test minFragCart2Polar

degreeFeatures <- lapply(inds_cond, 
    function(x) mean(circlize:::get.sector.data(x)[c("start.degree", "end.degree")]))
test_minFragCart2Polar <- function() {
    checkEquals(minFragCart2Polar(1,0,degreeFeatures), 106)
    checkEquals(minFragCart2Polar(0.1,0.9,degreeFeatures), 82)
    checkEquals(minFragCart2Polar(1,1, degreeFeatures), NA)
    checkEquals(minFragCart2Polar(1, 0, NULL), integer())
    checkException(minFragCart2Polar(NA, NA, degreeFeatures))
    checkException(minFragCart2Polar(1, NA, degreeFeatures))
    checkException(minFragCart2Polar(NA, 1, degreeFeatures))
}
## END unit test minFragCart2Polar


## START unit test cart2Polar
test_cart2Polar <- function() {
    checkEquals(cart2Polar(0, 0), list(r = 0, theta = 0))
    checkEquals(cart2Polar(1, 1), list(r = 1.414214, theta = 45), 
        tolerance = 0.00001)
    checkEquals(cart2Polar(0, 1), list(r = 1, theta = 90))
    checkEquals(cart2Polar(-1, 1), list(r = 1.414214, theta = 135),
        tolerance = 0.00001)
    checkEquals(cart2Polar(-1, -1), list(r = 1.414214, theta = 225), 
        tolerance = 0.00001)
    checkEquals(cart2Polar(1, -1), list(r = 1.414214, theta = 315), 
        tolerance = 0.00001)
    checkException(cart2Polar(NA, NA))
    checkException(cart2Polar(1, NA))
    checkException(cart2Polar(NA, 1))
}
## END unit test cart2Polar

## START unit test plotSpectra
## plot which does not fail
gg <- plotSpectra(spectra_tissue, subject="SPL_1", query="SPL_2")

test_plotSpectra <- function() {
    checkException(plotSpectra(NULL, "LIM_1", "LIM_1"))
    checkException(plotSpectra(spectra_tissue, NULL, "LIM_1"))
    checkException(plotSpectra(spectra_tissue, "LIM_1", NULL))
    checkException(plotSpectra(spectra_tissue, "LIM_1", "LIM_0"))
    checkEquals(is(gg), "gg")
}


## END unit test plotSpectra