## START unit test createLink0Matrix
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
link0df <- createLink0df(similarityMat, spectra_tissue, condition)
ndps <- as.numeric(link0df[,"similarity"])

test_createLink0df <- function() {
    checkEquals(dim(link0df)[2], 5)
    checkTrue(is.data.frame(link0df))
    checkTrue(all(
        colnames(link0df) == c("group1", "spectrum1", "group2", "spectrum2", "similarity")))
    checkTrue(
        all(unique(c(as.character(link0df$group1), as.character(link0df$group2))) %in% unique(group)))
    checkTrue(
        all(unique(c(as.character(link0df$name1), as.character(link0df$name2))) %in% unique(inds_cond)))
    checkTrue(all(0 < ndps & ndps <= 1))
}
## END unit test link0df

## START unit test thresholdLinkDf
test_thresholdLinkDf <- function() {
    checkEquals(dim(thresholdLinkDf(link0df, 0, 1)), dim(link0df))
    checkException(thresholdLinkDf(similarityMat, 0, 1))
    checkException(thresholdLinkDf(similarityMat, 0.6, 0.5))
    checkException(thresholdLinkDf(link0df, 1.05, 1.1))
    checkTrue(
        dim(thresholdLinkDf(link0df, 0.2, 1))[1] >= 
            dim(thresholdLinkDf(link0df, 0.3, 1))[1])
}
## END unit test thresholdLinkDf

## START unit test createLinkDf
tLinkDf1 <- thresholdLinkDf(link0df, 0.9, 1)
tLinkDf2 <- createLinkDf(similarityMat, spectra_tissue, condition, 0.9, 1)

test_createLinkDf <- function() {
    checkTrue(identical(tLinkDf1, tLinkDf2))
}
## END unit test createLinkDf

## START unit test cutLinkMatrix
cutLDFInter <- cutLinkDf(tLinkDf1, type="inter")
cutLDFIntra <- cutLinkDf(tLinkDf1, type="intra")

test_cutLinkDf <- function() {
    checkTrue(
        all(dim(cutLinkDf(tLinkDf1, type = "all")) == dim(tLinkDf2))
    )
    checkException(cutLinkDf(tLinkDf1, type = "foo"))
    checkTrue(
        all(unlist(lapply(1:dim(cutLDFInter)[1], 
                          function(x) cutLDFInter[x,1] != cutLDFInter[x,3])))
    )
    checkTrue(
        all(unlist(lapply(1:dim(cutLDFIntra)[1], 
                          function(x) cutLDFIntra[x,1] == cutLDFIntra[x,3])))
    )
}
## END unit test cutLinkDf
