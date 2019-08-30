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
linkDf <- createLinkDf(similarityMat, spectra_tissue, condition, 0.75, 1)
ind <- 18
linkDfInds <- getLinkDfIndices(inds_cond[ind], linkDf)

## MetCirc:::printInformationSelect(groupname = groupname, 
##  msp = NULL, ind = ind, lMatInd = linkMatInds, 
##  linkMatrixThreshold = linkMat_thr, highlight = TRUE, 
##  similarityMatrix = simMat)

## START unit test shinyCircos
test_shinyCircos <- function() {
    checkException(shinyCircos(1:2, spectra, condition))
    checkException(shinyCircos(similarityMat, NULL, condition))
    checkException(shinyCircos(similarityMat, spectra, "a"))
}
## END unit test shinyCircos


## START unit test printInformationSelect 
test_printInformationSelect <- function() {
    checkEquals(MetCirc:::printInformationSelect( 
        select=inds_cond[ind], spectra=spectra_tissue, linkDfInd=linkDfInds,  
        linkDf=linkDf, similarityMatrix=similarityMat),
    "LIM_18 (1398.71, 1018.98, , , , ) connects to  <br/>ANT_16 (0.836, 1052.94, 1020.96, , , , )<br/> ANT_17 (0.93, 1063.44, 1020.17, , , , )<br/>")
    checkException(MetCirc:::printInformationSelect( 
        select=NULL, spectra=spectra_tissue, linkDfInd=linkDfInds,  
        linkDf=linkDf, similarityMatrix=similarityMat))
    checkException(MetCirc:::printInformationSelect( 
        select=inds_cond[ind], spectra=NULL, linkDfInd=linkDfInds,  
        linkDf=linkDf, similarityMatrix=similarityMat))
    checkException(MetCirc:::printInformationSelect( 
        select=inds_cond[ind], spectra=spectra_tissue, linkDfInd=linkDfInds,  
        linkDf=NULL, similarityMatrix=similarityMat))
    checkException(MetCirc:::printInformationSelect( 
        select=inds_cond[ind], spectra=spectra_tissue, linkDfInd=linkDfInds,  
        linkDf=linkDf, similarityMatrix=NULL))
    checkEquals(MetCirc:::printInformationSelect( 
        select=inds_cond[ind], spectra=spectra_tissue, linkDfInd=numeric(),  
        linkDf=linkDf, similarityMatrix=similarityMat),
        "LIM_18 (1398.71, 1018.98, , , , ) does not connect to any feature ")
}
## END unit test printInformationSelect


