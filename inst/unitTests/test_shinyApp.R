## create objects which will be used in unit tests
data("spectra", package = "MetCirc")
## use only a selection
condition <- c("SPL", "LIM", "ANT", "STY")
spectra_tissue <- spectra_tissue[c(1:20, 29:48, 113:132, 240:259),]
similarityMat <- compare_Spectra(spectra_tissue, fun = normalizeddotproduct)
groupname <- rownames(similarityMat)
inds <- MetCirc:::spectraCond(spectra_tissue, condition = condition)
inds_match <- lapply(inds, function(x) {inds_match <- match(groupname, x)
inds_match <- inds_match[!is.na(inds_match)]; x[inds_match]})
inds_cond <- lapply(seq_along(inds_match),
    function(x) {
        if (length(inds_match[[x]]) > 0) {
            paste(condition[x], inds_match[[x]], sep = "_")
        } else character()
})
inds_cond <- unique(unlist(inds_cond))
group <- unlist(lapply(strsplit(inds_cond, "_"), "[", 1))

## create link0df
linkDf <- createLinkDf(similarityMat, spectra_tissue, condition, 0.75, 1)
ind <- 18
linkDfInds <- getLinkDfIndices(inds_cond[ind], linkDf)

## START unit test shinyCircos
test_shinyCircos <- function() {
    checkException(shinyCircos(1:2, spectra, condition))
    checkException(shinyCircos(similarityMat, NULL, condition))
    checkException(shinyCircos(similarityMat, spectra, "a"))
}
## END unit test shinyCircos


## START unit test printInformationSelect
test_printInformationSelect <- function() {
    checkException(MetCirc:::printInformationSelect(
        select = NULL, spectra = spectra_tissue, linkDfInd = linkDfInds,
        linkDf = linkDf, similarityMatrix = similarityMat))
    checkException(MetCirc:::printInformationSelect(
        select = inds_cond[ind], spectra = NULL, linkDfInd = linkDfInds,
        linkDf = linkDf, similarityMatrix = similarityMat))
    checkException(MetCirc:::printInformationSelect( 
        select = inds_cond[ind], spectra = spectra_tissue,
        linkDfInd = linkDfInds, linkDf = NULL,
        similarityMatrix = similarityMat))
    checkException(MetCirc:::printInformationSelect(
        select = inds_cond[ind], spectra = spectra_tissue,
        linkDfInd = linkDfInds, linkDf = linkDf, similarityMatrix = NULL))
    checkEquals(MetCirc:::printInformationSelect(
        select = inds_cond[ind], spectra = spectra_tissue,
        linkDfInd = numeric(),  linkDf = linkDf,
        similarityMatrix = similarityMat),
        "LIM_18 (1398.71, 1018.98, , , , ) does not connect to any feature ")
}
## END unit test printInformationSelect


## START unit test typeMatch_link0
MZ <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMat,
    spectra = spectra_tissue, type = "mz", condition = condition)
link0df_mz <- MZ[["link0df"]]
mz_match <- MZ[["type_match"]]
RT <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMat,
    spectra = spectra_tissue, type = "retentionTime", condition = condition)
link0df_rt <- RT[["link0df"]]
rt_match <- RT[["type_match"]]
Clust <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMat,
    spectra = spectra_tissue, type = "clustering", condition = condition)
link0df_clust <- Clust[["link0df"]]
clust_match <- Clust[["type_match"]]

test_typeMatch_link0 <- function() {
    checkEquals(length(mz_match), 106)
    checkEquals(length(mz_match), length(rt_match))
    checkEquals(length(mz_match), length(clust_match))
    checkEquals(dim(link0df_mz), c(5521, 5))
    checkEquals(dim(link0df_mz), dim(link0df_rt))
    checkEquals(dim(link0df_mz), dim(link0df_clust))
    checkTrue(is.character(mz_match))
    checkTrue(is.character(rt_match))
    checkTrue(is.character(clust_match))
    checkTrue(is.data.frame(link0df_mz))
    checkTrue(is.data.frame(link0df_rt))
    checkTrue(is.data.frame(link0df_clust))
    checkException(MetCirc:::typeMatch_link0(similarityMatrix = NULL,
        spectra = spectra_tissue, type = "mz", condition = condition))
    checkException(MetCirc:::typeMatch_link0(similarityMatrix = similarityMat,
        spectra = NULL, type = "mz", condition = condition))
    checkException(MetCirc:::typeMatch_link0(similarityMatrix = similarityMat,
        spectra = spectra_tissue, type = "a", condition = condition))
    checkException(MetCirc:::typeMatch_link0(similarityMatrix = similarityMat,
        spectra = spectra_tissue, type = "mz", condition = "abc"))
}
## END unit test typeMatch_link0


## START unit test recordPlotFill_degreeFeatures
plotFill_degree_mz <- MetCirc:::recordPlotFill_degreeFeatures(mz_match)
degree_mz <- plotFill_degree_mz[["degreeFeatures"]]
plotFill_mz <- plotFill_degree_mz[["plotFill"]]
plotFill_degree_rt <- MetCirc:::recordPlotFill_degreeFeatures(rt_match)
degree_rt <- plotFill_degree_rt[["degreeFeatures"]]
plotFill_rt <- plotFill_degree_rt[["plotFill"]]
plotFill_degree_clust <- MetCirc:::recordPlotFill_degreeFeatures(clust_match)
degree_clust <- plotFill_degree_clust[["degreeFeatures"]]
plotFill_clust <- plotFill_degree_clust[["plotFill"]]

test_recordPlotFill_degreeFeatures <- function() {
    checkTrue(is.list(degree_mz))
    checkTrue(is.list(degree_rt))
    checkTrue(is.list(degree_clust))
    is.numeric(unlist(degree_mz))
    is.numeric(unlist(degree_rt))
    is.numeric(unlist(degree_clust))
    checkEquals(length(degree_mz), length(mz_match))
    checkEquals(length(degree_rt), length(rt_match))
    checkEquals(length(degree_clust), length(clust_match))
    checkTrue(class(plotFill_mz) == "recordedplot")
    checkTrue(class(plotFill_rt) == "recordedplot")
    checkTrue(class(plotFill_clust) == "recordedplot")
    checkException(MetCirc:::recordPlotFill_degreeFeatures(c("1", "1", "1")))
}
## END unit test recordPlotFill_degreeFeatures


## START unit test recordPlotHighlight
highlightMz <- MetCirc:::recordPlotHighlight(mz_match)
highlightRt <- MetCirc:::recordPlotHighlight(rt_match)
highlightClust <- MetCirc:::recordPlotHighlight(clust_match)

test_recordPlotHighlight <- function() {
    checkTrue(class(highlightMz) == "recordedplot")
    checkTrue(class(highlightRt) == "recordedplot")
    checkTrue(class(highlightClust) == "recordedplot")
}
## END unit test recordPlotHighlight


## START unit test replayPlotOrder
plot_l <- list(highlightMz = highlightMz, highlightRT = highlightRt,
               highlightClust = highlightClust, fillMz = plotFill_mz,
               fillRT = plotFill_rt, fillClust = plotFill_clust)
plot_l_mock <- list(highlightMz = NULL, highlightRT = highlightRt,
               highlightClust = highlightClust, fillMz = plotFill_mz,
               fillRT = plotFill_rt, fillClust = plotFill_clust)
## plots that work
MetCirc:::replayPlotOrder(orderMatch = "mz", onCircle = FALSE,
    plot_l = plot_l, ind = NULL) ## fill mz
MetCirc:::replayPlotOrder(orderMatch = "retentionTime",
    onCircle = FALSE, plot_l = plot_l, ind = NULL) ## fill rt
MetCirc:::replayPlotOrder(orderMatch = "clustering",
    onCircle = FALSE, plot_l = plot_l, ind = NULL) ## fill clust
MetCirc:::replayPlotOrder(orderMatch = "mz", onCircle = FALSE,
    plot_l = plot_l, ind = 1) ## highlight mz
MetCirc:::replayPlotOrder(orderMatch = "retentionTime",
    onCircle = FALSE, plot_l = plot_l, ind = 1) ## highlight rt
MetCirc:::replayPlotOrder(orderMatch = "clustering",
    onCircle = FALSE, plot_l = plot_l, ind = 1) ## highlight clust

test_replayPlotOrder <- function() {
    checkException(MetCirc:::replayPlotOrder(orderMatch = "mz",
        onCircle = FALSE, plot_l = plot_l[-1], ind = 1))
    checkException(MetCirc:::replayPlotOrder(orderMatch = "mz",
        onCircle = FALSE, plot_l=plot_l_mock, ind=1))
    checkException(MetCirc:::replayPlotOrder(orderMatch = "mz",
        onCircle = "a", plot_l = plot_l, ind = 1))
    checkException(MetCirc:::replayPlotOrder(orderMatch = "a",
        onCircle = FALSE, plot_l = plot_l, ind = 1))
}
## END unit test replayPlotOrder


## START unit test replayPlotAdd
## order according to retention time
mz_match <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMat,
    spectra = spectra_tissue, type = "mz", condition = condition)
linkDf <- mz_match[["link0df"]]
mz_match <- mz_match[["type_match"]]
rt_match <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMat,
    spectra = spectra_tissue, type = "retentionTime", condition = condition)
rt_match <- rt_match[["type_match"]]
clust_match <- MetCirc:::typeMatch_link0(similarityMatrix = similarityMat,
    spectra = spectra_tissue, type = "clustering", condition = condition)
clust_match <- clust_match[["type_match"]]
circos.clear()
circos.initialize(factor(mz_match, levels = mz_match),
     xlim = matrix(rep(c(0,1), length(mz_match)), ncol = 2, byrow = TRUE))
circos.trackPlotRegion(factor(mz_match, levels = mz_match), ylim = c(0,1))
MetCirc:::replayPlotOrder(orderMatch = "mz",
        onCircle = FALSE, plot_l = plot_l, ind = NULL)
MetCirc:::replayPlotAdd(orderMatch = "mz", onCircle = FALSE, linkDf = linkDf,
     mz_match = mz_match, rt_match = rt_match, clust_match = clust_match,
     ind = 1, indMz = NULL, indRT = NULL, indCluster = NULL)

test_replayPlotAdd <- function() {
    checkException(MetCirc:::replayPlotAdd(orderMatch = "a", onCircle = FALSE,
        linkDf = link0df_mz, mz_match = mz_match, rt_match = rt_match,
        clust_match = clust_match, ind = 1, indMz = NULL, indRT = NULL,
        indCluster = NULL))
    checkException(MetCirc:::replayPlotAdd(orderMatch = "mz", onCircle = "a",
        linkDf = link0df_mz, mz_match = mz_match, rt_match = rt_match,
        clust_match = clust_match, ind = 1, indMz = NULL, indRT = NULL,
        indCluster = NULL))
    checkException(MetCirc:::replayPlotAdd(orderMatch = "mz", onCircle = TRUE,
        linkDf = NULL, mz_match = mz_match, rt_match = rt_match,
        clust_match = clust_match, ind = 1, indMz = NULL, indRT = NULL,
        indCluster = NULL))
    checkException(MetCirc:::replayPlotAdd(orderMatch = "mz", onCircle = FALSE,
        linkDf = link0df_mz, mz_match = NULL, rt_match = rt_match,
        clust_match = clust_match, ind = 1, indMz = NULL, indRT = NULL,
        indCluster = NULL))
    checkException(MetCirc:::replayPlotAdd(orderMatch = "retentionTime",
        onCircle = FALSE, linkDf = link0df_mz, mz_match = mz_match,
        rt_match = NULL, clust_match = clust_match, ind = 1, indMz = NULL,
        indRT = NULL, indCluster = NULL))
    checkException(MetCirc:::replayPlotAdd(orderMatch = "clustering",
        onCircle = FALSE, linkDf = link0df_mz, mz_match = mz_match,
        rt_match = rt_match, clust_match = NULL, ind = 1, indMz = NULL,
        indRT = NULL, indCluster = NULL))
    checkException(MetCirc:::replayPlotAdd(orderMatch = "mz", onCircle = TRUE,
        linkDf = link0df_mz, mz_match = mz_match, rt_match = rt_match,
        clust_match = clust_match, ind = "a", indMz = NULL, indRT = NULL,
        indCluster = NULL))
    checkException(MetCirc:::replayPlotAdd(orderMatch = "mz", onCircle = FALSE,
        linkDf = link0df_mz, mz_match = mz_match, rt_match = rt_match,
        clust_match = clust_match, ind = 1, indMz = "a", indRT = NULL,
        indCluster = NULL))
    checkException(MetCirc:::replayPlotAdd(orderMatch = "retentionTime",
        onCircle = FALSE, linkDf = link0df_mz, mz_match = mz_match,
        rt_match = rt_match, clust_match = clust_match, ind = 1, indMz = NULL,
        indRT = "a", indCluster = NULL))
    checkException(MetCirc:::replayPlotAdd(orderMatch = "clustering",
        onCircle = FALSE, linkDf = link0df_mz, mz_match = mz_match,
        rt_match = rt_match, clust_match = clust_match, ind = 1, indMz = NULL,
        indRT = NULL, indCluster = "a"))
}
## END unit test replayPlotAdd

## START unit test select
mz <- 1
rt <- 2
clust <- 3
test_select <- function() {
    checkEquals(MetCirc:::select(condition = "mz", mz = mz, rt = rt,
        clust = clust), 1)
    checkEquals(MetCirc:::select(condition = "retentionTime", mz = mz, rt = rt,
        clust = clust), 2)
    checkEquals(MetCirc:::select(condition = "clustering", mz = mz, rt = rt,
        clust = clust), 3)
    checkException(MetCirc:::select(condition = "a", mz = mz, rt = rt,
        clust = clust))
    checkException(MetCirc:::select(condition = "mz", rt = rt,
        clust = clust))
    checkException(MetCirc:::select(condition = "retentionTime", mz = mz,
        clust = clust))
    checkException(MetCirc:::select(condition = "clustering", mz = mz, rt = rt))
}
## END unit test select

