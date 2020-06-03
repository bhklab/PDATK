#' Get the top 1000 genes by meta-effect size within each meta-class.
#'
#' @param metaEstimateStatsDT A gene by meta-class \code{data.table} containing the weigthed mean
#'     of the meta-effectsize (calculated with `effsize::cohen.d`) for each
#'     predicted subtype across all cohorts and samples. As returned from the
#'     `calcClustMetaEstStats` function this this package.
#'
#' @return A \code{data.table} with the top 1000 genes per meta-class along
#'   with the associated weighted mean of the effect size and meta-class
#'   annotations.
#'
#' @import data.table
#' @export
rankMetaClassGenesByEffSize <- function(metaEstimateStatsDT) {
  meltDT <- melt(metaEstimateStatsDT, id.vars='gene', variable.name="class",
                 value.name="wMean")
  topGenes <- meltDT[order(-abs(wMean)), head(.SD, 1000), by=class]
  return(topGenes)
}

#' Read in all .gmt fiels within the specified directory as a list
#'
#' @param gmtFldir GMT file directory, the path to the directory where
#'    the .gmt files are located.
#'
#' @importFrom piano loadGSC
#' @import data.table
#' @export
loadGnSetGMTs <- function(gmtFldir) {
  pathwayGMTfl <- list.files(gmtFldir, pattern='.gmt', full.names=TRUE)
  pathwayL <- lapply(pathwayGMTfl, loadGSC)
  pathwayNames <- lapply(pathwayL, function(GSC) gsub('_.*$', '', GSC$addInfo[1, 1]))
  names(pathwayL) <- pathwayNames
  return(pathwayL)
}

#' Calculate the pathway score for a list of pathway specific gene-sets for
#'    the top 1000 genes from each meta-class.
#'
#' @param rankedMetaClassGenes A \code{data.table} with the top 1000 genes per
#'   meta-class along with the associated weighted mean of the effect size and
#'   meta-class annotations.
#' @param pathwayL A \code{list} of `GSC` objects from the `piano` package,
#'   as returned by the `loadGnSetGMTs` function in this package.
#' @param referenceGenes A \code{character} vector with the names of all common
#'   genes in the gene exoression cohorts.
#'
#' @return A \code{data.table} with the pathway score and p-value for each
#'   gene set in each pathway by mete-class.
#'
#' @importFrom piano runGSAhyper
#' @import data.table
#' @export
computePathwayScores <- function(rankedMetaClassGenes, pathwayL, referenceGenes) {
  splitByClass <- split(rankedMetaClassGenes, by='class')

  pathStats <- list()
  pvalsDT <- list()
  for (path in pathwayL) {
    pathwayName <- gsub('_.*$', '', path$addInfo[1, 1])
     pathwayScores <-
        lapply(splitByClass,
           function(DT, refGns, path)
               runGSAhyper(genes=DT$gene, gsc=path, adjMethod="fdr", universe=refGns),
           refGns=referenceGenes,
           path=path)
    pvals <- lapply(pathwayScores, `[[`, "pvalues")
    classScores <- lapply(pathwayScores,
                          function(path) qnorm(1 - (path$pvalues/2)))

    pathwayScoresDT <-
        as.data.table(do.call(cbind, classScores))[, pathway := names(classScores[[1]])]
    pvalsDT <-
        as.data.table(do.call(cbind, pvals))[, pathway := names(classScores[[1]])]

    pathScoreDT <-  melt(pathwayScoresDT, id.vars="pathway",
                          value.name="pathScore", variable.name="class",
                          variable.factor=FALSE
                          )
    pDT <- melt(pvalsDT, id.vars="pathway",
                      value.name="pval", variable.name="class",
                      variable.factor=FALSE)

    pathwayStats <- merge(pathScoreDT, pDT, by=c("class", "pathway"))
    pathStats[[pathwayName]] <- pathwayStats[, geneSet := rep(pathwayName, nrow(pathwayStats))]
  }
  pathStatsDT <- rbindlist(pathStats)[is.infinite(pathScore), pathScore := NA]
  return(pathStatsDT)
}

## TODO:: Implement palette selection for this function
#' Heatmap the gene expression signatures for each pathway in each meta-class
#'
#' @param pathStatsDT A \code{data.table} with the pathway score and p-value for each
#'   gene set in each pathway by mete-class. As returned by the
#'   `computePathwayScores` function in this package.
#' @param exclude A \code{character} vector of meta-class labels to exclude
#'   from the heatmap.
#' @param significance The alpha level to use as a cutoff for gene signatures
#'    to be included in the heatmap.
# @param palette The RColorBrewer palette to use when colouring the graph.
#'
#' @return A \code{list} of `ggplot`s, one for each pathway when calculating
#'    the pathway statistics.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom gplots heatmap.2
#' @import data.table
#' @export
heatmapPathwayScores <- function(pathStatsDT, exclude, significance=0.05) {

    if (!missing(exclude)) {
        pStatsDT <- pathStatsDT[pval < significance & !(class %in% exclude), ][, -'pval']

    } else {
        pStatsDT <- pathStatsDT[pval < significance, ][, -'pval']
    }

    pStatsDT <- dcast(pStatsDT, geneSet + pathway ~ class , value.var='pathScore')
    pStatsDT <- na.omit(pStatsDT)

    splitDT <- split(pStatsDT, by="geneSet")

    preprocMats <- lapply(splitDT, function(DT) {
        mat <- as.matrix(DT[, .SD, .SDcols=-c("geneSet", "pathway")])
        rownames(mat) <- DT$pathway
        mat
    })

    colorFun <- colorRampPalette(c("gold", "white", "black"))
    plotExpressions <- lapply(preprocMats,
                    function(mat)
                        as.expression(call('heatmap.2', x=mat, Rowv=TRUE, Colv=NA,
                                          col=colorFun(15), lhei=c(0.05, 0.20),
                                          scale="column", margins=c(5, 10),
                                          trace="none", dendrogram="row",
                                          cexRow=0.5, cexCol=1))
                    )

    plots <- lapply(plotExpressions, as.ggplot)
    return(plots)
}