#'
#' @param
#'
#'
#' @export
rankMetaClassGenesByEffSize <- function(metaEstimateStatsDT) {
  meltDT <- melt(metaEstimateStatsDT, id.vars='gene', variable.name="class",
                 value.name="wMean")
  topGenes <- meltDT[order(-abs(wMean)), head(.SD, 1000), by=class]
  return(topGenes)
}

#'
#'
#'
#'
#'
loadGnSetGMTs <- function(gmtFldir) {
  pathwayGMTfl <- list.files(gmtFldir, pattern='.gmt', full.names=TRUE)
  pathwayL <- lapply(pathwayGMTfl, loadGSC)
  pathwayNames <- lapply(pathwayL, function(GSC) gsub('_.*$', '', GSC$addInfo[1, 1]))
  names(pathwayL) <- pathwayNames
  return(pathwayL)
}

#'
#'
#'
#'
#'
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


#'
#'
#'
#'
#'
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