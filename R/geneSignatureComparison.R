#' Get the per cohort sample metaclusters and return as a data.table
#'
#' @param annotatedClusters A \code{list} of per cohort consensus clustering
#'    results annotated with sample metaclasses, as returned by the
#'    `annotateSampleMetaClasses` function in this package.
#'
#' @return A \code{data.table} of merge sample metaclass data with columns
#'    'cohort', 'samples' and 'metaClasses'
#'
#' @export
extractSampleMetaClasses <- function(annotatedClusters) {
  metaClasses <- lapply(annotatedClusters, `[[`, "metaClasses")
  sampleNames <- unlist(lapply(metaClasses, names))
  classLengths <- lapply(metaClasses, length)
  cohorts <- mapply(rep, names(annotatedClusters), classLengths, SIMPLIFY=FALSE)
  data.table("cohorts"=unlist(cohorts), "samples"=sampleNames,
             "metaClasses"=unlist(metaClasses))
}

#' Calculate the gene signature score
#'
#' @param cohortsDataL A \code{list} of gene by sample cohort gene expression
#'   matrixes.
#' @param signatureGenes A \code{character} vector of gene names in the
#'    expression signature.
#' @param sampleMetaClassDT A \code{data.table} containing the per cohort per
#'   sample meta-class predictions.
#'
#' @return A \code{data.table} of signature scores with the columns 'samples',
#'    'sigScores', 'cohorts' and 'metaClasses'.
#'
#' @export
computeSigScoreDT <- function(cohortsDataL, sampleMetaClassDT, signatureGenes) {
    cohortsSigGenes <- lapply(cohortsDataL, `[`, i=signatureGenes, j=TRUE)
    cohortsSigGenes <- lapply(cohortsSigGenes, na.omit)
    normalizedCohorts <- normalizeCohortsList(cohortsSigGenes)
    geneScoreList <- lapply(normalizedCohorts, colMeans)
    sampleNames <- unlist(lapply(geneScoreList, names))
    geneScoreDT <- data.table("samples"=sampleNames,
                              "sigScores"=unlist(geneScoreList))
    sampleClassDT <- annotateSampleMetaClassDT(sampleMetaClassDT,
                                         c('Basal', 'Exocrine', 'Classical'))
    sigScoreDT <- merge(geneScoreDT[!duplicated(samples)],
                        sampleClassDT[!duplicated(samples)],
                        on="samples")

}

#' Calculate the signature scores for a set of gene signatures
#'
#' @param geneSigL A \code{list} of gene signatures (character vector of gene
#'    names).
#' @param cohortsDataL A \code{list} of gene by sample cohort gene expression
#'   matrixes.
#' @param sampleMetaClassDT A \code{data.table} containing the per cohort per
#'   sample meta-class predictions.
#'
#' @return A \code{list} of gene signature score `data.table`s with each item
#'     returned from the `computeSigScoreDT` function from this package.
#'
#' @export
calcAllGeneSetSigScores <- function(geneSigL, cohortsDataL, sampleMetaClassDT) {
  mapply(computeSigScoreDT, geneSigL,
         MoreArgs=list(cohortsDataL=cohortsDataL,
                       sampleMetaClassDT=sampleMetaClassDT),
         SIMPLIFY=FALSE)
}

#' Z-score normalize a list of per cohort expression data
#'
#' @param cohortsList A \code{list} of gene by sample cohort gene expression
#'    matrixes.
#'
#' @export
normalizeCohortsList <- function(cohortsList) {
  .normalize <- function(cohort) t(scale(t(cohort)))
  lapply(cohortsList, .normalize)
}

#' Swap numeric metaclusters to names in the 'metaClasses' column
#'
#' @param sampleMetaClassDT A \code{data.table} containing per cohort per sample
#'     meta-class labels.
#' @param clusterLabels A \code{character} vector of meta-class labels matched
#'    according to vector index against the 'metaClasses' column.
#'
#' @return A \code{data.table} with 'metaClasses' as a character column labelled
#'    with `clusterLabels`.
#'
#' @import data.table
#' @export
annotateSampleMetaClassDT <- function(sampleMetaClassDT,
                                      clusterLabels=c("Basal", "Classical", "Exocrine")) {
  if (!is.data.table(sampleMetaClassDT)) {
    DT <- as.data.table(sampleMetaClassDT, keep.rownames='rownames')
  } else {
    DT <- copy(sampleMetaClassDT)
  }
  DT[, metaClasses := as.character(DT$metaClasses)]
  i <- 1
  for (val in na.omit(unique(DT$metaClasses))) {
    DT[metaClasses == val, metaClasses := clusterLabels[i]]
    i <- i + 1
  }
  return(na.omit(DT))
}


#' Draw a boxplot of the signature scores per sample per meta-class
#'
#' @param signatureScoreDT
#' @param comparisons A \code{list} of length two character vectors specifying
#'   the names of the metaclusters between which to compare means.
#' @param sigScoreName A \code{character} vector with the name of the gene
#'   signature being plotted. Becomes the y-axis label.
#' @param palette A \code{character} vector specifying the name of a palette
#'   from RColorBrewer for the plot. Passed as `palette` argument to
#'   `ggpubr:ggboxplot`. Defaults to 'Set1'.
#' @param saveDir An optional \code{character} vector specifying the path
#'    to the directory where the plot should be saved. If excluded, fileName
#'    will not work.
#' @param fileName An optional \code{character} vector specifying the
#'    name and extension of the file to save the plot it. This is passed to
#'    the `ggplot2::ggsave`.
#'
#' @importFrom ggpubr ggboxplot stat_compare_means
#' @importFrom ggplot2 ggsave
#' @export
plotSigScores <- function(sigScoreDT, comparisons, sigScoreName,
                          palette="Set1", saveDir, fileName) {
  plot <- ggboxplot(sigScoreDT[order(metaClasses)], x="metaClasses", y="sigScores",
                    color="metaClasses",
                    palette=palette,
                    ylab=sigScoreName,
                    add="jitter", xlab="", legend="none") +
          stat_compare_means(comparisons=comparisons)
  if (!missing(saveDir) && !missing(fileName)) {
    ggsave(plot, file=file.path(saveDir, fileName))
    return(plot)
  } else {
    return(plot)
  }
}

#' Plot a grid of gene signature score boxplots
#'
#' @param sigScoreL A \code{list} of signature scores as returned by the
#'    `calcAllGeneSigScores` function in this package.
#' @param comparison A \code{list} of character vectors containing the names
#'    of cohorts to compare. All comparisons must be pairwise (i.e., each
#'    character vector can have only two names).
#' @param palette A \code{character} vector specifying the name of a palette
#'   from RColorBrewer for the plot. Passed as `palette` argument to
#'   `ggpubr:ggboxplot`. Defaults to 'Set1'.
#' @param saveDir An optional \code{character} vector specifying the path
#'    to the directory where the plot should be saved. If excluded, fileName
#'    will not work.
#' @param fileName An optional \code{character} vector specifying the
#'    name and extension of the file to save the plot it. This is passed to
#'    the `ggplot2::ggsave`.
#'
#' @importFrom ggplot2 ggsave
#' @importFrom ggplotify as.ggplot
#' @export
plotSigScoreL <- function(sigScoreL, comparisons, palette="Set1", saveDir, fileName) {
  sigScorePlotL <- mapply(plotSigScores,
                          sigScoreDT=sigScoreL,
                          sigScoreName=names(sigScoreL),
                          MoreArgs=list(comparisons=comparisons,
                                        palette=palette),
                          SIMPLIFY=FALSE)

  sigScorePlots <- grid.arrange(grobs=sigScorePlotL)

  if (!missing(saveDir) && !missing(fileName)) {
    ggsave(plot, file=file.path(saveDir, fileName))
    return(as.ggplot(sigScorePlots))
  } else {
    return(as.ggplot(sigScorePlots))
  }
}