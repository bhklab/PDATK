#' Compute the biomarker signature score for a gene
#'
#' @param cohortsDataL A \code{list} per cohort gene by sample gene expression
#'   matrixes.
#' @param sampleMetaClassDT A \code{data.table} containing the per cohort per
#'   sample meta-class predictions.
#' @param geneName A \code{character} vector with the name of the biomarker
#'   gene to score.
#' @param clusterLabels An optional \code{character} vector to use when
#'   annotating the sample meta-classes in the biomarker score data.table.
#'   If excluded uses the unique meta-class names from sampleMetaClassDT,
#'   sorted lexically.
#'
#' @return A \code{data.table} with the columns 'cohorts', 'samples',
#'     'metaClasses', and 'expression' containing the expression value for each
#'     cohort for the sspecified gene.
#'
#' @export
computeGeneBiomarkerScores <- function(cohortsDataL, sampleMetaClassDT, geneName,
                                       clusterLabels) {
  if (missing(clusterLabels)) {
    clusterLabels <- unique(sampleMetaClassDT[order(metaClasses)]$metaClasses)
  }
  cohortsDataL <- lapply(cohortsDataL, na.omit)

  subsetExpressionL <- lapply(cohortsDataL, `[`, i=geneName, j=TRUE)
  subsetExpressionL <- normalizeCohortsList(subsetExpressionL)
  sampleMetaClassDT <- na.omit(sampleMetaClassDT)[!duplicated(samples), ]
  subsetCohortExpression <- unlist(subsetExpressionL)
  names(subsetCohortExpression) <- unlist(lapply(cohortsDataL, colnames))

  keep <- intersect(names(subsetCohortExpression), sampleMetaClassDT$samples)
  keepDT <- sampleMetaClassDT[samples %in% keep, ]
  keepDT$expression <- subsetCohortExpression[keep]

  annotateSampleMetaClassDT(keepDT, clusterLabels)
}

#' Boxplot the normalized per sample gene expression for the specified gene,
#'     grouped by meta-class.
#'
#' @param geneBiomarkerScores A \code{data.table} containing per sample expression
#'    values for the gene specified in `geneName`. As returned by the
#'    `computeGeneBiomarkerScores` function from this package.
#' @param comparisons A \code{list} of length two `character` vectors containing
#'    the meta-classes to calculate p-value for difference in means between.
#' @param geneName A \code{character} vector with the name of the gene being
#'    boxplotted. Used to label the y-axis
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
#' @return A \code{ggplot} object of the resulting boxplot.
#'
#' @importFrom ggpubr ggboxplot stat_comapre_means
#' @importFrom ggplot2 ggsave
#' @export
boxplotBiomarkerScores <- function(geneBiomarkerScores, comparisons, geneName,
                                   palette="Set1", saveDir, fileName) {
  plot <- ggboxplot(geneBiomarkerScores[order(metaClasses)], x="metaClasses",
                    y="expression", color = "metaClasses",palette=palette,
                    add = "jitter", xlab = "", ylab=geneName,
                    legend ="none") +
          stat_compare_means(label="p.format", comparisons=comparisons)

  if (!missing(saveDir) && !missing(fileName)) {
    ggsave(file.path(saveDir, fileName), plot)
    return(plot)
  } else {
    return(plot)
  }
}

#' Create a list of ggplots comparing the meta-class expression for each biomarker
#'
#' @param geneBiomarkerscoreL  A \code{list} of \code{data.table}s containing
#'    per sample expression values for the gene specified in `geneName`.
#' @param biomarkers A \code{character} vector of biomarker gene names, with
#'    each gene named for the meta-cluster it is a biomarker in.
#' @param comparison A \code{list} of character vectors containing the names
#'    of cohorts to compare. All comparisons must be pairwise (i.e., each
#'    character vector can have only two names).
#'
#' @export
boxplotBiomarkerScoresL <- function(geneBiomarkerScoreL, biomarkers, comparisons) {
  plotL <- mapply(boxplotBiomarkerScores,
                  geneBiomarkerScores=geneBiomarkerScoreL,
                  geneName=biomarkers,
                  MoreArgs=list(comparisons=comparisons),
                  SIMPLIFY=FALSE)
  names(plotL) <- names(biomarkersOfInterest)
  return(plotL)
}
