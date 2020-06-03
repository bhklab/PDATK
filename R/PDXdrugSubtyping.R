#' Preprocess PDX gene expression data
#'
#' @param PDXdata A \code{data.table}
#' @param classifModel A trained \code{randomForest} classifier, as returned
#'    by the `randomForest::randomForest` function.
#' @param topGenesDT A \code{data.table} with the columns 'topGenes1'
#'    and 'topGenes2' containing the top scoring pairs, 'genePairs' containing
#'    a label for each pair, and 'classes' indicating the meta-class each gene pair.
#'    is associated with.
#' @param trainData A gene-pair by sample binary \code{matrix} where 1 represents
#'    the predicted top scoring pair being correctly ordered, as returned
#'    by the `calcTopGenes1Gt2Matrix` function in this pacakge.
#' @param trainLabels A \code{factor} vector with the meta-class for each
#'    sample in the training data.
#'
#' @return A \code{data.table} of per sample per drug AAC values, labelled
#'    with sample name and predicted meta-class/subtype.
#'
#' @import data.table
#' @export
preprocPDXdata <- function(PDXdata, classifModel, topGenesDT, trainData,
                           trainLabels) {
  # Deal with matrix and df input
  if (is.matrix(PDXdata)) {
    PDXdata <- as.data.table(PDXdata)[, gene_name := rownames(PDXdata)]
  } else if (!is.data.table(PDXdata)) {
    if (!grepl('gene', colnames(PDXdata), ignore.case=TRUE))
      PDXdata$gene_name <- rownames(PDXdata)
    PDXdata <- as.data.table(PDXdata)
  }
  # Get remove annotation columns and find gene name column
  whichNumeric <- lapply(PDXdata, class) %in% c("numeric")
  IDcol <- grep('gene', colnames(PDXdata), value=TRUE, ignore.case=TRUE)[1]

  # Convert to matrix, log and average replicates
  dataMat <- as.matrix(PDXdata[, .SD, .SDcols=whichNumeric])
  rownames(dataMat) <- PDXdata[[IDcol]]
  normMat <- avereps(log2(dataMat + 1))

  # Predict subtypes
  subtypeDT <- predictSampleMetaClass(normMat, classifModel, topGenesDT,
                                      trainData, trainLabels)
  return(subtypeDT)
}

#' Boxplot per drug AAC values between subclasses
#'
#' @param PDXmergedDT A \code{data.table} of merge patient derived xenograph
#'   drug sensitivity and predicted meta-class/subtype.
#'
#' @return A \code{list} of ggplot objects, one for each drug in `PDXmergedDT`.
#'    This can be format to a plot grid using the `ggarrangePlotL` function.
#'
#' @importFrom ggpubr ggboxplot stat_compare_means
#' @import data.table
#' @export
boxplotPDXsubtypePerDrug <- function(PDXmergedDT) {
  splitOnDrug <- split(PDXmergedDT, by="drug")
  names(splitOnDrug) <- unlist(lapply(splitOnDrug, function(DT) unique(DT$drug)))
  plots <- lapply(splitOnDrug,
                  function(DT)
                    ggboxplot(DT, x="predClass", y="AAC", color="predClass",
                              add="jitter", ylab="AAC", xlab="Subtype",
                              pallette="Set1", title=unique(DT$drug),
                              legend="none") +
                              stat_compare_means(method="kruskal.test")
                              )
  return(plots)
}

#' Take a long list of plots and chunk it into a list of plot chunck x chunk
#'    plot grids.
#'
#' @param plots A \code{list} of grob or ggplot objects which can be
#'     passed to the `ggpubr::ggarrange` function.
#' @param chunk A \code{numeric} vector indicating the integer number
#'     of plots per plot grid. All returned grids have equal ncol and nrow (
#'     i.e., are ~ square).
#'
#' @importFrom BBmisc chunk
#' @export
ggarrangePlotL <- function(plotL, chunk) {
  nrow <- .ceilSqrt(chunk)
  plotChunks <- chunk(plotL, chunk.size=chunk)
  plotGrids <- lapply(plotChunks,
                      function(plots, ncol, nrow)
                          ggarrange(plotlist=plots, ncol=ncol, nrow=nrow),
                      ncol=ceiling(chunk/nrow),
                      nrow=nrow)
  return(plotGrids)
}