#'
#'
#'
#'
#'
#'
preprocPDXdata <- function(PDXdata, classifModel, topGenesDT, trainData, trainLabels) {
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
  # normDT <- as.data.table(t(normMat))[, sample := colnames(normMat)]
  # 
  # preprocDT <- merge(subtypeDT, normDT, by="sample")
  # meltPreproc
  return(subtypeDT)
}

#'
#'
#'
#'
#'
boxplotPDXsubtypePerDrug <- function(PDXmergedDT) {
  splitOnDrug <- split(PDXmergedDT, by="drug")
  names(splitOnDrug) <- unlist(lapply(splitOnDrug, function(DT) unique(DT$drug)))
  plots <- lapply(splitOnDrug, 
                  function(DT)
                    ggboxplot(DT, x="predClass", y="AAC", color="predClass", 
                              add="jitter", ylab="AAC", xlab="Subtype", 
                              pallette="jco", title=unique(DT$drug), 
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
  plotChunks <- chunk(plots, chunk.size=chunk)
  plotGrids <- lapply(plotChunks, 
                      function(plots, ncol, nrow) 
                          ggarrange(plotlist=plots, ncol=ncol, nrow=nrow),
                      ncol=ceiling(chunk/nrow),
                      nrow=nrow)
  return(plotGrids)
}