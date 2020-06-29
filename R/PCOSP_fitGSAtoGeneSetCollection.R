#' Fit GSEA model to Gene Set Collection
#'
##TODO:: HEEWON - documentation
#'
#' @examples
#' data(referenceGenes)
#' data(PCOSPgenes)
#' results <- fitGSEAtoGeneSetCollection()
#'
#' @param GSC A \code{gene set collection} of data to fit the GSEA model to
#' @param referenceGenes A \code{data.frame} containing the refence genes for the
#'   GSEA model
#' @param PCOSPgenes A \code{data.frame} containing the PCSOP genes for the
#'   GSEA model
#' @param adjMethod A \code{character} vector with the multiple testing
#'    adjustment type as specified in `runGDShyper`
#' @param referenceGenes A \code{data.frame}
#' @param filePath An optional \code{character} vector containing the path to
#'     a gene set collection in .gmt format. If this is excluded the \code{GSC}
#'     parameter will be used instead.
#'
#' @return A \code{data.frame} of significant results from the model fit
#'
#' @importFrom piano loadGSC runGSAhyper
#' @export
##TODO:: Add ... to allow setting piano method parameters
fitGSAtoGeneSetCollection <- function(PCOSPgenes, GSC, referenceGenes,
                                       adjMethod, filePath)
  {
  if (missing(GSC) && !missing(filePath)) {
    geneSetCollection <- loadGSC(filePath)
  } else if (!missing(data) && missing(filePath)) {
    geneSetCollection <- data
  } else {
    stop("Please either pass in a data object or specify a file path from
         which to load the data!")
  }
  results <- runGSAhyper(PCOSPgenes,
                         gsc=geneSetCollection,
                         adjMethod=adjMethod,
                         universe=referenceGenes)
  results <- results$resTab[results$resTab[, 2] < 0.05, ]
  return(results)
}
