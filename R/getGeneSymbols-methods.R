#' Extract gene symbols from a SummarizedExperiment object
#'
#' @param dataset A \code{SummarizedExperiment} object to extract gene symbols from
#' @return A \code{character} vector containing the gene symbol from the SummarizedExperiment
#'
#' @importMethodsFrom SummarizedExperiment rowData
#' @include generics.R
#' @export
setMethod("getGeneSymbols",
          signature('SummarizedExperiment'),
          function(dataset) {
      return(rowData(dataset)$gene)
})