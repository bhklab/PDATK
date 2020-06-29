#' Generic to retrieve gene symbols from an S4 object
#'
#' @param dataset The S4 object to extract gene symbols from
#' @param ... Allow addition of new parameters to this generic function
#' @return A \code{character} vector with
#'
#' @export
setGeneric('getGeneSymbols', function(dataset, ...) standardGeneric('getGeneSymbols'))