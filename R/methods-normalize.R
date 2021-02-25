#' Normalize the `assays` in a `SummarizedExperiment` Object
#'
#' @param object A `SummarizedExperiment` object with assays to normalize.
#' @param MARGIN An `integer` indicating if rows (1) or columns (2) should be normalized. Defaults to
#'   columns. Defaults to 2.
#' @param FUN A function to normalize your data with. The function should be
#'   vectorized and 
#' @param ... Fall through parameters to FUN.
#' @param whichAssays A `numeric` or `character` vector specifying the indices
#'   of the assays to normalize.
#' 
#' @return The `SummarizedExperiment` with one or more of the matrices in `assays` normalized.
#' 
#' @importMethodsFrom BiocGenerics normalize
#' 
#' @md
#' @export
setMethod('normalize', signature(object='SummarizedExperiment'), 
    function(object, MARGIN=2, FUN='scale', ..., whichAssays=seq_len(assays(object))) 
{
    funContext <- .context(1)

    FUN_NAME <- as.character(substitute(FUN))

    if (is.character(FUN)) FUN <- get(FUN)
    if (!is.function(FUN)) stop(.errorMsg(funContext, 'The argument to the FUN ',
        'parameter is not a function... Please ensure you pass a function or
        the name of a function to calculate row summaries with!'))
    
    if (MARGIN == 1) FUN <- function(x, ...) t(FUN(t(x), ...))

    assays(object)[whichAssays] <- endoapply(assays(object)[whichAssays], FUN, ...)
    return(object)
})

#'
#' 
#' 
#' 
#' 
#' 
