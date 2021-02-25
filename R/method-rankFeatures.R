#' Rank the Features in a `S4` Object
#'
#' @param object An `S4` object where rows represent features.
#' @param ... Allow new parameters to be defined for this generic.
#'
#' @return The `object` with the features per value and ranking of the features
#'   in the `rowData` slow of the object.
#'
#' @examples
#' data(sampleICGCmicro)
#' rankFeatures(sampleICGCmicro)
#'
#' @md
#' @export
setGeneric('rankFeatures', function(object, ...) standardGeneric('rankFeatures'))
#'
#' Rank the Features in a `SummarizedExperiment` Object
#'
#' @param object A `SummarizedExperiment` to rank the features in.
#' @param FUN A row-wise summary function, such as `rowVars`, `rowMads`, etc.
#'   defaults to `MatrixGenerics::rowMads`.
#' @param RANK_FUN A ranking function, such as `rank` or `dense_rank`. Defaults
#'   to `dplyr::dense_rank`.
#' @param ... Fall through arguments to `FUN`, such as `na.rm=TRUE`.
#' @param descending Should your rank function be called with `-` before the
#'   values from `FUN`. Defaults to `TRUE`.
#' @param assay `integer` assay to use for the ranking, as passed to the
#'   `SummarizedExperiment::assay` function. Defaults to the first assay.
#'
#' @return The `SummarizedExperiment` with the column `feature_score` and
#'   `feature_rank` in the `rowData` slot. Information about which functions
#'   where used for each column can be found in the object `mcols` in the
#'   `calculated_with` column.
#'
#' @examples
#' data(sampleICGCmicro)
#' rankFeatures(sampleICGCmicro, FUN='rowMads', RANK_FUN='dense_rank')
#'
#' @importFrom MatrixGenerics rowMads
#' @importFrom dplyr dense_rank
setMethod('rankFeatures', signature(object='SummarizedExperiment'),
    function(object, FUN='rowMads', RANK_FUN='dense_rank', ..., descending=TRUE,
        assay=1)
{
    funContext <- .context(1)

    FUN_NAME <- as.character(substitute(FUN))
    RANK_FUN_NAME <- as.character(substitute(RANK_FUN))

    if (is.character(FUN)) FUN <- get(FUN)
    if (is.character(RANK_FUN)) RANK_FUN <- get(RANK_FUN)
    if (!is.function(FUN)) stop(.errorMsg(funContext, 'The argument to the FUN ',
        'parameter is not a function... Please ensure you pass a function or
        the name of a function to calculate row summaries with!'))
    if (!is.function(RANK_FUN)) stop(.errorMsg(funContext, 'The argument to the',
        'RANK_FUN parameter is not a function... Please ensure you pass a ',
        'a function or the name of a function to calculate the ranks with.'))

    rowData(object)$feature_score <- FUN(assay(object, assay), ...)
    rowData(object)$feature_rank <-
        if (descending) RANK_FUN(-rowData(object)$feature_score) else
            RANK_FUN(rowData(object)$feature_score)

    mcols(rowData(object))$calculated_with <- NA
    mcols(rowData(object))$which_assay <- NA
    mcols(rowData(object))[c('feature_score', 'feature_rank'), ]$calculated_with <-
        c(FUN_NAME, RANK_FUN_NAME)
    mcols(rowData(object))[c('feature_score', 'feature_rank'), ]$which_assay <-
        names(assays(object))[assay]

    return(object)
})