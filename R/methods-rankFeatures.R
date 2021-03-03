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
#' 
#' @md
#' @export
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

#' Rank the Features in a `MultiAssayExperiment` Object
#'
#' @param object A `MultiAssayExperiment` to rank the features in.
#' @param FUN A vectorized feature scoring function, such as `var` or `mad`. 
#'   Defaults to `mad` from the `stats` package.
#' @param RANK_FUN A ranking function, such as `rank` or `dense_rank`. Defaults
#'   to `dense_rank` from `dplyr`.
#' @param ... Fall through arguments to `FUN`, such as `na.rm=TRUE`.
#' @param descending Should your rank function be called with `-` before the
#'   values from `FUN`. Defaults to `TRUE`, which should be used if high values
#'   returned from `FUN` are good.
#' @param weights A named `numeric` weighting vector with a weight for each
#'   experiment in the `MultiAssayExperiment` object. Names must match the
#'   `names(experiments(object))`. Passed to `matrixStats::weightedMedian` when
#'   aggregating feature scores per assay. Defaults to the sample size of an
#'   assay relative to the largest sample size when this paramter is missing.
#'
#' @return The `MultiAssayExperiment` with the item `featureRanks` in the object
#'   metadata, which stores a `DataFrame` containing ranks accross all assays for 
#'   each unique feature and the additional columns `feature_score` and `feature_rank`, 
#'   as calculated with `FUN` and `RANK_FUN`, respectively. Information 
#'   about which functions were used for each column can be found in the object 
#'   `mcols` in the `calculated_with` column.
#' 
# @examples ## TODO:: Add example MAE to the package
# data(sampleICGCmicro)
# rankFeatures(sampleICGCmicro, FUN='mads', RANK_FUN='dense_rank')
#'
#' @seealso [`stats::mad`], [`dplyr::dense_rank`], 
#'   [`matrixStats::weightedMedian`]
#' 
#' @importFrom dplyr dense_rank
#' @importFrom stats mad
#' @importFrom matrixStats weightedMedian
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder set
#' 
#' @md
#' @export
setMethod('rankFeatures', signature(object='MultiAssayExperiment'),
    function(object, FUN='mad', RANK_FUN='dense_rank', ..., 
        descending=TRUE, weights)
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

    # convert all assays to long format
    assayList <- lapply(assays(object), as.data.table, keep.rownames='feature')
    meltedAssayList <- lapply(assayList, melt.data.table, id.vars='feature', 
        variable.name='sample', value.name='value')
    for (i in seq_along(meltedAssayList)) {
        set(meltedAssayList[[i]], j='assay', value=names(assayList)[i])
    }
    longDT <- rbindlist(meltedAssayList)

    featureAssayDT <- longDT[, .(feature_score=FUN(value, ...)), 
        by=.(feature, assay)]
    if (missing(weights)) {
        sampleSizes <- vapply(experiments(object), FUN=ncol, numeric(1))
        weights <- sampleSizes / max(sampleSizes)
    }
    assayHasWeights <- names(assayList) %in% names(weights)
    if (!all(assayHasWeights)) {
        stop(.errorMsg(funContext, 'The assays ', 
            paste0(missignAssays, collapse=', '),
            ' do not have weights, in the weights vector. Please specify', 
            ' a weight for all assays in your MultiAssayExperiment or ', 
            'exclude this weights argument to have them automatically ',
            'calculated.'))
    }
    featureDT <- featureAssayDT[, 
        .(feature_score=weightedMedian(feature_score, w=weights[assay], 
            na.rm=TRUE)), 
        by='feature']

    featureDT[, 
        feature_rank := if (descending) RANK_FUN(-feature_score) else 
            RANK_FUN(feature_score)
            ]

    featureDF <- DataFrame(featureDT[order(feature_rank)])
    mcols(featureDF)$calculated_with <- c(NA, FUN_NAME, RANK_FUN_NAME)

    metadata(object)$featureRanks <- featureDF
    return(object)
})