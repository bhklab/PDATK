#' Get the Top Predictive Features from an S4 Object
#'
#' @param object An `S4` object to get top scoring features from.
#' @param ... Allow additional parameters to be defined for this generic.
#'
#' @return A `character` vector of top predictive features.
#'
#' @examples
#' data(sampleTrainedPCOSPmodel)
#'
#' # Get the top features
#' topFeatures <- getTopFeatures(sampleTrainedPCOSPmodel, numModels=2)
#'
#' @md
#' @export
setGeneric('getTopFeatures',
    function(object, ...) standardGeneric('getTopFeatures'))
#'
#' Get the top features for classification using a PCOSP model object.
#'
#' @param object A `PCOSP` model object which has been trained with `trainModel`.
#' @param numModels An `integer` specifying the number of top models to
#'   use features from. Defaults to top 10% of KTSPs models.
#'
#' @return A `character` vector of gene names representing the unique genes
#'   from the top `numModels` KTSPs models in the model `object`.
#'
#' @examples
#' data(sampleTrainedPCOSPmodel)
#'
#' # Get the top features
#' topFeatures <- getTopFeatures(sampleTrainedPCOSPmodel, numModels=5)
#'
#' @md
#' @export
setMethod('getTopFeatures', signature(object='PCOSP'),
    function(object, numModels)
{
    if (missing(numModels)) numModels <- floor(length(models(object))*0.1)

    TSPs <- lapply(models(object), `[[`, 'TSPs')
    topFeatures <- unique(unlist(TSPs))
    return(topFeatures)
})

#' Get the Top Ranked Features in a `SummarizedExperiment` object
#'
#' @param object A `SummarizedExperiment` to extract top features from
#' @param numFeats An `integer` number of top ranked features to extract.
#' @param ... Fall through arguments to `rankFeatures`, which is used to
#'   calculate the ranks if `rankCol` is not present the object `rowData`.
#' @param rankCol The name of the column containing the integer feature ranks.
#'   Defaults to `feature_rank`, as calculated with `rankFeatures`, but users
#'   can alternatively specify their own custom rank column if desired.
#'
#' @return A `character` vector of top ranked features, with the features in
#'   order of rank ascending.
#'
#' @examples
#' data(sampleICGCmicro)
#' getTopFeatures(sampleICGCmicro)
#'
#' @md
#' @export
setMethod('getTopFeatures', signature(object='SummarizedExperiment'),
    function(object, numFeats, ..., rankCol='feature_rank')
{
    funContext <- .context(1)
    if (!(rankCol %in% colnames(rowData(object)))) {
        if (rankCol != ' feature_rank') {
            warning(.warnMsg(funContext, 'The ',
                'column ', rankCol, ' is missing from the data. Calculating new ',
                'ranks with `rankFeature` and using the "feature_rank" column ',
                'to get the top features'))
            rankCol <- 'feature_rank'
        }
        object <- rankFeatures(object, ...)
    }
    if (!is.numeric(numFeats)) stop(.errorMsg(funContext, 'The numFeats ',
        'parameter is not numeric!'))

    sortedRowData <- rowData(object)[order(rowData(object)[[rankCol]]), ]
    if (is.null(rownames(sortedRowData))) stop(.errorMsg('This function ',
        'requires rownames to work, please assign your feature names to the ',
        'rownames of rowData!'))
    topFeatures <- rownames(sortedRowData)[seq_len(numFeats)]

    return(topFeatures)
})