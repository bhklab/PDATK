#' Get the Top Predictive Features from an S4 Object
#'
#' @param object An `S4` object to get top scoring features from.
#' @param ... Allow additional parameters to be defined for this generic.
#'
#' @return A `character` vector of top predictive features.
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