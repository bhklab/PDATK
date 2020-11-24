#' Predict Classes for New Data Based on a Train Classifier Model
#'
#' @param object An `S4` object containing data
#' @param model An `S4` object containing one or more trained classification
#'   model.
#'
#' @md
#' @export
setGeneric('predictClasses',
    function(object, model, ...) standardGeneric('predictClasses'))
#'
#' @param object A `SurvivalExperiment` object to predict classes for.
#' @param model A trained `PCOSP` model to use for predicting classes.
#' @param ... Fall through arguments to `BiocParallel::bplapply` for configuring
#'   parallelization settings.
#'
#' @seealso BiocParallel::bplapply switchBox::SWAP.KTSP.Classify
#'
#' @return A `SurvivalExperiment` with the predictions in its metadata and
#'   a column in colData, `prob_good_survival`, which contains the proportion
#'   of models which predicted good prognosis for each sample.
#'
#' @md
#' @export
setMethod('predictClasses', signature(object='SurvivalExperiment',
    model='PCOSP'), function(object, model, ...)
{
    modelList <- models(model)
    if (length(modelList) < 1)
        stop(.errorMsg(.context(), 'There are no models in the PCOSP model ',
            'passed as the object argument. Please ensure you train your model',
            ' with `trainModel` before attempting to predict classes with it.'))

    assayData <- assays(object)
    if (length(assayData) > 1)
        warning(.warnMsg(.context(), 'Please ensure your prediction ',
            'data only has one assay! Ignoring ',
            'all but the first assay!'))

    assayMatrix <- assayData[[1]]

    predictionList <- bplapply(modelList, SWAP.KTSP.Classify,
        inputMat=assayMatrix)
    # convert factor to character
    predictionListChar <- lapply(predictionList, as.character)
    predictions <- do.call(rbind, predictionListChar)
    colnames(predictions) <- colnames(assayData)

    metadata(object)$PCOSP_predictions <- predictions

    colData(object)$proportion_good_prognosis <-
        colSums(predictions == 'good') / nrow(prediction)

    return(object)
})
#'
#' @param object A `CohortList` with `SurvivalExperiment`s to predict classes
#'   for.
#' @param model A trained `PCOSP` model to use for predicting classes.
#' @param ... Fall through arguments to `BiocParallel::bplapply` for configuring
#'   parallelization settings.
#'
#' @return A `CohortList` with the model predictions attached to each
#'   `SurvivalExperiment` in the metadata slot and the `prob_good_survival`
#'   column added to the colData slot.
#'
#' @seealso BiocParallel::bplapply switchBox::SWAP.KTSP.Classify
#'
#' @md
#' @export
setMethod('predictClasses', signature(object='CohortList',
    model='PCOSP'), function(object, model)
{
    predictionResults <- endoapply(object, predictClasses, model=model, ...)
    return(predictionResults)
})