#' @export
setClassUnion('PCOSP_or_RLS_or_RGA', c('PCOSP', 'RLSModel', 'RGAModel'))

#' Predict Classes for New Data Based on a Train Classifier Model
#'
#' @param object An `S4` object containing data to predict classes from.
#' @param model An `S4` object containing one or more trained classification
#'   models.
#' @param ... Allow further parameters to be defined on this generic.
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
#' @seealso [`BiocParallel::bplapply`], [`switchBox::SWAP.KTSP.Classify`]
#'
#' @return A `SurvivalExperiment` with the predictions in its metadata and
#'   a column in colData, `prob_good_survival`, which contains the proportion
#'   of models which predicted good prognosis for each sample.
#'
#' @md
#' @include
#' @export
setMethod('predictClasses', signature(object='SurvivalExperiment',
    model='PCOSP_or_RLS_or_RGA'), function(object, model, ...)
{
    # drop NA samples, they mess with calculating statistics
    keepSamples <- rownames(na.omit(colData(object)))
    if (!all(colnames(object) %in% keepSamples)) {
        warning(.warnMsg(.context(), 'Dropped sampels with NA survival data!'))
    }
    object <- object[, keepSamples]

    modelList <- models(model)
    if (length(modelList) < 1)
        stop(.errorMsg(.context(), 'There are no models in the PCOSP model ',
            'passed as the model argument. Please ensure you train your model',
            ' with `trainModel` before attempting to predict classes with it.'))

    assayData <- assays(object)
    if (length(assayData) > 1)
        warning(.warnMsg(.context(), 'Please ensure your prediction ',
            'data only has one assay! Ignoring ',
            'all but the first assay!'))

    assayMatrix <- assayData[[1]]

    predictionList <- bplapply(modelList, FUN=SWAP.KTSP.Classify,
        inputMat=assayMatrix)
    # convert factor to character
    predictionListChar <- lapply(predictionList, as.character)
    predictions <- do.call(rbind, predictionListChar)
    colnames(predictions) <- colnames(assayMatrix)

    metadata(object)$PCOSPpredictions <- predictions
    metadata(object)$PCOSPparams <- metadata(model)$modelParams

    colData(object)$PCOSP_prob_good <-
        colSums(predictions == 'good') / nrow(predictions)
    colData(object)$prognosis <-
        ifelse(colData(object)$days_survived > 365,'good', 'bad')

    return(object)
})



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
#' @importFrom S4Vectors endoapply
#' @md
#' @export
setMethod('predictClasses', signature(object='CohortList',
    model='PCOSP_or_RLS_or_RGA'), function(object, model, ...)
{
    predictionResults <- endoapply(object, predictClasses, model=model, ...)
    mcols(predictionResults)$hasPredictions <- TRUE
    metadata(predictionResults)$predictionModel <- model
    return(predictionResults)
})
