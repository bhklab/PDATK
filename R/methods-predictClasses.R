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

## TODO:: refactor to use a hook for the prediction method

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

#'
#' @param object A `SurvivalExperiment` object with the correct columns in
#'   `colData` to match the formula for the `ClinicalModel` object.
#' @param model A trained `ClinicalModel` object, as return by `trainModel`.
#' @param ... Fall through parameters to [`stats::predict`].
#' @param na.action The `na.action` paramter passed to [`stats::predict.glm`].
#' @param type The `type` parameter passed to [`stats::predict.glm`]
#'
#' @return A `SurvivalExperiment` with the model predictions in the colData
#'   slot as clinical_prob_good.
#'
#' @md
#' @importFrom stats predict glm
#' @export
setMethod('predictClasses', signature(object='SurvivalExperiment',
    model='ClinicalModel'), function(object, model, ..., na.action='na.exclude',
        type='response')
{

    # check that the formula is valid and the variables are in the training data
    formula <- as.formula(metadata(model)$modelParams$formula)
    formulaCols <- as.character(formula[seq(2, 3)])
    # split the formula into a vector where each variable is an item
    formulaCols <- unlist(strsplit(formulaCols,
        split='[\\s]*[\\+\\-\\~\\=\\*][\\s]*', perl=TRUE))
    hasFormulaCols <- formulaCols %in% colnames(colData(object))
    if (!all(hasFormulaCols))
        stop(.errorMsg(.context(), 'The columns ', formulaCols[!hasFormulaCols],
            ' are missing from the colData slot of the training data',
            'Please only specify valid column names in colData to the formula!'))
    if (length(models(model)) > 1) warning(.warnMsg(.context(), 'There is more
        than one model in your ClinicalModel. Only using the first one...'))

    # Skip rows with levels that aren't in the model; prevents predict.glm for
    #   breaking if there are new levels prediction data
    modelFactorLevels <- models(model)$glm$xlevels
    keepRows <- rep(TRUE, nrow(colData(object)))
    for (name in names(modelFactorLevels)) {
        keep <- colData(object)[[name]] %in% modelFactorLevels[[name]]
        keepRows <- keepRows & keep
    }
    if (!all(keepRows))
        warning(.warnMsg(.context(1), 'Rows ', paste0(which(!keepRows), collapse=', '),
            ' have levels that are not in the model, skipping these rows...'))

    # Calculate survival probabiltiies
    predictions <- predict(models(model)[[1]],colData(object)[keepRows, ], ...,
        na.action=na.action, type=type)

    # Update the `SurvivalExperiment` object with the predicted probabilities
    colData(object)$clinical_prob_good <- NA
    colData(object)$clinical_prob_good[keepRows] <- predictions

    metadata(object)$GLMpredictions <- matrix(
        ifelse(colData(object)$clinical_prob_good > 0.5, 'good', 'bad'),
        byrow=TRUE, nrow=1, dimnames=list('glm', rownames(colData(object))))
    metadata(object)$GLMparams <- metadata(model)$modelParams

    return(object)
})
#'
#' Use a Clinical GLM to Predict Classes for a `CohortList` of
#'   `SurvivalExperment` Objects.
#'
#' @param object A `CohortList` with `SurvivalExperiment`s to predict classes
#'   for. The colData slot in ALL `SurvivalExperiment`s must have column names
#'   which match the formula in the model object.
#' @param model A trained `ClinicalModel` object, as return by `trainModel`.
#' @param ... Fall through parameters to [`stats::predict`].
#' @param na.action The `na.action` paramter passed to [`stats::predict.glm`].
#' @param type The `type` parameter passed to [`stats::predict.glm`]
#'
#' @return A `CohortList` with the model predictions in the colData
#'   slot as clinical_prob_good for each `SurvivalExperiment`, and the
#'   model in the metadata as predictionModel.
#'
#' @md
#' @export
setMethod('predictClasses', signature(object='CohortList',
    model='ClinicalModel'), function(object, model, ..., na.action='na.exclude',
        type='response')
{
    predictionResults <- endoapply(object, predictClasses, model=model, ...,
        na.action=na.action, type=type)
    mcols(predictionResults)$hasPredictions <- TRUE
    metadata(predictionResults)$predictionModel <- model
    return(predictionResults)
})

#'
#'
#'
#' @importFrom genefu sig.score
#' @export
setMethod('predictClasses', signature(object='SurvivalExperiment',
    model='GeneFuModel'), function(object, model, ..., annot=NA)
{

    # Calculate survival releative risks (not probabiltiies)
    predictions <- lapply(models(model), FUN=genefu::sig.score,
        data=t(assay(object, 1)), annot=annot, ...)

    # Add a column to colData for each model in the GeneFuModel
    for (name in names(models(model))) {
        colData(object)[[paste('genefu', name, 'score', sep='_')]] <-
            predictions[[name]]$score
    }

    .ifElseScoreLt <- function(item, cutoff) ifelse(item$score < cutoff,
        'good', 'bad')
    metadata(object)$genefuPredictions <- lapply(predictions, .ifElseScoreLt,
        cutoff=2)
    metadata(object)$genefuParams <- lapply(predictions, `[`, -1)

    return(object)
})


#' Use a Gene Signature Based Prediciton Model from the `genefu` Package
#'   to Predict Signature Scores
#'
#' @param object A `CohortList` with `SurvivalExperiment`s to predict classes
#'   for.
#' @param model A trained `GeneFuModel` object.
#' @param ... Fall through parameters to [`genefu::sig.score`].
#' @param  annot The `annot` parameter passed to [`genefu::sig.score`].
#'   Defaults to NA, which assumes your assay rowname match the gene labels
#'   in the model.
#'
#' @return A `CohortList` with the model predictions in the colData
#'   slot as genefu_<model_name>_score for each `SurvivalExperiment`, and the
#'   model in the metadata as predictionModel.
#'
#' @md
#' @export
setMethod('predictClasses', signature(object='CohortList',
    model='GeneFuModel'), function(object, model, ..., annot=NA)
{
    predictionResults <- endoapply(object, predictClasses, model=model, ...,
        annot=annot)
    mcols(predictionResults)$hasPredictions <- TRUE
    metadata(predictionResults)$predictionModel <- model
    return(predictionResults)
})