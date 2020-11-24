#' Perform Validation on an `S4` Object Respresenting a Trained Model
#'
#' @param object An `S4` object
#'
#' @md
#' @export
setGeneric('validateModel', function(object, validationData, ...)
    standardGeneric('validateModel'))
#'
#' Evaluate the Performance of a List of Trained KTSP Models from a PCOSP
#'   Object
#'
#' @param object A `PCOSP` model which has been trained using `trainModel`.
#' @param validationData A `CohortList` containing one or more
#'   `SurvivalExperiment`s. The first assay in each `SurvivalExperiment` will
#'   be classified using all top scoring KTSP models in `models(object)`.
#' @param ... Fallthrough arguments to `BiocParallel::bplapply`, use this to
#'   configure the parallelization settings for this function. For example
#'   to specify BPARAM.
#'
#' @seealso BiocParallel::bplapply switchBox::SWAP.KTSP.Classify
#'
#' @return
#'
#' @importFrom BiocParallel bplapply
#' @importFrom switchBox SWAP.KTSP.Classify
#' @md
#' @export
setMethod('validateModel', signature(object='PCOSP',
    validationData='CohortList'), function(object, validationData, ...)
{
    # Extract the best KTSP models
    modelList <- models(object)
    if (length(modelList) < 1)
        stop(.errorMsg(.context(), 'There are no models in the PCOSP model ',
            'passed as the object argument. Please ensure you train your model',
            ' before attempting to validate it.'))

    # Extract the assay matrices from the validation data
    testMatrixList <- lapply(validationData, assays)

    # Get the first assay and warn if there are more than one assays
    for (i in seq_along(testMatrixList)) {
        if (length(testMatrixList[[i]]) > 1)
            warning(.warnMsg(.context(), 'Please ensure your validation ',
                'data only has one assay per SurvivalExperiment! Ignoring ',
                'all but the first assay!'))

        testMatrixList[[i]] <- testMatrixList[[i]][[1]]
    }

    # Classify the validation data using the models
    .classifyWithModels <- function(modelList, inputMat,
        FUN=SWAP.KTSP.Classify)
    {
        predictions <- lapply(modelList, FUN=FUN, inputMat=inputMat)
        predictions <- lapply()
    }

    predictions <- lapply(testMatrixList, FUN=.classifyWithModels,
        modelList=modelList)


})