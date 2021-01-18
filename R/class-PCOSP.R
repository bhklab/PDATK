#' Pancreatic Cancer Overall Survival Predictor (PCOSP) Class
#'
#' @slot trainCohorts A `CohortList` of `SurvivalExperiment`s for training the
#'   PCOSP model.
#' @slot model `` The trained PCOSP model.
#' @slot .intern An `environment` storing internal object metadata. For
#'   developer use.
#'
#' @md
#' @include class-SurvivalExperiment.R
#' @export
.PCOSP <- setClass("PCOSP", contains='SurvivalExperiment',
    slots=list(models='SimpleList', validationStats='data.frame',
        validationData='CohortList'))

#' Pancreatic Cancer Overall Survival Predictor (PCOSP) Constructor
#'
#' @details This function assumes there is only 1 assay per `SurvivalExperiment`.
#'
#' @param trainCohorts A `CohortList`s for training
#'   the PCOSP model, with.
#' @param minDaysSurvived An `integer` indicating the minimum number of day
#'   required to be in the 'high' survival group. Any patients below this
#'   cut-off will be considered low survival. Default is 365 days.
#' @param ... Force subsequent parameters to be named. This parameter is not
#'   used.
#' @param randomSeed An `integer` random seed to use for training the model. If
#'   missing, this defaults to the seed 1234.
#'
#' @return A `PCOSP` object with training data in the assays slot, concatenating
#'   together the molecular data types and labelling the genes with the data
#'   type to ensure the results are easily interpretable.
#'
#' @md
#' @export
PCOSP <- function(trainCohort, minDaysSurvived=365, ..., randomSeed) {

    if (!is(trainCohort, 'SurvivalExperiment'))
        stop(.errorMsg(.context(), 'The trainCohorts argument is not a ',
            'SurvivalExperiment object. Please convert it before building a ',
            'PCOSP model!'))

    assaysL <- assays(trainCohort)

    modelMatrix <- do.call(rbind, as.list(assaysL))
    colData(trainCohort)$prognosis <-  # split into high and low survival
        ifelse(colData(trainCohort)$days_survived >= minDaysSurvived,
            'good', 'bad')

    rowData <- rbind(rowData(trainCohort), rowData(trainCohort))
    rownames(rowData) <- rownames(modelMatrix)

    survExp <- SurvivalExperiment(colData=colData(trainCohort),
        rowData=rowData, assays=SimpleList(trainMatrix=modelMatrix),
        metadata=metadata(trainCohort))

    PCOSPmodel <- .PCOSP(survExp)
    metadata(PCOSPmodel)$modelParams <-
        list(randomSeed=if (!missing(randomSeed)) randomSeed else 1234,
            RNGkind=RNGkind())
    return(PCOSPmodel)
}

##'
##' @importFrom settings is_setting clone_and_merge
##' @export
#setMethod('manageOptions', 'PCOSP', function(where=NULL, ...) {
#  if (settings::is_setting(...)) {
#    where@options <- settings::clone_and_merge(where@options(), ...)
#    where
#  } else {
#    where@options(...)
#  }
#})

#' @noRd
setValidity('PCOSP', function(object) {
    hasSurvivalGroup <- 'prognosis' %in% colnames(colData(object))
    if (hasSurvivalGroup) hasSurvivalGroup else .errorMsg(.context(), 'The ',
        '`survival_group` column is missing. Please use the PCOSP constructor',
        ' to initialize your model and ensure the minDaysSurvived parmater ',
        'is a valid integer!')
})


#' Accessor for the models slot of an `S4` object
#'
#' @param object An `S4` object to retrieve models from.
#' @param ... Allow
#'
#' @export
setGeneric('models', function(object, ...) standardGeneric('models'))
#'
#' Getter for the models slot of a `PCOSP` object
#'
#' @param object A `PCOSP` model object to retrieve the models slot from.
#'
#' @return A `SimpleList` of top scoring KTSPmodels
#'
#' @export
setMethod('models', signature('PCOSP'), function(object) {
    object@models
})

#' Generic for Setting the Models Slot for an S4 Object
#'
#' @param object An `S4` object to set the models slot for
#' @param ... Allow new parameters to be added
#' @param value A model or list of models to assign to the object
#'
#' @return None, updates the object.
#'
#' @export
setGeneric('models<-',
    function(object, ..., value) standardGeneric('models<-'))
#' Setter for the models slot of a `PCOSP` object
#'
#' @param object A `PCOSP` object to update
#' @param value A `SimpleList` of trained KTSP models
#'
#' @return None, updates the object.
#'
#' @export
setReplaceMethod('models', signature(object='PCOSP', value='SimpleList'),
    function(object, value)
{
    object@models <- value
    return(object)
})

#' Accessor for the `validationStats` slot of an `S4` object
#'
#' @param object An `S4` object
#' @param ... Allow definition of new arguments to this generic.
#'
#' @return A `data.frame` of validation statistics for the validation cohorts
#'   provided to `validateModel` function for a given `PCOSP` object.
#'
#' @export
setGeneric('validationStats', function(object, ...)
    standardGeneric('validationStats'))
#'
#' @param object A `PCOSP` object.
#'
#' @export
setMethod('validationStats', signature(object='PCOSP'), function(object) {
    object@validationStats
})

#' Generic for setting the `validationStats` slot on an `S4` object
#'
#' @param object An `S4` object.
#' @param ... Allow definition of additional parameters to this generic.
#' @param value A `data.frame` of validation statistics.
#'
#' @return None, updates the object
#'
#' @export
setGeneric('validationStats<-', function(object, ..., value)
    standardGeneric('validationStats<-'))
#'
#' @param object A `PCOSP` model.
#' @param value A `data.frame` of
#'
#' @export
setReplaceMethod('validationStats', signature(object='PCOSP', value='data.frame'),
    function(object, value)
{
    object@validationStats <- value
    return(object)
})

#' Accessor for the `validationData` slot of an `S4` object
#'
#' @param object An `S4` object
#' @param ... Allow definition of new arguments to this generic.
#'
#' @return A `CohortList` with the validation data used for the `PCOSP` model,
#'   or nothing if the model has not be validated.
#'
#' @export
setGeneric('validationData', function(object, ...)
    standardGeneric('validationData'))
#'
#' @param object A `PCOSP` object.
#'
#' @export
setMethod('validationData', signature(object='PCOSP'), function(object) {
    object@validationData
})

#' Generic for setting the `validationData` slot on an `S4` object
#'
#' @param object An `S4` object.
#' @param ... Allow definition of additional parameters to this generic.
#' @param value A `data.frame` of validation statistics.
#'
#' @return None, updates the object
#'
#' @export
setGeneric('validationData<-', function(object, ..., value)
    standardGeneric('validationData<-'))
#'
#' @param object A `PCOSP` model.
#' @param value A `CohortList` of validation cohorts for the `PCOSP` model.
#'
#'
#' @export
setReplaceMethod('validationData', signature(object='PCOSP', value='CohortList'),
    function(object, value)
{
    object@validationData <- value
    return(object)
})