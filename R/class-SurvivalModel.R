#' A Virtual Generic Container for Storing Mathematical Models of
#'   SurvivalExperiments
#'
#' @description A virtual S4 class with a number of predefined methods for
#'   accessing slots relavant to a surival model. More specific model types
#'   will inherit from this class for their accessor methods.
#'
#' @slot colData A `DataFrame` with merge column-level metadata for the model
#'   training data in the `assays` slot.
#' @slot rowData A `DataFrame` with merge row-level metadata for the model
#'   training data in the `assays` slot.
#' @slot assays A `SimpleList` of matrices, one for each cohort in the training
#'   data for a `SurivalModel`.
#' @slot models A `SimpleList` containing one or more model object.
#' @slot validationData A `CohortList` containing one or more `SurvivalExperiment`
#'   objects used to validate the model. This slot is populated by the when the
#'   `validateModel` method is called on a `SurvivalModel` object.
#' @slot validationStats A `data.frame` object containing validation statistics
#'   calculated by the `validateModel` method.
#' @slot modelParams A `list` containing configurable parameters related to
#'   training the `SurvivalModel`
#' @slot metadata A `SimpleList` of metadata related to the PCOSP model. This
#'   slot will contain a modelParams item, which has the relevant parameters
#'   used when training the model.
#'
#' @describeIn SurivalExperiment.Rd
#' @md
#' @include class-SurvivalExperiment.R
#' @export
setClass("SurvivalModel", contains=c('SurvivalExperiment', 'VIRTUAL'),
    slots=list(models='SimpleList', validationData='CohortList',
        validationStats='data.frame'
        ))

#' Accessor for the models slot of an `S4` object
#'
#' @param object An `S4` object to retrieve models from.
#' @param ... Allow
#'
#' @export
setGeneric('models', function(object, ...) standardGeneric('models'))
#'
#' Getter for the models slot of a `SurvivalModel` object
#'
#' @param object A `SurvivalModel` model object to retrieve the models slot from.
#'
#' @return A `SimpleList` of top scoring KTSPmodels
#'
#' @export
setMethod('models', signature('SurvivalModel'), function(object) {
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
#' Setter for the models slot of a `SurvivalModel` object
#'
#' @param object A `SurvivalMdeol` object to update
#' @param value A `SimpleList` of trained KTSP models
#'
#' @return None, updates the object.
#'
#' @export
setReplaceMethod('models', signature(object='SurvivalModel',
    value='SimpleList'), function(object, value)
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
#'   provided to `validateModel` function for a given `SurvivalModel` object.
#'
#' @export
setGeneric('validationStats', function(object, ...)
    standardGeneric('validationStats'))
#'
#' @param object A `` object.
#'
#' @export
setMethod('validationStats', signature(object='SurvivalModel'), function(object) {
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
#' @param object A `SurvivalModel` model.
#' @param value A `data.frame` of
#'
#' @export
setReplaceMethod('validationStats', signature(object='SurvivalModel',
    value='data.frame'), function(object, value)
{
    object@validationStats <- value
    return(object)
})

#' Accessor for the `validationData` slot of an `S4` object
#'
#' @param object An `S4` object
#' @param ... Allow definition of new arguments to this generic.
#'
#' @return A `CohortList` with the validation data used for the `SurvivalModel` model,
#'   or nothing if the model has not be validated.
#'
#' @export
setGeneric('validationData', function(object, ...)
    standardGeneric('validationData'))
#'
#' @param object A `SurvivalModel` object.
#'
#' @export
setMethod('validationData', signature(object='SurvivalModel'), function(object) {
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
#' @param object A `SurvivalModel` model.
#' @param value A `CohortList` of validation cohorts for the `SurvivalModel` model.
#'
#'
#' @export
setReplaceMethod('validationData', signature(object='SurvivalModel',
    value='CohortList'), function(object, value)
{
    object@validationData <- value
    return(object)
})

#' @noRd
setValidity('SurvivalModel', function(object) {
    hasSurvivalGroup <- 'prognosis' %in% colnames(colData(object))
    if (hasSurvivalGroup) hasSurvivalGroup else .errorMsg(.context(), 'The ',
        '`prognosis` column is missing. Object inheriting from the ',
        '`SurvivalModel` virtual class must contain this column!')
    hasModelParams <- 'modelParams' %in% colnames(metadata(object))
    if (hasModelParams) hasModelParams else .errorMsg(.context(), 'The ',
        '`modelParams` item is missing from the from the `SurvivalModel`',
        'metadata slot. Please ensure you are using the constructor to build ',
        'your `SurvivalModel` object!')
})