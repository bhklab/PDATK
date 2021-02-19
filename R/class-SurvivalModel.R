#' A Generic Container for Storing Mathematical Models of
#'   SurvivalExperiments
#'
#' @description An S4 class with a number of predefined methods for
#'   accessing slots relavant to a surival model. More specific model types
#'   will inherit from this class for their accessor methods and constructor.
#'
#' @slot models A `SimpleList` containing one or more model object.
#' @slot validationData A `CohortList` containing one or more `SurvivalExperiment`
#'   objects used to validate the model. This slot is populated by the when the
#'   `validateModel` method is called on a `SurvivalModel` object.
#' @slot validationStats A `data.frame` object containing validation statistics
#'   calculated by the `validateModel` method.
#'
#' @examples
#' data(sampleICGCmicro)
#' set.seed(1987)
#' survModel <- SurvivalModel(sampleICGCmicro, minDaysSurvived=385,
#'   randomSeed=1987)
#'
#' @md
#' @include class-SurvivalExperiment.R
#' @include class-CohortList.R
#' @export
.SurvivalModel <- setClass("SurvivalModel",
    contains='SurvivalExperiment', slots=list(models='SimpleList',
    validationData='CohortList', validationStats='data.frame'
    ))

#' Constructor for a SurvivalModel Object.
#'
#' @param trainCohorts A 'SurvivelExperiment' containing training data for
#'   the `SurvivalModel` object.
#' @param minDaysSurvived An `integer` minimum number of days survived to be
#'   classified as a 'good' prognosis.
#' @param ... Force subsequent paramters to be named. Not used.
#' @param randomSeed An `integer` randomSeed that was used to train the model.
#'   Users should specify this when initializing a model to ensure
#'   reproducibilty.
#'
#' @return A `SurvivalModel` object.
#'
#' @examples
#' data(sampleICGCmicro)
#' set.seed(1987)
#' survModel <- SurvivalModel(sampleICGCmicro, minDaysSurvived=365,
#'   randomSeed=1987)
#'
#' @md
#' @import BiocGenerics
#' @importFrom CoreGx .errorMsg .warnMsg
#' @import S4Vectors
#' @export
SurvivalModel <- function(trainCohorts, minDaysSurvived=365, ...,
    randomSeed)
{
    funContext <- .context(1)

    if (missing(randomSeed)) stop(.errorMsg(funContext, 'No random seed was ',
        'specied for your model. Please include the value used for set.seed ',
        'when training this model! This ensures other can reproduce your ',
        'results.'))

    if (!is(trainCohorts, 'SurvivalExperiment')) {
        if (is(trainCohorts, 'CohortList')) {
            message(funContext, 'Merging trainCohorts `CohortList` into a',
                ' `SurvivalExperiment` with shared genes and samples...')
            .mergeWithNames <- function(x, y) merge(x, y,
                cohortNames=names(trainCohorts))
            trainCohorts <- Reduce(.mergeWithNames, trainCohorts)
            trainCohorts <- dropNotCensored(trainCohorts, minDaysSurvived)

        } else {
            stop(.errorMsg(funContext,
                'The trainCohorts argument is not a CohortList or ',
                'SurvivalExperiment object. Please convert it before building',
                ' a SurvivalModel object!'))
        }
    }

    if (!('prognosis' %in% colnames(colData(trainCohorts)))) {
        colData(trainCohorts)$prognosis <-  # split into high and low survival
            ifelse(colData(trainCohorts)$days_survived >= minDaysSurvived,
                'good', 'bad')
    }

    metadata(trainCohorts)$modelParams <-
        list(randomSeed=randomSeed, RNGkind=RNGkind(),
            minDaysSurvived=minDaysSurvived)

    SurvModel <- .SurvivalModel(trainCohorts, models=SimpleList(),
        validationData=CohortList(mDataTypes=""), validationStats=data.frame())
    return(SurvModel)
}

#' Accessor for the models slot of an `S4` object
#'
#' @param object An `S4` object to retrieve models from.
#' @param ... Allow new parameters to be defined for this generic.
#'
#' @return An R object representing a model.
#'
#' @examples
#' data(samplePCOSPmodel)
#' models(samplePCOSPmodel)
#'
#' @md
#' @export
setGeneric('models', function(object, ...) standardGeneric('models'))
#'
#' Getter for the models slot of a `SurvivalModel` object
#'
#' @param object A `SurvivalModel` model object to retrieve the models slot from.
#'
#' @return A `SimpleList` of top scoring KTSPmodels
#'
#' @examples
#' data(samplePCOSPmodel)
#' models(samplePCOSPmodel)
#'
#' @md
#' @export
setMethod('models', signature('SurvivalModel'), function(object) {
    object@models
})


#' Generic for Setting the Models Slot for an S4 Object
#'
#' @param object An `S4` object to set the models slot for
#' @param ... Allow new parameters to be added.
#' @param value A model or list of models to assign to the object
#'
#' @return None, updates the object.
#'
#' @examples
#' data(samplePCOSPmodel)
#' models(samplePCOSPmodel) <- SimpleList(model1=NA)
#'
#' @md
#' @export
setGeneric('models<-',
    function(object, ..., value) standardGeneric('models<-'))
#'
#' Setter for the models slot of a `SurvivalModel` object
#'
#' @param object A `SurvivalModel` object to update
#' @param value A `SimpleList` of trained KTSP models
#'
#' @return None, updates the object.
#'
#' @examples
#' data(samplePCOSPmodel)
#' models(samplePCOSPmodel) <- SimpleList(model1=NA)
#'
#' @md
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
#' @examples
#' data(samplePCOSPmodel)
#' validationStats(samplePCOSPmodel)
#'
#' @md
#' @export
setGeneric('validationStats', function(object, ...)
    standardGeneric('validationStats'))
#'
#' Accessor for the `validationStats` slot of a `SurvivalModel` object.
#'
#' @param object A `SurvivalModel` object to get validation statistics from.
#'
#' @return A `data.table` of validation statistics for the `SurvivalModel`
#'   object.
#'
#' @examples
#' data(samplePCOSPmodel)
#' validationStats(samplePCOSPmodel)
#'
#' @md
#' @export
setMethod('validationStats', signature(object='SurvivalModel'),
    function(object)
{
    object@validationStats
})


#' Setter for the `validationStats` slot on an `S4` object
#'
#' @param object An `S4` object.
#' @param ... Allow definition of additional parameters to this generic.
#' @param value A `data.frame` of validation statistics.
#'
#' @return None, updates the object
#'
#' @examples
#' data(samplePCOSPmodel)
#' validationStats(samplePCOSPmodel) <- data.frame()
#'
#' @md
#' @export
setGeneric('validationStats<-', function(object, ..., value)
    standardGeneric('validationStats<-'))
#'
#' Setter for the validationStats slot of a `SurvivalModel` object with a
#'   `data.frame`
#'
#' @param object A `SurvivalModel` model.
#' @param value A `data.frame` of validation statistics for a `SurvivelModel`
#'   object.
#'
#' @return None, updated the object.
#'
#' @examples
#' data(samplePCOSPmodel)
#' validationStats(samplePCOSPmodel) <- data.frame()
#'
#' @md
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
#' @examples
#' data(samplePCOSPmodel)
#' validationData(samplePCOSPmodel)
#'
#' @md
#' @export
setGeneric('validationData', function(object, ...) standardGeneric('validationData'))
#'
#' Accessor for the `validationData` slot of a `SurvivalModel` object.
#'
#' @param object A `SurvivalModel` object.
#'
#' @return A `CohortList` object containing the datasets used to compute
#'   validation statistics for this model.
#'
#' @examples
#' data(samplePCOSPmodel)
#' validationData(samplePCOSPmodel)
#'
#' @md
#' @export
setMethod('validationData', signature(object='SurvivalModel'), function(object)
{
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
#' @examples
#' data(samplePCOSPmodel)
#' validationData(samplePCOSPmodel) <- validationData(samplePCOSPmodel)
#'
#' @md
#' @export
setGeneric('validationData<-', function(object, ..., value)
    standardGeneric('validationData<-'))
#'
#' Setter for the `validationData` slot of a `SurvivalModel` object with a
#'  `CohortList`.
#'
#' @param object A `SurvivalModel` model.
#' @param value A `CohortList` of validation cohorts for the `SurvivalModel` model.
#'
#' @return None, updates the object.
#'
#' @examples
#' data(samplePCOSPmodel)
#' validationData(samplePCOSPmodel) <- validationData(samplePCOSPmodel)
#'
#' @md
#' @export
setReplaceMethod('validationData', signature(object='SurvivalModel',
    value='CohortList'), function(object, value)
{
    object@validationData <- value
    return(object)
})

#' Generic for retrieving the randomSeed parameter from a `SurvivalModel` object.
#'
#' This should be used to set the seed before model training to ensure
#' reproducible results.
#'
#' @param object An `S4` object to get the seed from.
#'
#' @return An `integer` seed to be used when training the a `SurivvalModel`.
#'
#' @examples
#' data(sampleSurvivalModel)
#' getModelSeed(sampleSurivalModel)
#'
#' @md
#' @export
setGeneric('getModelSeed', function(object) standardGeneric('getModelSeed'))
#'
#' Method for retrieving the random seed used for training a specific survival
#' model object.
#'
#' This should be used to set the seed before model training to ensure
#' reproducible results.
#'
#' @param object A `SurvivalModel` object to get the seed from.
#'
#' @return An `integer` seed to be used when training the a `SurivalModel`.
#'
#' @examples
#' data(sampleSurvivalModel)
#' getModelSeed(sampleSurivalModel)
#'
#' @md
#' @importFrom S4Vectors metadata
#' @export
setMethod('getModelSeed', signature(object='SurvivalModel'),
    function(object) metadata(object)$modelParams$randomSeed)


#' Class Validity Method for the SurvivalModel Class
#'
#' @param object A `SurvivalModel` object to verify class validity of.
#'
#' @md
#' @importFrom CoreGx .errorMsg
setValidity('SurvivalModel', function(object) {
    funContext <- .context(1)

    hasSurvivalGroup <- 'prognosis' %in% colnames(colData(object))
    if (!hasSurvivalGroup) .errorMsg(funContext, 'The ',
        '`prognosis` column is missing. Object inheriting from the ',
        '`SurvivalModel` class must contain this column!')
    hasModelParams <- 'modelParams' %in% names(metadata(object))
    if (!hasModelParams) .errorMsg(funContext, 'The ',
        '`modelParams` item is missing from the from the `SurvivalModel`',
        'metadata slot. Please ensure you are using the constructor to build ',
        'your `SurvivalModel` object!')

    return(hasSurvivalGroup && hasModelParams)
})
