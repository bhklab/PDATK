setClassUnion('List_or_list_or_NULL', c('List', 'list', NULL))
#' @importClassesFrom data.table data.table
setClassUnion('DFrame_or_data.frame_data.table_or_NULL',
    c('DFrame', 'data.frame', 'data.table', NULL))

#' An `S4` Virtual Class For the Concept of a Statistical or ML Model
#'
#' @slot trainData An object inheriting from `List` or `list` representing
#'   the training data for the model.
#' @slot modelParams An object inherting from `List` or `list` representing
#'   the parameters needed to train the model.
#' @slot models An object inheriting from `List` or `list` representing the
#'   trained models.
#' @slot validationStats An object inheriting from `DFrame` or `data.frame`
#'   and storing statistics assessing model performance.
#' @slot validationData An object inheriting `List` or `list` representing
#'   the data used to validate or evaluate the performance of a model.
#' @slot elementMetadata A `DataFrame` or 'data.frame' of item metadata for
#'   the `models` slot.
#' @slot metadata A `List` or `list` of model level metadata.
#'
#' @md
#' @aliases S4Model-class
#' @export
.S4Model <- setClass('S4Model',
    slots=c(
        trainData='ANY',
        modelParams='List_or_list_or_NULL',
        models='List_or_list_or_NULL',
        validationStats='DFrame_or_data.frame_data.table_or_NULL',
        validationData='List_or_list_or_NULL'
        ),
    contains=c('Vector', 'VIRTUAL'))

#' Generic for Accessing the Training Data of an `S4` Object
#'
#' @param object An `S4` object to retrieve training data from.
#' @param ... Allow new parameters to be defined for this generic.
#'
#' @return The training data for an S4 object.
#'
#' @export
setGeneric('trainData', function(object, ...) standardGeneric('trainData'))
#'
#' Accessor for the Training Data in a `S4Model` Object
#'
#' @param object An `S4Model` object to retrieve training data from.
#'
#' @return The training data for an `S4Model` Object.
#'
#' @md
#' @export
setMethod('trainData', signature(object='S4Model'), function(object)
{
    return(object@trainData)
})

#' Generic for Accessing the Training Data of an `S4` Object
#'
#' @param object An `S4` object to retrieve training data from.
#' @param ... Allow new parameters to be defined for this generic.
#' @param value An object to place in the objects training data slot.
#'
#' @return None, updates the object.
#'
#' @md
#' @export
setGeneric('trainData<-', function(object, ..., value)
    standardGeneric('trainData<-'))
#'
#' Accessor for the Training Data in a `S4Model` Object
#'
#' @param object An `S4Model` object to retrieve training data from.
#' @param value An object to put into the model training data.
#'
#' @return The training data for an `S4Model` Object.
#'
#' @md
#' @keywords internal
#' @export
setReplaceMethod('trainData', signature(object='S4Model'),
    function(object, value)
{
    funContext <- .context(1)
    if (sys.nframe() == 1) warning(.warnMsg(funContext, 'The training data of
        a ', class(object), ' object should not be modified manually, please
        use the constructor instead.'))

    object@trainData <- value
    return(object)
})

#' Generic for Accessing the Model Parameters of an `S4` Object
#'
#' @param object A `S4` Object.
#' @param ... Allow additional arguments to be defined for this generic.
#'
#' @return A `List`- or `list`-like object containing all the parameters needed
#'   to reproduce a specific model training run, including environmental
#'   settings like the random seed and RNGkind.
#'
#'
#' @md
#' @export
setGeneric('modelParams', function(object, ...) standardGeneric('modelParams'))
#'
#' Accessor for the Model Parameters of an `S4Model` Object
#'
#' @param object A `S4Model` Object
#'
#' @return A `List`- or `list`-like object containing all the parameters needed
#'   to reproduce a specific model training run, including environmental
#'   settings like the random seed and RNGkind.
#'
#' @md
#' @export
setMethod('modelParams', signature(object='S4Model'), function(object)
{
    return(object@modelParams)
})

#' Generic for Setting the Model Parameters of An `S4` Object
#'
#' @param object An `S4` Object
#' @param ... Allow additional parameters to be defined for this generic.
#' @param value A `List`- or `list`-like object containing the parameters
#'   for the model.
#'
#' @return None, modifies the object.
#'
#' @md
#' @keywords internal
#' @export
setGeneric('modelParams<-', function(object, ..., value)
    standardGeneric('modelParams<-'))
#'
#' Setter for the `modelParams` of an `S4Model`
#'
#' @param object An `S4Model`
#' @param value A `List`- or `list`-like object containing the parameters
#'   for the model.
#'
#' @details
#' This method if not intended for interactive use and will throw an warning
#'   when used interactively.
#'
#' @return None, modifies the object.
#'
#' @md
#' @keywords internal
#' @export
setReplaceMethod('modelParams', signature(object='S4Model',
    value='List_or_list_or_NULL'), function(object, value)
{
    funContext <- .context(1)
    if (sys.nframe() == 1) warning(.warnMsg(funContext, 'The model parameters of
        a ', class(object), ' object should not be modified manually, please
        use the constructor instead.'))

    object@modelParams <- value
    return(object)
})

#' Accessor for the models slot of an `S4` object
#'
#' @param object An `S4` object to retrieve models from.
#' @param ... Allow new parameters to be defined for this generic.
#'
#' @return An `S4` object representing a model.
#'
#' @md
#' @export
setGeneric('models', function(object, ...) standardGeneric('models'))
#'
#' Getter for the models slot of an `S4Model` Object
#'
#' @param object An `S4Model` object to retrieve models from.
#'
#' @return An `List`- or `list`-like object representing a model.
#'
#' @md
#' @export
setMethod('models', signature(object='S4Model'), function(object)
{
    return(object@models)
})

#' Accessor for the models slot of an `S4` object
#'
#' @param object An `S4` object to retrieve models from.
#' @param ... Allow new parameters to be defined for this generic.
#' @param value A `List`- or `list`-like object.
#'
#' @return None, updates the object.
#'
#' @md
#' @keywords internal
#' @export
setGeneric('models<-', function(object, ..., value) standardGeneric('models<-'))
#'
#' Setter for the Models Slot of an `S4Model` Object
#'
#' @param object An `S4Model` object to retrieve models from.
#' @param value A `List`- or `list`-like object.
#'
#' @return An `List`- or `list`-like object representing a model.
#'
#' @md
#' @keywords internal
#' @export
setReplaceMethod('models', signature(object='S4Model',
    value='List_or_list_or_NULL'), function(object, value)
{
    funContext <- .context(1)
    if (sys.nframe() == 1) warning(.warnMsg(funContext, 'The model parameters of
        a ', class(object), ' object should not be modified manually, please
        use the constructor instead.'))

    object@models <- value
    return(object)
})

#' Accessor for the `validationStats` slot of an `S4` object
#'
#' @param object An `S4` object
#' @param ... Allow definition of new arguments to this generic.
#'
#' @return A `data.frame` of validation statistics for the validation data
#'   provided to `validateModel` function for a given `S4` object.
#'
#' @md
#' @export
setGeneric('validationStats', function(object, ...)
    standardGeneric('validationStats'))
#'
#' Acessor for the `validationStats` slot of an `S4Model` Object
#'
#' @param object An `S4` object
#'
#' @return A `data.frame` of validation statistics for the validation data
#'   provided to `validateModel` function for a given `S4Model` object.
#'
#' @md
#' @export
setMethod('validationStats', signature(object='S4Model'),
    function(object)
{
    return(object@validationStats)
})

#' Setter for the `validationStats` slot on an `S4` object
#'
#' @param object An `S4` object.
#' @param ... Allow definition of additional parameters to this generic.
#' @param value A `DataFrame`- or `data.frame`-like object of validation
#'   statistics.
#'
#' @return None, updates the object
#'
#' @md
#' @export
setGeneric('validationStats<-', function(object, ..., value)
    standardGeneric('validationStats<-'))
#'
#' Setter for the `validationStats` slot on an `S4Model` object
#'
#' @param object An `S4Model` object.
#' @param value A `DataFrame`- or `data.frame`-like object of validation
#'   statistics.
#'
#' @return None, updates the object
#'
#' @md
#' @export
setReplaceMethod('validationStats', signature(object='S4Model',
    value='DFrame_or_data.frame_data.table_or_NULL'), function(object, value)
{
    object@validationStats <- value
    return(object)
})


#' Accessor for the `validationData` slot of an `S4` object
#'
#' @param object An `S4` object
#' @param ... Allow definition of new arguments to this generic.
#'
#' @return A `List`- or `list`-like object containing one or more
#'   sets of validation data.
#'
#' @md
#' @export
setGeneric('validationData', function(object, ...)
    standardGeneric('validationData'))
#'
#' Acessor for the `validationData` slot of an `S4Model` object
#'
#' @param object An `S4Model` object
#'
#' @return A `List`- or `list`-like object containing one or more
#'   sets of validation data.
#'
#' @md
#' @export
setMethod('validationData', signature(object='S4Model'),
    function(object)
{
    return(object@validationData)
})


#' Generic for setting the `validationData` slot on an `S4` object
#'
#' @param object An `S4` object.
#' @param ... Allow definition of additional parameters to this generic.
#' @param value A `data.frame` of validation statistics.
#'
#' @return None, updates the object
#'
#' @md
#' @export
setGeneric('validationData<-', function(object, ..., value)
    standardGeneric('validationData<-'))
#'
#' Setter Method for the `validationData` of an `S4Model` Object.
#'
#' @param object An `S4Model` object.
#' @param value A `List`- or `list`-like object containing
#'   the new validation data for the `S4Model` object.
#'
#' @return None, updates the object.
#'
#' @md
#' @keywords internal
#' @export
setReplaceMethod('validationData', signature(object='S4Model',
    value='List_or_list_or_NULL'), function(object, value)
{
    object@validationData <- value
    return(object)
})

#' Show method for Classes Inheriting from `S4Model`
#'
#' @param object A `S4Model` derivative to show.
#'
#' @return None, prints to console.
#'
#' @importFrom CoreGx .collapse
#' @import S4Vectors
#' @export
setMethod('show', signature(object='S4Model'), function(object) {
    # -- Class
    cat('<', class(object)[1], '>', '\n')

    # -- trainData
    data <- trainData(object)
    cat('trainData: ', .collapse(dim(trainData(object))), '\n')
    cat('\t ', summary(data), '\n')
    if (length(assays(data)) > 0)
        cat('\t ', .collapse(names(assays(data))), '\n')

    # -- modelParams
    modParams <- modelParams(object)
    cat('modelParams: ', .collapse(dim(modParams)), '\n')
    cat(paste0('\t', capture.output(show(modParams)), collapse='\n\t'), '\n')

    # -- models
    cat('models: \n')
    mods <- models(object)
    cat(paste0('\t', capture.output(show(mods)), collapse='\n\t'), '\n')

    # -- validationStats
    cat('validationStats: \n')
    valStats <- validationStats(object)
    cat(paste0('\t', capture.output(head(valStats)), collapse='\n\t'), '\n')

    # -- validationData
    cat('validationData: \n')
    valData <- validationData(object)
    cat(paste0('\t', capture.output(show(valData)), collapse='\n\t'), '\n\n')

})