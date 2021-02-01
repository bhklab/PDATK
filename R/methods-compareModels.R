#' Compare Two Mathematical Models Represented as `S4` Objects
#'
#' @param model1 A `S4` object representing some kind of mathematical model.
#' @param model2 A `S4` object representing some kind of mathematical model.
#'
#' @return A `S4` object with statistics about the performance of each model.
#'
#' @md
#' @export
setGeneric('compareModels', function(model1, model2, ...)
    standardGeneric('compareModels'))
#'
#' Compare Two `SurivalModel` Objects, Returing A `ModelComparison` Object
#'   With Statistics Comparing the Performance of Each Model.
#'
#' @param model1 An object inherting from the `SurvivalModel` class.
#' @param model2 Another object inherting from the `SurvivalModel` class
#' @param modelNames Optional character vector with a name for each model.
#'   Defaults to the class of the model plus 1 and 2 if missing.
#'
#' @return A `ModelComparison` object with statistics comparing the two models.
#'
#' @md
#' @export
setMethod('compareModels', signature(model1='SurvivalModel',
    model2='SurvivalModel'), function(model1, model2, modelNames)
{
    # deal with making model names unique when comparing two models of the
    # same type
    if (missing(modelNames)) {
        if (class(model1) != class(model2)) {
            modelNames <- c(class(model1), class(model2))
        } else {
            modelNames <- c(paste0(class(model1), '_', 1),
                paste0(class(model2), '_', 2))
        }
    }

    validationStats(model1)$model_name <- modelNames[1]
    validationStats(model2)$model_name <- modelNames[2]

    ModelComparison(model1, model2)
})

#'
#'
#'
#' @md
#' @export
setMethod('compareModels', signature(model1='ModelComparison',
    model2='SurvivalModel'), function(model1, model2, model2Name)
{
    if (missing(model2Name)) {
        model2Name <- class(model2)
        modelCompDT <- as.data.table(model1)
        # count the number of existing model name to ensure
        # the class + number combination is unique
        if (model2Name %in% modelCompDT$models)
            model2Name <- paste0(model2Name, '_',
                DT[model == model2Name, length(unique(model_name))])
    }

    validationStats(model2)$model_name <- model2Name

    ModelComparison(model1, model2)
})