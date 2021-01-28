#'
#' @param model1 An `S4` object representing some kind of mathematical model.
#' @param model2 An `S4` object representing some kind of mathematical model.
#'
#' @return A `S4` object with statistics about the performance of each model.
#'
#' @md
#' @export
setGeneric('compareModels', function(model1, model2, ...)
    standardGeneric('compareModels'))
#'
#' @param model1 An object inherting from the `SurvivalModel` class.
#' @param model2 Another object inherting from the `SurvivalModel` class
#'
#' @return A `ModelComparison` object with statistics comparing the two models.
#'
#' @md
#' @export
setMethod('compareModels', signature(model1='SurvivalModel',
    model2='SurvivalModel'), function(model1, model2)
{
    ModelComparison(model1, model2)
})