#' ModelComparison Class Definition
#'
#' @md
#' @importClassesFrom S4Vectors DataFrame
#' @export
.ModelComparison <- setClass('ModelComparison', contains='DataFrame')


#' ModelComparison Constructor
#'
#' @param model1 An object with a `validationStats` method which returns a
#'   `data.table`. Probably this object should inherit from `SurvivalModel`.
#' @param model2 An object with a `validationStats` method which returns a
#'   `data.table`. Probably this object should inherit from `SurvivalModel`.
#' @param ... Not used.
#'
#' @return A `ModelComparison` object, which is a soft wrapper for `DataFrame`
#'   which is used for method dispatch.
#'
#' @md
#' @importFrom methods is
#' @importFrom S4Vectors DataFrame
#' @import data.table
#' @export
ModelComparison <- function(model1, model2, ...) {


    ## TODO:: Is it better to define a validationStats method for a
    ##   ModelComparsion? Then can't do class for model column.
    if (is(model1, 'ModelComparison')) {
        model1StatsDT <- as.data.table(model1)
    } else {
        model1StatsDT <- validationStats(model1)
        model1StatsDT[, model := class(model1)]
    }

    if (is(model2, 'ModelComparison')) {

    } else {
        model2StatsDT <- validationStats(model2)
        model2StatsDT[, model := class(model2)]
    }

    sharedCohorts <- intersect(model1StatsDT$cohort, model2StatsDT$cohort)
    modelComparisonDT <- rbind(model1StatsDT, model2StatsDT)
    modelComparisonDT <- modelComparisonDT[cohort %in% sharedCohorts, ]

    .ModelComparison(DataFrame(modelComparisonDT))
}
