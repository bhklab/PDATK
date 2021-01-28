#' @inherit DataFrame
#'
#' @md
#' @export
.ModelComparison <- setClass('ModelComparison', contains='DataFrame')


#' @inherit DataFrame
#'
#' @param model1 An object with a `validationStats` method which returns a
#'   `data.table`. Probably this object should inherit from `SurvivalModel`.
#' @param model2 An object with a `validationStats` method which returns a
#'   `data.table`. Probably this object should inherit from `SurvivalModel`.
#'
#' @return A `ModelComparison` object, which is a soft wrapper for `DataFrame`
#'   which is used for method dispatch.
#'
#' @md
#' @export
ModelComparison <- function(model1, model2, ...) {

    model1StatsDT <- validationStats(model1)
    model2StatsDT <- validationStats(model2)

    model1StatsDT[, model := class(model1)]
    model2StatsDT[, model := class(model2)]

    sharedCohorts <- intersect(model1StatsDT$cohort, model2StatsDT$cohort)
    modelComparisonDT <- rbind(model1StatsDT, model2StatsDT)
    modelComparisonDT <- modelComparisonDT[cohort %in% sharedCohorts, ]

    .ModelComparison(DataFrame(modelComparisonDT))

}
