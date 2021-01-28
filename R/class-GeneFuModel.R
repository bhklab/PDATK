#' A `SurvivalModel` Sub-class Designed to Hold A Survival Model Generated
#'   Using the `genefu` R package.
#'
#' @inherit SurvivalModel-class
#'
#' @md
#' @export
.GeneFuModel <- setClass('GeneFuModel', contains='SurvivalModel')


#'
#'
#'
#' @md
#' @export
GeneFuModel <- function(trainCohorts,
    minDaysSurvived=365, ..., randomSeed)
{
    survModel <- SurivalModel(trainCohorts, minDaysSurvived=minDaysSurvived,
        ..., randomSeed=randomSeed)

    return(.GeneFuModel(survModel))
}