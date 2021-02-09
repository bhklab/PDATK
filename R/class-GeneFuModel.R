#' @title A `SurvivalModel` Sub-class Designed to Hold A Survival Model Generated
#'   Using the `genefu` R package.
#'
#' @md
#' @export
.GeneFuModel <- setClass('GeneFuModel', contains='SurvivalModel')

#' `GeneFuModel` Constructor Method
#'
#' @param trainCohorts A `CohortList` or `SurvivalExperiment` containing
#'   training data for the genefu model. If you don't have training data,
#'   but have a trained model this will default to an empty `SurvivalExperiment`.
#'   You can then assign the model using the `models` setter method.
#' @param minDaysSurvived An `integer` specifying the minimum days survived
#'   to be considered in the 'good' survival prognosis group.
#' @param ... Fall through paramater to `SurvivalModel` constructor.
#' @param randomSeed A random seed to use for any sampling done in when
#'   training or evaluating your model.
#'
#' @return A `GeneFuModel` object, with model parameters in the
#'
#' @examples
#' geneFuModel <- GeneFuModel(randomSeed=1987)
#'
#' @md
#' @export
GeneFuModel <- function(trainCohorts=SurvivalExperiment(),
    minDaysSurvived=365, ..., randomSeed)
{

    survModel <- SurvivalModel(trainCohorts, minDaysSurvived=minDaysSurvived,
        ..., randomSeed=randomSeed)

    return(.GeneFuModel(survModel))
}