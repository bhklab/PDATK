#'
#'
#'
#' @inherit .SurvivalModel
#'
#' @export
.RLSModel <- setClass('RLSModel', contains='SurvivalModel')
#'
#'
#'
#'
#' @alias RLSModel
#' @export
RandomLabelShuffingModel <- function(trainCohorts, minDaysSurvived=365, ...,
    randomSeed)
{
    RLSModel <- .RLSModel(SurvivalModel(trainCohorts, minDaysSurvived,
        randomSeed=randomSeed))
    return(RLSModel)
}
#' @export
RLSModel <- RandomLabelShuffingModel