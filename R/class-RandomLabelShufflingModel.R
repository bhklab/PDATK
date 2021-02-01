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
#' @aliases RLSModel
#' @export
RandomLabelShufflingModel <- function(trainCohorts, minDaysSurvived=365, ...,
    randomSeed)
{
    RLSModel <- .RLSModel(SurvivalModel(trainCohorts, minDaysSurvived,
        randomSeed=randomSeed))
    return(RLSModel)
}
#' @export
RLSModel <- RandomLabelShufflingModel