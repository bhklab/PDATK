#'
#'
#'
#' @inherit .SurvivalModel
#'
#' @export
.RGAmodel <- setClass('RGAModel', contains='SurvivalModel')
#'
#'
#'
#'
#' @alias RGAModel
#' @export
RandomGeneAssignmentModel <- function(trainCohorts, minDaysSurvived=365, ...,
    randomSeed)
{
    RGAmodel <- .RGAmodel(SurvivalModel(trainCohorts, minDaysSurvived,
        randomSeed=randomSeed))
    return(RGAmodel)
}
#' @export
RGAModel <- RandomGeneAssignmentModel