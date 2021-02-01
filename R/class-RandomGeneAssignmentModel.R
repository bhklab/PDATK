#' @inherit SurvivalModel
#'
#' @export
.RGAModel <- setClass('RGAModel', contains='SurvivalModel')
#'
#'
#' @inherit SurvivalModel
#'
#' @aliases RGAModel
#' @export
RandomGeneAssignmentModel <- function(trainCohorts, minDaysSurvived=365, ...,
    randomSeed)
{
    RGAmodel <- .RGAModel(SurvivalModel(trainCohorts, minDaysSurvived,
        randomSeed=randomSeed))
    return(RGAmodel)
}
#' @export
RGAModel <- RandomGeneAssignmentModel