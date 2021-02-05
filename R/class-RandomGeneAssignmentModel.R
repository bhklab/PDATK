#' RGAModel Class Definition
#'
#' @export
.RGAModel <- setClass('RGAModel', contains='SurvivalModel')
#'
#' RandomGeneAssignmentModel Constructor
#'
#' @inherit SurvivalModel
#'
#' @examples
#' data(sampleICGCmicro)
#' RGAmodel <- RGAModel(sampleICGCmicro, minDaysSurvived=365, randomSeed=1987)
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