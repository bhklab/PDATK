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
#' @alias RGAmodel
#' @export
RandomGeneAssignmentModel <- function(trainCohorts, minDaysSurvived=365, ...,
    randomSeed)
{
    RGAmodel <- .RGAmodel(trainCohorts, minDaysSurvived, randomSeed=randomSeed)
    return(RGAmodel)
}
RGAModel <- RandomGeneAssignmentModel


    #
    # .randomSampleRows <- function(x, size, replace=FALSE) {
    #     x[sample.int(nrow(x), size=size, replace=replace), ]
    # }
    # dataMatrix <- assay(trainCohorts)
    #
    # .oldSeed <- .Random.seed
    # on.exit(set.seed(.oldSeed))
    # set.seed(randomSeed)
    #
    # trainData <- lapply(seq_len(numModels), function(i, ...)
    #     .randomSampleRows(...), x=dataMatrix, size=numGenes, replace=FALSE)
    #
    # assays(trainCohorts) <- trainData