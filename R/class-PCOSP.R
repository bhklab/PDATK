#' Pancreatic Cancer Overall Survival Predictor (PCOSP) Class
#'
#' @md
#' @include class-SurvivalModel.R
#' @export
.PCOSP <- setClass("PCOSP", contains='SurvivalModel')

#' Pancreatic Cancer Overall Survival Predictor (PCOSP) Constructor
#'
#' @details This function assumes there is only 1 assay per
#'  `SurvivalExperiment`.
#'
#' @param trainCohorts A `CohortList` or `SurivalExperiment` containing the
#'   training data for the `PCOSP` model.
#' @param minDaysSurvived An `integer` indicating the minimum number of day
#'   required to be in the 'good' survival group. Any patients below this
#'   cut-off will be considered low survival. Default is 365 days.
#' @param ... Force subsequent parameters to be named. This parameter is not
#'   used.
#' @param randomSeed An `integer` randomSeed that was used to train the model.
#'   Users should specify this when initializing a model to ensure
#'   reproducibilty.
#'
#' @return A `PCOSP` object with training data in the assays slot, concatenating
#'   together the molecular data types and labelling the genes with the data
#'   type to ensure the results are easily interpretable.
#'
#' @examples
#' data(sampleICGCmicro)
#' PCOSPmodel <- PCOSP(sampleICGCmicro, minDaysSurvived=365, randomSeed=1987)
#'
#' @md
#' @export
PCOSP <- function(trainCohorts, minDaysSurvived=365, ..., randomSeed) {

    PCOSPmodel <- .PCOSP(SurvivalModel(trainCohorts, minDaysSurvived,
        randomSeed=randomSeed))

    return(PCOSPmodel)
}
