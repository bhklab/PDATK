#' Pancreatic Cancer Overall Survival Predictor (PCOSP) Class
#'
#' @slot trainCohorts A `CohortList` of `SurvivalExperiment`s for training the
#'   PCOSP model.
#' @slot models `SimpleList`
#' @slot
#'
#' @md
#' @include class-SurvivalModel.R
#' @export
.PCOSP <- setClass("PCOSP", contains='SurvivalModel')

#' Pancreatic Cancer Overall Survival Predictor (PCOSP) Constructor
#'
#' @details This function assumes there is only 1 assay per `SurvivalExperiment`.
#'
#' @param trainCohorts A `CohortList`s for training
#'   the PCOSP model, with.
#' @param minDaysSurvived An `integer` indicating the minimum number of day
#'   required to be in the 'high' survival group. Any patients below this
#'   cut-off will be considered low survival. Default is 365 days.
#' @param ... Force subsequent parameters to be named. This parameter is not
#'   used.
#' @param randomSeed An `integer` random seed to use for training the model. If
#'   missing, this defaults to the seed 1234.
#'
#' @return A `PCOSP` object with training data in the assays slot, concatenating
#'   together the molecular data types and labelling the genes with the data
#'   type to ensure the results are easily interpretable.
#'
#' @md
#' @export
PCOSP <- function(trainCohorts, minDaysSurvived=365, ..., randomSeed) {

    if (!is(trainCohort, 'SurvivalExperiment'))
        stop(.errorMsg(.context(), 'The trainCohorts argument is not a ',
            'SurvivalExperiment object. Please convert it before building a ',
            'PCOSP model!'))

    assaysL <- assays(trainCohort)

    modelMatrix <- do.call(rbind, as.list(assaysL))
    colData(trainCohort)$prognosis <-  # split into high and low survival
        ifelse(colData(trainCohort)$days_survived >= minDaysSurvived,
            'good', 'bad')

    rowData <- rbind(rowData(trainCohort), rowData(trainCohort))
    rownames(rowData) <- rownames(modelMatrix)

    survExp <- SurvivalExperiment(colData=colData(trainCohort),
        rowData=rowData, assays=SimpleList(trainMatrix=modelMatrix),
        metadata=metadata(trainCohort))

    PCOSPmodel <- .PCOSP(survExp)
    metadata(PCOSPmodel)$modelParams <-
        list(randomSeed=if (!missing(randomSeed)) randomSeed else 1234,
            RNGkind=RNGkind())
    return(PCOSPmodel)
}