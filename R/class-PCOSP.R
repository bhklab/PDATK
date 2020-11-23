#' Pancreatic Cancer Overall Survival Predictor (PCOSP) Class
#'
#' @slot trainCohorts A `CohortList` of `SurvivalExperiment`s for training the
#'   PCOSP model.
#' @slot model `` The trained PCOSP model.
#' @slot .intern An `environment` storing internal object metadata. For
#'   developer use.
#'
#' @md
#' @include class-SurvivalExperiment.R
#' @export
## TODO:: Can I subclass a BiocConductor object for this?
setClass("PCOSP", contains='SurvivalExperiment',
    slots=list(model='ANY'))

#' Pancreatic Cancer Overall Survival Predictor (PCOSP) Constructor
#'
#' @details This function assumes there is only 1 assay per `SurvivalExperiment`.
#'
#' @param trainCohorts A `CohortList` of `SurvivalExperiment`s for training
#'   the PCOSP model, where each cohort has a unique molecular data type. It
#'   is recommended to remove early deaths using `dropNotCensored` function
#'   before starting to build your PCOSP model.
#'
#' @return A `PCOSP` object with training data in the `trainCohorts` slot. The
#'   model slot will be empty until `trainModel` is called on the object.
#'
#' @md
#' @export
PCOSP <- function(trainCohort) {

    if (!is(trainCohorts, 'SurvivalExperiment'))
        stop(.errorMsg(.context(), 'The trainCohorts argument is not a ',
            'SurvivalExperiment object. Please convert it before building a ',
            'PCOSP model!'))

    assaysL <- assays(trainCohort)
    for (i in seq_along(assaysL)) {
        rownames(assaysL[[i]]) <- paste0(names(assaysL)[i], rownames(assaysL),
            sep='.')
    }
    modelMatrix <- do.call(rbind, assaysL)
    colData(trainCohort)$survival_group <-  # split into high and low survival
        ifelse(colData(trainCohort)$days_survived >= 365, 1, 0)

    new('PCOSP', colData=,
        rowData(rowData(trainCohort)), assays=modelMatrix,
        metadata=metadata(trainCohort))
}

##'
##' @importFrom settings is_setting clone_and_merge
##' @export
#setMethod('manageOptions', 'PCOSP', function(where=NULL, ...) {
#  if (settings::is_setting(...)) {
#    where@options <- settings::clone_and_merge(where@options(), ...)
#    where
#  } else {
#    where@options(...)
#  }
#})
