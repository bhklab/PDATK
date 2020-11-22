#' Pancreatic Cancer Overall Survival Predictor (PCOSP) Class
#'
#'
#' @slot trainCohorts A `CohortList` of `SurvivalExperiment`s for training the
#'   PCOSP model.
#' @slot model `` The trained PCOSP model.
#' @slot .intern An `environment` storing internal object metadata. For
#'   developer use.
#'
#' @md
#' @include options.R
#' @export
## TODO:: Can I subclass a BiocConductor object for this?
setClass("PCOSP", slots=list(trainCohorts='CohortList',
    modelMatrix='matrix',  model='ANY', metadata='list',
    .intern='environment'))

#' Pancreatic Cancer Overall Survival Predictor (PCOSP) Constructor
#'
#' @param trainCohorts A `CohortList` of `SurvivalExperiment`s for training
#'   the PCOSP model.
#'
#' @return A `PCOSP` object with training data in the `trainCohorts` slot. The
#'   model slot will be empty until `trainModel` is called on the object.
#'
#' @md
#' @export
PCOSP <- function(trainCohorts) {

    if (!is(trainCohorts, 'CohortList'))
        stop(.errorMsg(.context(), 'The trainCohorts argument is not a ',
            'CohortList object. Please convert it before building a PCOSP ',
            'model!'))

    ## TODO:: Ensure the SurvivalExperiments have different mDataTypes?

    # ensure the SurivalExperiments have common genes
    commonGenes <- findCommonGenes(trainCohorts)
    actualGenes <- Reduce(intersect, lapply(trainCohorts, rownames))

    # subset to common genes
    if (!all(actualGenes %in% commonGenes)) {
        trainCohorts <- subset(trainCohorts, subset=commonGenes)
        warning(.warnMsg(.context(), 'The training cohorts did not have only ',
            'common genes. Subsetting to common genes...'))
    }

    # ensure the SurvivalExperiments have common samples
    commonSamples <- findCommonSamples(trainCohorts)
    actualSamples <- Reduce(intersect, lapply(trainCohorts, colnames))

    # subset to common samples
    if (!all(actualSamples %in% commonSamples)) {
        trainCohorts <- subset(trainCohorts, select=commonSamples)
        warning(.warnMsg(.context(), 'The training cohorts did not have only ',
                         'common samples Subsetting to common samples...'))
    }
    new('PCOSP', trainCohorts=trainCohorts)
}

#'
#' @importFrom settings is_setting clone_and_merge
#' @export
setMethod('manageOptions', 'PCOSP', function(where=NULL, ...) {
  if (settings::is_setting(...)) {
    where@options <- settings::clone_and_merge(where@options(), ...)
    where
  } else {
    where@options(...)
  }
})

## TODO:: Validity method?

#' #'
#' #'
#' #'
#' #'
#' #'
#' setValidity('PCOSP', function(object) {
#'
#' })