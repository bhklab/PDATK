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
#' @export
setClass("PCOSP", slots=list(trainCohorts='CohortList', model='ANY', ##TODO:: What class is model?
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

    # ensure the SurivalExperiments have common genes
    commonGenes <- findCommonGenes(trainCohorts)
    actualGenes <- Reduce(intersect, lapply(trainCohorts, rownames))

    # subset to common genes
    if (!all(actualGenes %in% commonGenes)) {
        trainCohorts <- subset(trainCohorts, subset=commonGenes)
        warning(.warnMsg(.context(), 'The training cohorts did not have only ',
            'common genes. Subsetting to common genes...'))
    }

    # remove

    new('PCOSP', trainCohorts=trainCohorts)
}

## TODO:: Validity method?

#' #'
#' #'
#' #'
#' #'
#' #'
#' setValidity('PCOSP', function(object) {
#'
#' })