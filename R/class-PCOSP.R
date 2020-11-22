#' Pancreatic Cancer Overall Survival Predictor (PCOSP) Class
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
PCOSP <- function(trainCohorts) {

    if (!is(trainCohorts, 'CohortList'))
        stop(.errorMsg(.context(), 'The trainCohorts argument is not a ',
            'CohortList object. Please convert it before building a PCOSP ',
            'model!'))

    ## TODO:: Ensure the SurvivalExperiments have different mDataTypes?
    if (length(unique(mcols(trainCohorts)$mDataType)) <
        nrow(mcols(trainCohorts)))
    {
        stop(.errorMsg(.context(), 'The training cohorts have the same ',
            'molecular data type. Please manually merge all cohorts with the ',
            'same data type using some summary statistic such as mean, median, ',
            'etc., before passing into a PCOSP model.'))
    }

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

    # concatenate the datatypes to label genes
    cohortAssays <- lapply(trainCohorts, assay, 1)
    for (i in seq_along(cohortAssays)) {
        rownames(cohortAssays[[i]]) <-
            paste(mcols(trainCohorts)$mDataType[i],
                rownames(cohortAssays[[i]]), sep='.')
    }
    modelMatrix <- do.call(rbind, cohortAssays)

    .getSurvivalData <-
        function(x) c(colData(x)$days_survived, colData(x)$is_deceased)

    survivalData <- lapply(trainCohorts, .getSurvivalData)

    if (!all(survivalData[[1]] == survivalData[[2]]))
        stop(.errorMsg(.context(), 'The training cohorts for a PCOSP model ',
            'should be on the same samples, but the survival data in the ',
            'colData slot of the SurvivalExperiments are different. Please ',
            'ensure the molecular data in each survival experiment is the same!'))

    # get the survived data
    survivalData <- data.table(
        sample=
    )

    new('PCOSP', trainCohorts=trainCohorts, modelMatrix=modelMatrix)
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