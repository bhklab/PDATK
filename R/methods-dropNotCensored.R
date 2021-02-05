#' Remove Censored Patient Samples from An `S4` Object.
#'
#' @param object An `S4` object containing survival data which needs to have
#'   patients who were not censored before some criteria.
#' @param ... Allow new parmeters to be defined for this generic.
#'
#' @return `S4` The object subset to only those patients which pass the
#'   censoring criteria.
#'
#' @md
#' @export
setGeneric('dropNotCensored',
    function(object, ...) standardGeneric('dropNotCensored'))
#'
#' Remove Censored Patients from A `SurvivalExperiment` Object
#'
#' @param object A `SurvivalExperiment` to censor.
#' @param minDaysSurvived An `integer` specifying the minimum number of days
#'   a patient needs to have survived to be included in the cohort.
#'
#' @details
#' Censored means no event before end of measurement. Since we want not
#'   censored, we keep patients who had an event before minDaysSurvived.
#'   Therefore we keep individuals surviving > `minDaysSurvived`, or who had an
#'   event (died) before minDaysSurvived.
#'
#' @return The `SurvivalExperiment` with censored samples removed.
#'
#' @examples
#' data(sampleICGCmicro)
#' ICGCmicro <- dropNotCensored(sampleICGCmicro)
#'
#' @md
#' @importFrom S4Vectors na.omit
#' @importFrom SummarizedExperiment colData colData<-
#' @export
setMethod('dropNotCensored', signature('SurvivalExperiment'),
    function(object, minDaysSurvived=365)
{
    # drop NA rows
    object <- object[, !is.na(colData(object)$days_survived)]

    days_survived <- colData(object)$days_survived
    is_deceased <- colData(object)$is_deceased

    notCensoredBefore <- days_survived <= minDaysSurvived & is_deceased == 1
    notYearOne <- days_survived > minDaysSurvived

    keepPatients <- notCensoredBefore | notYearOne
    object <- object[, keepPatients]
    colData(object)$prognosis <-
        ifelse(days_survived[keepPatients] > 365, 'good', 'bad')
    return(object)
})
#'
#' Remove Censored Patients from Each `SurvivalExperiemnt` in a `CohortList`
#'
#' @param object A `CohortList` for which to drop patients who died before
#'   each `SurvivalExperiment` item a specified date.
#' @param minDaysSurvived An `integer` specifying the minimum number of days
#'   a patient needs to have survived to be included in the cohort.
#'
#' @return The `CohortList` with censored samples removed.
#'
#' @examples
#' data(sampleCohortList)
#' valCohortList <- dropNotCensored(sampleCohortList)
#'
#' @md
#' @importFrom S4Vectors endoapply
#' @export
setMethod('dropNotCensored', signature('CohortList'),
    function(object, minDaysSurvived=365)
{
    endoapply(object, dropNotCensored, minDaysSurvived=minDaysSurvived)
})