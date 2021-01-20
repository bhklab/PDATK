#' Remove Deaths Due to Disease Severity of Illness from An S4 Object Where
#'   Columns Are Samples.
#'
#' @param object An `S4` object containing survival data which needs to have
#'   patients who were not censored before some criteria.
#'
#' @return `S4` The object subset to only those patients which pass the
#'   censoring criteria.
#'
#' @export
setGeneric('dropNotCensored',
    function(object, ...) standardGeneric('dropNotCensored'))
#'
#' @param object A `SurvivalExperiment` to censor.
#' @param minDaysSurvived An `integer` specifying the minimum number of days
#'   a patient needs to have survived to be included in the cohort.
#'
#' @export
setMethod('dropNotCensored', signature('SurvivalExperiment'),
    function(object, minDaysSurvived=365)
{
    # drop NA rows
    object <- object[, rownames(na.omit(colData(object)))]

    days_survived <- colData(object)$days_survived
    is_deceased <- colData(object)$is_deceased

    notCensoredBefore <- days_survived <= minDaysSurvived & is_deceased == 1
    notYearOne <- days_survived > minDaysSurvived

    keepPatients <- notCensoredBefore | notYearOne
    return(object[, keepPatients])
})
#'
#' @param object A `CohortList` for which to drop patients who died before
#'   each `SurvivalExperiment` item a specified date.
#' @param minDaysSurvived An `integer` specifying the minimum number of days
#'   a patient needs to have survived to be included in the cohort.
#'
#' @importFrom S4Vectors endoapply
#' @md
#' @export
setMethod('dropNotCensored', signature('CohortList'),
    function(object, minDaysSurvived=365)
{
    endoapply(object, dropNotCensored, minDaysSurvived=minDaysSurvived)
})