#' Remove Deaths Due to Disease Severity from An S4 Object Where Columns Are
#'   Samples.
#'
#' @param object An `S4` object containing survival data to which needs to be
#'   censored.
#'
#' @return `S4` The object subset to only those patients which pass the censoring
#'   criteria.
#'
#' @export
setGeneric('censorSurvival',
    function(object, ...) standarndGeneric('censorSurvival'))
#'
#' @param object A `SurvivalExperiment` to censor.
#' @param minDaysSurvived The `integer` minimum number of days a patient needs
#'   to survived after treatment to not be censored. Default is 385.
#' @param removeDeceased Should all deceased patients be removed? Defaults to
#'   TRUE, since we are trying to predict survival.
#'
#' @export
setMethod('censorSurvival', signature('SurvivalExperiment'),
    function(object, minDaysSurvived=365, removeDeceased=TRUE)
{
    keepDays <- colData(object)$days_survived >= minDaysSurvived
    isDeceased <- colData(object)$is_deceased == 1

    keepPatients <- if (removeDeceased) keepDays & !isDeceased else keepDays
    return(object[, keepPatients])
})
#'
#' @param object A `CohortList` for which to sensor each `SurvivalExperiment`
#'   item.
#'
#' @importFrom S4Vectors endoapply
#' @md
#' @export
setMethod('censorSurvival', signature('CohortList'),
    function(object, minDaysSurvived=365, removeDeceased=TRUE)
{
    endoapply(object, censorSurvival, minDaysSurvived=minDaysSurvived,
        removeDeceased=removeDeceased)
})