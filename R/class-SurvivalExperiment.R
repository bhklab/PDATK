#' Survival Experiment Class
#'
#' A SummarizedExperiment with mandatory survival metadata in the `colData`
#'   slot.
#'
#' @inherit SummarizedExperiment::SummarizedExperiment
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.SurvivalExperiment <- setClass('SurvivalExperiment',
    contains='SummarizedExperiment')

#' Constructor for `SurvivalExperiment` Class
#'
#' Builds a `SurvivalExperiment` object, which is just a wrapper for a
#'   `SummarizedExperiment` with mandatory survival metadata in the `colData`
#'   slot.
#'
#' @param ... `pairlist` Fall through arguments to the `SummarizedExperiment`
#'   constructor. These are ignored if `sumExp` is specified
#' @param days_survived A `character` vector indicating the column name in
#'   `colData` which contains the `integer` number of days a patient
#'   has survived since treatment at the time of data collection. If
#'   `is_deceased` is 1/TRUE, then this is the number of days the patient lived.
#' @param is_deceased A `character` vector indicating the column name in
#'   `colData` which contains `logical` or `integer` values where 0/FALSE means
#'   a patient is alive and 1/TRUE means a patient is deceased.
#' @param sumExp An optional `SummarizedExperiment` object to coerce to a
#'   `SurvivalExperiment`. If this parameter is included, all arguments to `...`
#'   are ignored.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
SurvivalExperiment <- function(..., days_survived='days_to_death',
    is_deceased='vital_status', sumExp)
{
    SE <- if (missing(sumExp)) SummarizedExperiment(...) else sumExp

    if (!is(SE, 'SummarizedExperiment'))
        stop(CoreGx::.errorMsg(.context(),
            'sumExp is not a `SummarizedExperiment`!'))

    renameVector <- c('days_survived', 'is_deceased')
    names(renameVector) <- c(days_survived, is_deceased)

    colData(SE) <- rename(colData(SE), renameVector)

    if (!is.integer(colData(SE)$is_deceased)) {
        is_deceased <- colData(SE)$is_deceased
        switch(typeof(is_deceased),
            'logical'={ colData(SE)$is_deceased <- as.integer(is_deceased) },
            'character'={ tryCatch({
                colData(SE)$is_deceased <-
                    ifelse(is_deceased == 'deceased', 1, 0)
                },
                error=function(e) stop(.errorMsg(.context(), 'Tried to coerce ',
                    'is_deceased from character to integer, but failed.')))
            },
            stop(.errorMsg(.context(), 'The is_deceased column is not logical
              or integer, please convert this column such that deceased is 1
              or TRUE and alive is 0 or FALSE before retrying the conversion!'))
        )
    }

    return(.SurvivalExperiment(SE))
}

#' Coerce Method from SummarizedExperiment to SurvivalExperiment
#'
#' @param from A ` SummarizedExperiment` to coerce to a `SurvivalExperiment`.
#'
#' @md
#' @export
setAs('SummarizedExperiment', 'SurvivalExperiment',
    function(from) SurvivalExperiment(sumExp=from))
#' @export
setAs('RangedSummarizedExperiment', 'SurvivalExperiment',
    function(from) SurvivalExperiment(sumExp=from))

#' Check that a SurvivalExperiment object is valid
#'
#' @md
#' @export
setValidity('SurvivalExperiment', function(object) {

    validateSummarizedExperiment <-
        getValidity(getClassDef('SummarizedExperiment'))
    isValidSummarizedExperiment <- validateSummarizedExperiment(object)

    survivalColNames <- c("days_survived", "is_deceased")
    hasSurvivalCols <- survivalColNames %in% colnames(colData(object))
    if (!all(hasSurvivalCols))
        CoreGx::.errorMsg(.context(), 'Mandatory columns ',
            paste0(surivalColNames[!hasSurvivalCols], collapse=', '),
            ' are missing from colData. Please add them or double check
            the column names are spelled correctly.')
    else
        isValidSummarizedExperiment
})