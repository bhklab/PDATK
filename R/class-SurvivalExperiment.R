#' SurvivalExperiment Class
#'
#' A SummarizedExperiment with mandatory numeric survival metadata columns
#'   `days_survived` and `is_deceased`.
#'
#' @inherit SummarizedExperiment::SummarizedExperiment
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.SurvivalExperiment <- setClass('SurvivalExperiment',
    contains='SummarizedExperiment')

#' Constructor for `SurvivalExperiment` Class
#'
#' Builds a `SurvivalExperiment` object, which is just a wrapper for a
#'   `SummarizedExperiment` with mandatory survival metadata numeric columns
#'   `days_survived` and `is_deceased`.
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
SurvivalExperiment <- function(..., days_survived='days_survived',
    is_deceased='is_deceased', sumExp)
{
    ## TODO:: Clean up constructor logic
    SE <- if (missing(sumExp)) SummarizedExperiment(...) else sumExp

    if (!is(SE, 'SummarizedExperiment'))
        stop(CoreGx::.errorMsg(.context(),
            'sumExp is not a `SummarizedExperiment`!'))

    renameVector <- c('days_survived', 'is_deceased')
    names(renameVector) <- c(days_survived, is_deceased)

    if (all(names(renameVector) %in% colnames(colData(SE))))
        colData(SE) <- rename(colData(SE), renameVector)

    # allow empty SurivalExperiments to exist
    if (nrow(colData(SE)) == 0) {
        if (!all(renameVector %in% colnames(colData(SE)))) {
            colData(SE) <- cbind(colData(SE),
                DataFrame(days_survived=integer(), is_deceased=integer()))
        }
    }

    if (!is.integer(colData(SE)$is_deceased)) {
        is_deceased_col <- colData(SE)$is_deceased
        switch(class(is_deceased_col),
            'logical'={ colData(SE)$is_deceased <- as.integer(is_deceased_col) },
            'character'={
                if (!('deceased' %in% is_deceased_col))
                    stop(.errorMsg(.context(), 'The string deceased is not in ',
                        'the is_deceased column. Please convert this column to
                        integer manually, where 1 is deceased and 0 is alive.'))
                colData(SE)$is_deceased <-
                    as.integer(is_deceased_col == 'deceased')
            },
            stop(.errorMsg(.context(), 'The is_deceased column is not logical
              or integer, please convert this column such that deceased is 1
              or TRUE and alive is 0 or FALSE before retrying the conversion!'))
        )
    }
    if (!is.integer(colData(SE)$days_survived)) {
        days_survived <- colData(SE)$days_survived
        switch(class(is_deceased),
            'numeric'={ colData(SE)$days_survived <- as.integer(days_survived) },
            'character'={ tryCatch({
                colData(SE)$is_deceased <- as.integer(is_deceased)
                },
                warning=function(w) stop(.errorMsg(.context(), 'Tried to ',
                    'coerce days_survived from character to integer, but ',
                    'failed.')))
            },
            stop(.errorMsg(.context(), 'The days_survived column is not numeric
              or integer, please convert this column before retrying the
              conversion'))
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
        isValidSummarizedExperiment &&
            is.numeric(colData(object)$is_deceased) &&
            is.numeric(colData(object)$days_survived)
})
