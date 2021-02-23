#' SurvivalExperiment Class
#'
#' @description
#' A SummarizedExperiment with mandatory numeric survival metadata columns
#'   `survival_time` and `event_occurred`.
#'
#' @md
#' @export
.SurvivalExperiment <- setClass('SurvivalExperiment',
    contains='SummarizedExperiment')

#' @title Constructor for `SurvivalExperiment` Class
#'
#' Builds a `SurvivalExperiment` object, which is just a wrapper for a
#'   `SummarizedExperiment` with mandatory survival metadata numeric columns
#'   `survival_time` and `event_occurred`.
#'
#' @param ... `pairlist` Fall through arguments to the `SummarizedExperiment`
#'   constructor. If the first argument to dots is a `SummarizedExperiment`,
#'   that object is used instead.
#' @param survival_time A `character` vector indicating the column name in
#'   `colData` which contains the `integer` number of days a patient
#'   has survived since treatment at the time of data collection. If
#'   `event_occurred` is 1/TRUE, then this is the number of days the patient lived.
#' @param event_occurred A `character` vector indicating the column name in
#'   `colData` which contains `logical` or `integer` values where 0/FALSE means
#'   a patient is alive and 1/TRUE means a patient is deceased.
#'
#' @return A `SurvivalExperiment` object.
#'
#' @examples
#' data(sampleICGCmicro)
#'
#' # build a SurvivalExperiment from raw data
#' ICGCmicro <- SurvivalExperiment(assays=assays(sampleICGCmicro),
#'   rowData=rowData(sampleICGCmicro), colData=colData(sampleICGCmicro),
#'   metadata=metadata(sampleICGCmicro), survival_time='survival_time',
#'   event_occurred='event_occurred')
#'
#' # build a SurvivalExperiment from an existig SummarizedExperment
#' ICGCmicroSumExp <- as(sampleICGCmicro, 'SummarizedExperiment')
#' ICGCmicro <- SurvivalExperiment(sumExp=ICGCmicroSumExp,
#'   survival_time='survival_time', event_occurred='event_occurred')
#'
#' @md
#' @importFrom SummarizedExperiment SummarizedExperiment colData colData<-
#' @importFrom S4Vectors rename
#' @importFrom CoreGx .errorMsg .warnMsg
#' @export
SurvivalExperiment <- function(..., survival_time='survival_time',
                               event_occurred='event_occurred')
{
    funContext <- .context(1)

    ## TODO:: Clean up constructor logic
    dots <- list(...)
    SE <- if (is(dots[[1]], 'SummarizedExperiment')) dots[[1]] else
        SummarizedExperiment(...)

    if (!is(SE, 'SummarizedExperiment'))
        stop(CoreGx::.errorMsg(funContext,
                               'sumExp is not a `SummarizedExperiment`!'))

    renameVector <- c('survival_time', 'event_occurred')
    names(renameVector) <- c(survival_time, event_occurred)

    # allow empty SurivalExperiments to exist
    if (nrow(colData(SE)) == 0) {
        if (!all(renameVector %in% colnames(colData(SE)))) {
            colData(SE) <- cbind(colData(SE),
                                 DataFrame(survival_time=integer(), event_occurred=integer()))
        }
    }

    hasColumnsToRename <- names(renameVector) %in% colnames(colData(SE))
    if (all(hasColumnsToRename)) {
        colData(SE) <- rename(colData(SE), renameVector)
    } else {
        stop(.errorMsg(funContext, 'The columns ',
                       names(renameVector)[!hasColumnsToRename], ' are not present in ',
                       'the object colData, please ensure you specify existing column',
                       'names to the days_surived and event_occurred parameters!'))
    }

    if (!is.integer(colData(SE)$event_occurred)) {
        event_occurred_col <- colData(SE)$event_occurred
        switch(class(event_occurred_col),
               'logical'={ colData(SE)$event_occurred <- as.integer(event_occurred_col) },
               'character'={
                   if (!('deceased' %in% event_occurred_col))
                       stop(.errorMsg(funContext, 'The string deceased is not in ',
                                      'the event_occurred column. Please convert this column to ',
                                      'integer manually, where 1 is deceased and 0 is alive.'))
                   colData(SE)$event_occurred <-
                       as.integer(event_occurred_col == 'deceased')
               },
               stop(.errorMsg(funContext, 'The event_occurred column is not logical ',
                              'or integer, please convert this column such that deceased is 1 ',
                              'or TRUE and alive is 0 or FALSE before retrying the conversion!'))
        )
    }
    if (!is.integer(colData(SE)$survival_time)) {
        survival_time <- colData(SE)$survival_time
        switch(class(survival_time),
               'numeric'={ colData(SE)$survival_time <- as.integer(survival_time) },
               'character'={ tryCatch({
                   colData(SE)$survival_time <- as.integer(survival_time)
               },
               warning=function(w) stop(.errorMsg(funContext, 'Tried to ',
                                                  'coerce survival_time from character to integer, but ',
                                                  'failed.')),
               error=function(e) stop(.errorMsg(funContext, 'Tried to ',
                                                'coerce survival_time from character to integer, but ',
                                                'failed.')))
               },
               stop(.errorMsg(funContext, 'The survival_time column is not logical',
                              ' or integer, please convert this column before retrying the ',
                              'conversion'))
        )
    }

    return(.SurvivalExperiment(SE))
}

#' Coerce Method from SummarizedExperiment to SurvivalExperiment
#'
#' @param from A `SummarizedExperiment` to coerce to a `SurvivalExperiment`.
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
#' @param object A `SurvivalExperiment` object to verify class validity of.
#'
#' @md
#' @importFrom CoreGx .errorMsg
#' @export
setValidity('SurvivalExperiment', function(object) {
    funContext <- .context(1)
    validateSummarizedExperiment <-
        getValidity(getClassDef('SummarizedExperiment'))
    isValidSummarizedExperiment <- validateSummarizedExperiment(object)

    survivalColNames <- c("survival_time", "event_occurred")
    hasSurvivalCols <- survivalColNames %in% colnames(colData(object))
    if (!all(hasSurvivalCols))
        .errorMsg(funContext, 'Mandatory columns ',
            paste0(surivalColNames[!hasSurvivalCols], collapse=', '),
            ' are missing from colData. Please add them or double check
            the column names are spelled correctly.')
    else
        isValidSummarizedExperiment &&
            is.numeric(colData(object)$event_occurred) &&
            is.numeric(colData(object)$survival_time)
})
