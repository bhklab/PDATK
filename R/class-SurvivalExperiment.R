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
#'   constructor. These are ignored if `SumExp` is specified
#' @param survivalColNames `character` A named character vector specifying
#'   which columns in colData contain the days_to_death and
#' @param sumExp
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
SurvivalExperiment <- function(..., sumExp,
    survivalColNames=c(days_to_death='days_to_death',
        vital_status='vital_status'))
{
    SE <- if (missing(sumExp)) SummarizedExperiment(...) else sumExp

    if (!is(SE, 'SummarizedExperiment'))
        stop(CoreGx::.errorMsg(.context(),
            'sumExp is not a `SummarizedExperiment`!'))

    hasSurvivalCols <- survivalColNames %in% colnames(colData(SE))
    if (!all(hasSurvivalCols))
        stop(.errorMsg(.context(), 'Mandatory columns ',
            paste0(surivalColNames[!hasSurvivalCols], collapse=', '),
            ' are missing from colData. Please add them or double check
            the column names are spelled correctly.'))

    survExp <- .SurvivalExperiment(SE)
    return(survExp)
}

#' Coerce Method from SummarizedExperiment to SurvivalExperiment
#'
#'
#' @export
setAs('SummarizedExperiment', 'SurvivalExperiment',
    function(from) SurvivalExperiment(sumExp=from))