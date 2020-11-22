#' Subset method for a `CohortList`
#'
#' Works using endoapply of `[` over the list `SurvivalExperiment`s
#'
#' @param x A `CohortList` object
#' @param subset The row query. Defaults to TRUE, i.e., select all.
#' @param select The column query. Defaults to TRUE, i.e., select all.
#' @param invert A `logical` vector indicating if the matches should be
#'   inverted. Default is FALSE.
#'
#' @return A `CohortList` containing only the rows and columns selected in i
#'   and j, respectively.
#'
#' @importMethodsFrom S4Vectors subset
#' @importFrom S4Vectors endoapply
#' @md
#' @export
setMethod('subset', signature(x='CohortList'),
    function(x, subset=TRUE, select=TRUE, invert=FALSE)
{
    if (invert) {
        rows <- lapply(x, rownames)
        keepRows <- lapply(rows, setdiff, y=subset)
        cols <- lapply(x, colnames)
        keepCols <- lapply(cols, setdiff, y=select)
        return(mendoapply(`[`, x=x, i=keepRows, j=keepCols))
    }

    return(endoapply(x, `[`, i=subset, j=select))
})