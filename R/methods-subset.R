#' Subset method for a `CohortList`
#'
#' Works using lapply of `[` over the list `SurvivalExperiment`s
#'
#' @param x A `CohortList` object
#' @param i The row query. Defaults to TRUE, i.e., select all.
#' @param j The column query. Defaults to TRUE, i.e., select all.
#'
#' @return A `CohortList` containing only the rows and columns selected in i
#'   and j, respectively.
#'
#' @importMethodsFrom S4Vectors subset
#' @export
setMethod('subset', signature(x='CohortList'),
    function(x, i=TRUE, j=TRUE)
{
    return(CohortList(lapply(x, `[`, i=i, j=j)))
})