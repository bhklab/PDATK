#' `CohortList` Class
#'
#' A list containing only `SurvivalExperiment` objects.
#'
#' @keywords internal
.CohortList <- setClass('CohortList', contains='list')


#' Constructor for the `CohortList` class, a specialized list for stored two or
#'   more `SurvivalExperiment` objects.
#'
#' @param ... One or more `SurvivalExperiment` objects.
#'
#' @md
#' @export
CohortList <- function(...) {
    .CohortList(...)
}


#'
#'
#' @md
#' @export
setAs('list', 'CohortList', function(from) CohortList(from))


#' Class Validity Method for CohortList
#'
#' Ensure the CohortList abides the proper structure. This requires that all
#'   objects in it be of class `SurvivalExperiment.`
#'
#' @param object A `CohortList` to check the validity of.
#'
#' @md
#' @export
setValidity('CohortList', function(object) {
    isSurvivalExperiment <- vapply(object, FUN=is,
        class2='SurvivalExperiment', FUN.VALUE=logical(1))

    if (all(isSurvivalExperiment))
        return(TRUE)
    else  ## TODO:: Ensure that this returns the correct execution context
        return(.errorMsg(.context(3), 'The items at indexes ',
            paste0(which(!isSurvivalExperiment), collapse=', '),
            ' are not `SurvivalExperiment`s. A `CohortList` can only ',
            'contain objects of that class!'))
})