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


#' Coerce a `list` to a `CohortList`
#'
#' @param from A `list` of `SurvivalExperiment` objects to coerce to a
#'   `CohortList`.
#'
#' @md
#' @export
setAs('list', 'CohortList', function(from) CohortList(from))

#'
#'
#' @param from A `CohortList` to coerce to a `list`.
#'
#' @md
#' @export
setAs('CohortList', 'list', function(from)
    structure(from@.Data, .Names=names(from)))
#' @export
as.list.CohortList <- function(object) return(as(object, 'list'))

#' Class Validity Method for CohortList
#'
#' Ensure the CohortList abides the proper structure. This requires that all
#'   objects in it be of class `SurvivalExperiment`.
#'
#' @param object A `CohortList` to check the validity of.
#'
#' @md
#' @export
setValidity('CohortList', function(object) {
    isSurvivalExperiment <- vapply(object, FUN=is,
        class2='SurvivalExperiment', FUN.VALUE=logical(1))

    if (all(isSurvivalExperiment))
        TRUE
    else  ## TODO:: Ensure that this returns the correct execution context
        .errorMsg(.context(), 'The items at indexes ',
            paste0(which(!isSurvivalExperiment), collapse=', '),
            ' are not `SurvivalExperiment`s. A `CohortList` can only ',
            'contain objects of that class!')
})