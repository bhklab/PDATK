#' `CohortList` Class
#'
#' A list containing only `SurvivalExperiment` objects.
#'
#' @importClassesFrom S4Vectors SimpleList
#' @keywords internal
setClass('CohortList', contains='SimpleList',
    prototype=prototype(elementType='SurvivalExperiment'))


#' Constructor for the `CohortList` class, a specialized list for stored two or
#'   more `SurvivalExperiment` objects.
#'
#' @param ... One or more `SurvivalExperiment` objects.
#'
#' @md
#' @importFrom S4Vectors metadata metadata<- mcols mcols<-
#' @export
CohortList <- function(...) {
    new('CohortList', ...)
}


#' Coerce a `list` to a `CohortList`
#'
#' @param from A `list` of `SurvivalExperiment` objects to coerce to a
#'   `CohortList`.
#'
#' @md
#' @export
setAs('list', 'CohortList', function(from) CohortList(from))

#' Coerce a `SimpleList` to a `CohortList`
#'
#' @param from A `SimpleList` object
#'
#' @md
#' @importClassesFrom S4Vectors SimpleList
#' @export
setAs('SimpleList', 'CohortList', function(from) CohortList(from))

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

    if (all(isSurvivalExperiment)) {
        TRUE
    } else {
        .errorMsg(.context(), 'The items at indexes ',
            paste0(which(!isSurvivalExperiment), collapse=', '),
            ' are not `SurvivalExperiment`s. A `CohortList` can only ',
            'contain objects of that class!')
    }
})