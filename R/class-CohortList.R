#' `CohortList` Class
#'
#' A list containing only `SurvivalExperiment` objects.
#'
#' @importClassesFrom S4Vectors SimpleList
#' @keywords internal
.CohortList <- setClass('CohortList', contains='SimpleList',
    prototype=prototype(elementType='SurvivalExperiment'))


#' Constructor for the `CohortList` class, a specialized list for storing
#'   `SurvivalExperiment` objects.
#'
#' @param ... One or more `SurvivalExperiment` objects.
#' @param mDataTypes A `character` vector with the same length as ... which
#'   indicates the molecular data type in each cohort. This will be assigned
#'   to `mcols` for the `CohortList` as well as to the metadata of each
#'   item in the list as `mDataType`.
#'
#' @md
#' @importFrom S4Vectors metadata metadata<- mcols mcols<- mendoapply
#' @importClassesFrom S4Vectors SimpleList
#' @export
CohortList <- function(..., mDataTypes) {
    cohortList <- .CohortList(...)

    # Use existing mDataTypes if they exist
    if (missing(mDataTypes) && !is.null(mcols(cohortList)$mDataType) &&
            !any(is.na(mcols(cohortList)$mDataType)))
    {
        mDataTypes <- mcols(cohortList)$mDataType
    }

    # Try to get mDataTypes from the item metadata
    tryCatch({
        metadata <- lapply(cohortList, metadata)
        mData <- unlist(lapply(metadata, `[[`, i='mDataType'))
        if (!(is.null(mData) || length(mData) < 1 || any(is.na(mData)))) {
            if (missing(mDataTypes)) mDataTypes <- mData
        }
    })

    # warning if we couldn't get the mDataTypes
    ## FIXME:: .context() doesn't find the correct execution context
    if (missing(mDataTypes)) {
        warning(.warnMsg(.context(), 'The mDataTypes for each cohort have ',
            'not been set. Please set them manually using ',
            'mcols(cohortList)$mDataType <- c("rna_seq", "rna_micro", etc.)'))
    } else {
        mcols(cohortList)$mDataType <- mDataTypes
        if (!(all(mData == mDataTypes))) {
            for (i in seq_along(cohortList))
                metadata(cohortList[[i]]) <- mDataTypes[i]
        }
    }
    return(cohortList)
}


#' Coerce a `SimpleList` to a `CohortList`
#'
#' @param from A `SimpleList` object
#'
#' @md
#' @importClassesFrom S4Vectors SimpleList
#' @export
setAs('SimpleList', 'CohortList', function(from) CohortList(from))

## TODO:: Add mcols and metadata as attributes to returned list so that it can
## be coerce back to a  CohortList, then defined another as method

#'
#' @param from A `CohortList` object.
#'
#' @md
#' @export
setAs('CohortList', 'list', function(from) from@listData)
#' @export
as.list.CohortList <- function(from) as(from, 'list')


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

#'
#'
#'
#'
