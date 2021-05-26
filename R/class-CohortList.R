#' @name CohortList-class
#' `CohortList` Class Definition
#'
#' @description A list containing only `SurvivalExperiment` objects.
#'
#' @md
#' @importClassesFrom S4Vectors SimpleList
#' @export
.CohortList <- setClass('CohortList', contains='SimpleList',
    prototype=prototype(elementType='SurvivalExperiment'))

##FIXME:: CohortList constructor doesn't work if passed named items
## to dots... Review SimpleList class to see how they accomplish that.

#' @name CohortList
#' Constructor for the `CohortList` class, a specialized list for storing
#'   `SurvivalExperiment` objects.
#'
#' @param ... One or more `SurvivalExperiment` objects.
#' @param mDataTypes A `character` vector with the same length as ... which
#'   indicates the molecular data type in each cohort. This will be assigned
#'   to `mcols` for the `CohortList` as well as to the metadata of each
#'   item in the list as `mDataType`.
#'
#' @return A `CohortList` object containing one or more `SurvivalExperiment`
#'   objects.
#'
#' @examples
#' data(sampleICGCmicro)
#' set.seed(1987)
#' cohortList <- CohortList(list(survExp1=sampleICGCmicro,
#'   survExp2=sampleICGCmicro), mDataTypes=c('rna_micro', 'rna_micro'))
#'
#' @md
#' @importFrom S4Vectors metadata metadata<- mcols mcols<- mendoapply
#' @importClassesFrom S4Vectors SimpleList
#' @export
CohortList <- function(..., mDataTypes) {

    funContext <- .context(1)
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
    if (missing(mDataTypes)) {
        warning(.warnMsg(funContext, 'The mDataTypes for each cohort have ',
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
#' @return A `CohortList`.
#'
#' @md
#' @importClassesFrom S4Vectors SimpleList
#' @export
setAs('SimpleList', 'CohortList', function(from) CohortList(from))

## TODO:: Add mcols and metadata as attributes to returned list so that it can
## be coerce back to a  CohortList, then defined another as method

#' Coerce a `CohortList` to a `list`
#'
#' @param from A `CohortList` object.
#'
#' @return A `list` containing the contents of `from`.
#'
#' @md
#' @export
setAs('CohortList', 'list', function(from) from@listData)

#' Class Validity Method for `CohortList`
#'
#' Ensure the CohortList abides the proper structure. This requires that all
#'   objects in it be of class `SurvivalExperiment`.
#'
#' @param object A `CohortList` to check the validity of.
#'
#' @return TRUE if the class if valid, otherwise errors with a message.
#'
#' @md
#' @export
setValidity('CohortList', function(object) {
    isSurvivalExperiment <- vapply(object, FUN=is,
        class2='SurvivalExperiment', FUN.VALUE=logical(1))

    funContext <- .context(1)
    if (all(isSurvivalExperiment)) {
        TRUE
    } else {
        .errorMsg(funContext, 'The items at indexes ',
            paste0(which(!isSurvivalExperiment), collapse=', '),
            ' are not `SurvivalExperiment`s. A `CohortList` can only ',
            'contain objects of that class!')
    }
})
