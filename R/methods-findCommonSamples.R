#' Find Common Samples in a List-like S4 Object where The Columns of Each Item
#'   Represent Samples
#'
#' @param object A `S4` object, where the columns of each element represent
#'   samples.
#' @param ... Allow new parameters to be defined for this generic.
#'
#' @return A `character` vector of common sample names.
#'
#' @examples
#' data(sampleCohortList)
#' commonSamples <- findCommonSamples(sampleCohortList)
#' head(commonSamples)
#'
#' @md
#' @export
setGeneric('findCommonSamples',
    function(object, ...) standardGeneric('findCommonSamples'))
#'
#' Find Common Samples in a CohortList Object where The Columns of Each Item
#'   Represent Samples
#'
#' @param object A `CohortList` for which we want to find common samples
#'   between all `SurvivalExperiment` objects.
#'
#' @return A `character` vector of common sample names.
#'
#' @examples
#' data(sampleCohortList)
#' commonSamples <- findCommonSamples(sampleCohortList)
#' head(commonSamples)
#'
#' @md
# @importFrom CoreGx .errorMsg
#' @export
setMethod('findCommonSamples', signature(object='CohortList'),
    function(object)
{
    if (is.null(colnames(object[[1]])))
        stop(.errorMsg(.context(), 'There are no colnames for the first item ',
            'in the object list-like. This method required column names!'))

    Reduce(intersect, lapply(object, colnames))
})