#' Find Common Samples in a List-like S4 Object where The Columns of Each Item
#'   Represent Samples
#'
#' @param object A `S4` object, where the columns of each element represent
#'   samples.
#'
#' @md
#' @export
setGeneric('findCommonSamples',
    function(object, ...) standardGeneric('findCommonSamples'))
#'
#' @param object A `CohortList` for which we want to find common samples
#'   between all `SurvivalExperiment` objects.
#'
#' @export
setMethod('findCommonSamples', signature('CohortList'),
    function(object)
{
    if (is.null(colnames(object[[1]])))
        stop(.errorMsg(.context(), 'There are no colnames for the first item ',
            'in the object list-like. This method required column names!'))

    Reduce(intersect, lapply(object, colnames))
})