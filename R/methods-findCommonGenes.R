#' Intersect Gene Names for All `SurvivalExperiments` in a `CohortList`
#'
# @examples
# data(sampleCohorts)
# commonGenes <- findCommonGenes(sampleCohorts)
#'
#' @param object An `S4` object to find common genes for.
#'
#' @return A `character` vector of common gene names.
#'
#' @md
#' @export
setGeneric('findCommonGenes', function(object, ...)
    standardGeneric('findCommonGenes'))
#' @export
setMethod('findCommonGenes', signature='CohortList', function(object)
{
    Reduce(intersect, lapply(object, rownames))
})

