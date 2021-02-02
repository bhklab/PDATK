#' Find the common genes in an `S4` object.
#'
#' @param object An `S4` object to find common genes for.
#' @param ... Allow new parameters to be defined for this generic.
#'
#' @return A `character` vector of common gene names.
#'
#' @md
#' @export
setGeneric('findCommonGenes', function(object, ...)
    standardGeneric('findCommonGenes'))
#'
#' Intersect Gene Names for All `SurvivalExperiments` in a `CohortList`
#'
#' @param object A `CohortList` of `SurvivalExperiment`s to find common genes
#'   between.
#'
#' @return A `character` vector of genes common to all `SurvivalExperiment`s
#'   in the `CohortList`.
#'
#' @md
#' @export
setMethod('findCommonGenes', signature='CohortList', function(object)
{
    Reduce(intersect, lapply(object, rownames))
})

