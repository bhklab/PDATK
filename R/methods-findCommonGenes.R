#' Find the common genes in an `S4` object.
#'
#' @param object An `S4` object to find common genes for.
#' @param ... Allow new parameters to be defined for this generic.
#'
#' @return A `character` vector of common gene names.
#'
#' @examples
#' data(sampleCohortList)
#' commonGenes <- findCommonGenes(sampleCohortList)
#' head(commonGenes)
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
#' @examples
#' data(sampleCohortList)
#' commonGenes <- findCommonGenes(sampleCohortList)
#' head(commonGenes)
#'
#' @md
#' @export
setMethod('findCommonGenes', signature(object='CohortList'), function(object)
{
    Reduce(intersect, lapply(object, rownames))
})

#' Intersect Gene Names for All `experiments` in a `MultiAssayExperiment`
#'
#' @param object A `MultiAssayExperiment` where rownames represent genes.
#'
#' @return A `character` vector of genes common to all `experiments`
#'   in the `MutliAssayExperiment`.
#'
#' @md
#' @export
setMethod('findCommonGenes', signature(object='MultiAssayExperiment'),
    function(object)
{
    Reduce(intersect, rownames(object))
})
