#' Run Gene Set Enrichment Analysis
#'
#' @examples
#' \donttest{
#'   data(samplePCOSPmodel)
#'   geneSet <- msigdbr()
#'   GSEAresults <- runGSEA(samplePCOSPmodel, geneSet)
#' }
#'
#' @param object An `S4` object to conduct Gene Set Enrichment Analysis (GSEA)
#'   with.
#' @param geneSet An object representing a gene set, such as a `data.frame`.
#' @param ... Allow additional parameters to be defined for this generic.
#'
#' @return A `data.frame` containing the significantly enriched gene sets.
#'
#' @md
#' @export
setGeneric('runGSEA', function(object, geneSet, ...)
    standardGeneric('runGSEA'))
#'
#' Run Gene Set Enrichment Analysis On A `PCOSP` Model Object.
#'
#' @param object A `PCOSP` model which has been trained with `trainModel`.
#' @param geneSet A `data.frame` with two columns, the first being the name
#'   of the gene and the second the gene set. The gene names must match the
#'   rownames of `object`. Additional columns will be dropped.
#' @param numModels The number of models to use when selecting the top features
#'   from the PCOSP model in `object`. If missing will default to the top 10%
#'   of models.
#' @param ... Force subsequent parameters to be named. Not used.
#' @param adjMethod An optional parameter specifying the multiple testing
#'   correction to use in [`piano::runGSAhyper`]. This parameter must be named.
#' @param allResults Return the full results from [`piano::runGSAhyper`] instead
#'   of a `data.frame` of significant results? Default is FALSE. This parameter
#'   must be named.
#'
#' @return A `data.table` containing the significantly enriched gene sets.
#'
#' @md
#' @importFrom piano runGSAhyper loadGSC
#' @export
setMethod('runGSEA', signature(object='PCOSP', geneSet='data.frame'),
    function(object, geneSet, numModels, ..., adjMethod='fdr', allResults=FALSE)
{
    if (missing(numModels))
        numModels <- floor(length(models(object))*0.1)
    # get the features of interest
    topFeatures <- getTopFeatures(object, numModels)

    if (!any(rownames(object) %in% geneSet[[1]]))
        stop(.errorMsg(.context(), "No gene names in `geneSet` are in ",
            "the rownames of `object`"))

    # provide a background of uninteresting features
    refFeatures <- setdiff(rownames(object), topFeatures)

    if (ncol(geneSet) > 2) {
        # NOTE:: seq_len does not work for column subsetting in data.table
        GSC <- loadGSC(file=geneSet[, 1:2],
            addInfo=geneSet[, 3:ncol(geneSet)])
    } else {
        GSC <- loadGSC(file=geneSet)
    }

    GSEAresults <- runGSAhyper(genes=topFeatures, gsc=GSC, adjMethod=adjMethod)

    gseaDT <- as.data.table(GSEAresults$resTab, keep.rownames='gene_set')

    return(if (allResults) GSEAresults else gseaDT[`Adjusted p-value` < 0.05, ])
})