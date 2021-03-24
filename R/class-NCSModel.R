#' S4Model Class for Metaclustering Using Network Community Search
#' 
#' @inherit S4Model
#' 
#' @md
#' @include class-S4Model.R
#' @export
.NCSModel <- setClass('NCSModel',
    contains='S4Model'
)

#' Constructor for a NetworkCommunitySearchModel (NCSModel)
#' 
#' Use `igraph::fastgreedy.community` cross-cohort metaclusters for
#'   a `ConsensusMetaclusteringModel` which has been validated with 
#'   `validateModel`.
#' 
#' @param model A validated `ConsensusMetaclusteringModel` object.
#' 
#' @md
#' @importFrom CoreGx .errorMsg
#' @importFrom igraph graph_from_edgelist layout_with_fr as.undirected
#'     fastgreedy.community
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder setnames
#' @importFrom MultiAssayExperiment experiments experiments<-
#' @importFrom S4Vectors endoapply mendoapply merge
#' @aliases NCSModel
#' @export
NetworkCommunitySearchModel <- function(model) 
{
    funContext <- .context(1)

    if (!is(model, 'ConsensusMetaclusteringModel'))
        stop(errorMsg(funContext, ' The model argument is a ', class(model),
            ' but must be a ConsensusClusteringModel!'))

    cohortMAE <- c(trainData(model), validationData(model)$experiments)
    corTestStats <- validationStats(model)[metric == 'cor.test', ]
    clusterReproStats <- validationStats(model)[metric == 'clusterRepro', ]

    corThresholdDT <- corTestStats[, .SD[which.max(estimate), ], 
        by=.(comparison, centroid_K)]

    if (unique(corThresholdDT[, length(unique(assay_K)), 
        by=.(comparison, centroid_K)]$V1) != 1) 
    {
        stop(.errorMsg(funContext, ' More than one assay cluster has been',
            'selected per centroid cluster. Something has gone wrong, please',
            'check the validationStats of your ', class(model)[1] ,' object.'))
    }

    # Prepare correlation threshold DT for join
    setnames(corThresholdDT, 'estimate', 'cor_threshold', skip_absent=TRUE)
    corThresholdDT[, c('metric', 'p_value', 'assay_N') := NULL]
    
    # Prepare cluster reproducibility stats DT for join
    setnames(clusterReproStats, 'estimate', 'ingroup_proportion', 
        skip_absent=TRUE)
    clusterReproStats[, c('metric', 'assay_K') := NULL]

    # Join to filter for significantly reproducible clusters
    clusterEdgeDT <- merge.data.table(corThresholdDT, clusterReproStats, 
        by=intersect(colnames(corThresholdDT), colnames(clusterReproStats)))

    .NCSModel(
        trainData=cohortMAE,
        modelParams=SimpleList(),
        models=SimpleList(list(
            networkEdges=clusterEdgeDT
        )),
        validationStats=data.table(),
        validationData=SimpleList(),
        metadata=list(ConsensusMetaclusteringModel_modelParams=modelParams(model))
    )

}
#' @export
NCSModel <- NetworkCommunitySearchModel
