#' Generic for Predicting Metaclusters Based on A Trained S4Model
#' 
#' @param object An `S4` object to predict classes for.
#' @param model An `S4` object representing a classification or regression model.
#'   Should inherit from `S4Model`.
#' 
#' @return The `object` argument with metacluster predictions in the column
#'   metadata.
#' 
#' @md
#' @export
setGeneric('predictMetaclasses', function(object, model, ...) 
    standardGeneric('predictMetaclasses'))

#' predictMetaclasses Method for ConsensusMetaclusteringModel
#' 
#' Use `igraph::fastgreedy.community` cross-cohort metaclusters for
#'   a `ConsensusMetaclusteringModel` which has been validated with 
#'   `validateModel`.
#' 
#' @param object Not used, should be omitted.
#' @param model A validated `ConsensusMetaclusteringModel` object.
#' @param alpha A `float` specifying the significance level for cluster
#'   reproducibility. Default is 0.05.
#' @param minRepro A `float` specifying the minimum in-group proportion for
#'   a cluster to be included in the metacluster labels. Default is 0.5.
#' @param minCor A `float` specifying the minimum correlation between a 
#'   centroid and assay cluster to be included in the metacluster labels. 
#'   Default is 0.0.
#' 
#' @md
#' @importFrom CoreGx .errorMsg
#' @importFrom igraph graph_from_edgelist layout_with_fr as.undirected
#'     fastgreedy.community
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder
#' @export
setMethod('predictMetaclasses', signature(object='missing', 
    model='ConsensusMetaclusteringModel'), function(object, model, alpha=0.05, 
        minRepro=0.5, minCor=0.0) 
{
    funContext <- .context()

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

    # Prepare corrleation threshold DT for join
    setnames(corThresholdDT, 'estimate', 'cor_threshold', skip_absent=TRUE)
    corThresholdDT[, c('metric', 'p_value', 'assay_N') := NULL]
    
    # Prepare cluster reproducibility stats DT for join
    setnames(clusterReproStats, 'estimate', 'ingroup_proportion', 
        skip_absent=TRUE)
    clusterReproStats[, c('metric', 'assay_K') := NULL]

    # Join to filter for significantly reproducible clusters
    clusterEdgeDT <- merge.data.table(corThresholdDT, clusterReproStats, 
        by=intersect(colnames(corThresholdDT), colnames(clusterReproStats)))

    # Filter to signficant cluster edges based on input parameters
    signifClusterEdgeDT <- clusterEdgeDT[
        p_value <= alpha &
        ingroup_proportion >= minRepro &
        cor_threshold >= minCor,
    ]

    # Extract the edge labels for network community search
    clusterEdges <- as.matrix(signifClusterEdgeDT[, 
        .(centroid_cluster=paste0(centroid_cohort, '-', centroid_K), 
            assay_cluster=paste0(assay_cohort, '-', assay_K))
    ])

    # Construct network graph and identify metaclusters using network community
    #   search
    graph <- graph_from_edgelist(clusterEdges)
    coords <- layout_with_fr(graph)
    ugraph <- as.undirected(graph)
    metaClusters <- fastgreedy.community(ugraph, 
        weights=signifClusterEdgeDT$cor_threshold)

})
