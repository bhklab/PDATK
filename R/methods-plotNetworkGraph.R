#' A Generic for Plotting a Network Graph From an S4 Object
#'
#' @param object An `S4` object with a valid plotNetworkGraph method set.
#' @param ... Allow additional arguments to be defined for this generic.
#'
#' @return A network plot, either as an object or via side effects.
#'
#' @md
#' @export
setGeneric('plotNetworkGraph',
    function(object, ...) standardGeneric('plotNetworkGraph'))

#' Plot a Network Graph for a Classified NCSModel Object
#'
#' Visualize metaclusters predicted using network community search on the
#'   consensus clustering results for a MultiAssayExperiment of patient
#'   cohorts.
#'
#' @param object A classified `NCSModel` object, as returned by the
#'   `predictClasses` method.
#' @param ... Not used, force subsequent arguments to be named.
#' @param palette `character` A valid pallete for use in `RColourBrewer::brewer.pal`
#' @param clusterLabels A `character` vector of names for the metaclusters.
#'   Defaults to the cluster number.
#'
#' @return A `ggplot` object containing the network graph, showing the
#'   relative edge distances between each cluster in each cohort along with
#'   the predicted metacluster label.
#'
#' @md
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplotify as.ggplot base2grob
#' @export
setMethod('plotNetworkGraph', signature(object='NCSModel'),
    function(object, ..., palette="Set1", clusterLabels)
{
    funContext <- .context(1)

    ## TODO:: Add state labels to S4Model object to determine if it has
    #  been trained, classified, etc.

    ## TODO:: Improve
    if (!('graphData' %in% names(models(object))))
        stop(.errorMsg(funContext, 'There is not graph data in the models ',
            'slot of object. Please ensure you have classified the model ',
            'using predictClasses before attempting to plot the metacluster' ,
            'network grpah!'))

    # -- Extract the relevant graph data
    graphData <- models(object)$graphData
    ugraph <- graphData$ugraph
    coords <- graphData$coords
    metaclusters <- graphData$metaclusters
    colours <- brewer.pal(n=8, palette)

    # -- Plot using igraph::plot.igraph
    if (missing(clusterLabels)) {
        clusterLabels <- as.character(sort(unique(membership(metaclusters))))
    }
    plotNetwork <- call('.plotNetwork', ugraph=ugraph,
        metaClusters=metaclusters, coords=coords, colours=colours,
        clusterLabels=clusterLabels)

    # -- Coerce to ggplot and return
    networkGgplot <- as.ggplot(base2grob(as.expression(plotNetwork)))
    return(networkGgplot)
})

# A helper function for plotting network graphs
#' @export
#' @importFrom igraph plot.igraph membership layout_with_dh
#' @keywords internal
.plotNetwork <- function(ugraph, metaClusters, coords, colours, clusterLabels) {
    igraph::plot.igraph(ugraph,
         vertex.color=colours[membership(metaClusters)],
         vertex.shape="sphere",
         vertex.size=6,
         edge.arrow.size=0.5,
         vertex.label.cex=0.8,
         vertex.label.dist=2,
         edge.curved=0.1,
         vertex.color=colours[membership(metaClusters)],
         edge.arrow.size=0.4,
         layout=layout_with_dh(ugraph),
         layout=coords)
    legend('topleft',
           legend=clusterLabels,
           pt.cex=1.8, pch=21, pt.bg=colours, bty='n',
           col=colours)
}