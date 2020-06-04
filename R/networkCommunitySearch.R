#' Draw a plot showing the network graph of the significant cluster edges
#'
#'
#' @param clusterEdges A \code{data.table} containing the statistics for significant
#'    inter-cohort cluster comparisons, as returned by the `compareClusters`
#'    function in this package.
#' @param seed The seed to set for the graph community search.
#' @param palette The name of an RColorBewer palette to apply to the plot,
#'     defaults to 'Set1'.
#' @param clusterLabs An optional \code{character} vector of cluster names
#'     corresponding to the numbered metaclusterings.
#' @param saveDir An optional \code{character} vector with the path to
#'    the directory in which the plot should be saved.
#' @param fileName An optional \code{character} vector specifying the name
#'    and extension of the file. Only used if `saveDir` is also specified.
#'    Saving is done via the `ggsave` function from `ggplot2`.
#'
#' @importFrom igraph graph_from_edgelist layout_with_fr as.undirected
#'     fastgreedy.community membership
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggsave
#' @importFrom grid grid.draw
#' @importFrom ggplotify base2grob as.ggplot
#' @export
plotClusterNetwork <- function(clusterEdges, seed=NULL, palette="Set1",
                               clusterLabels, saveDir, fileName, ...) {
  
 # A helper function for plotting network graphs
  .plotNetwork <- function(ugraph, metaClusters, coords, colours, clusterLabels) {
      plot(ugraph,
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
          pt.cex=1.8, pch=21, pt.bg=colours, bty='n', col=colours)
  }
  
  # Set seed for reproducible results
  if (!is.null(seed)) set.seed(seed)

  # Prepare the edge labels
  cohort1 <- vapply(strsplit(clusterEdges$comparison, '-'), `[`,i= 1, character(1))
  cohort2 <- vapply(strsplit(clusterEdges$comparison, '-'), `[`, i=2, character(1))
  edges <- cbind(
    'cluster1'=unlist(mapply(paste0, cohort1, '-', clusterEdges$c1Clust, SIMPLIFY=FALSE)),
    'cluster2'=unlist(mapply(paste0, cohort2, '-', clusterEdges$c2Clust, SIMPLIFY=FALSE))
  )

  # Format the network graph
  graph <- graph_from_edgelist(edges)
  coords <- layout_with_fr(graph)
  ugraph <- as.undirected(graph)
  metaClusters <- fastgreedy.community(ugraph, weights=clusterEdges$threshold)
  colours <-  brewer.pal(n=8, name=palette)

  if (missing(clusterLabels)) {
    clusterLabels <- as.character(sort(unique(membership(metaClusters))))
  }
  plotNetwork <- call('.plotNetwork', ugraph, metaClusters, coords, colours,
                      clusterLabels)

  # Plot the network graph
  plot <- as.ggplot(base2grob(as.expression(plotNetwork)))

  if(!missing(saveDir) && !missing(fileName)) {
    ggsave(file.path(saveDir, fileName), plot)
    message(paste0("Saved to ", file.path(saveDir, fileName)))
  }
  return(plot)
}





#' Get the survival data from a list of expression
#'
#' @param expressionData A \code{list} of per cohort sample by gene expression
#'    matrixes annotated witht he columns 'OS' for overall survival time and
#'    'OS_Status' for survival event.
#' @param survivalData A \code{data.frame} containing additional sample
#'    survival data.
#' @param survivalDataLabel A \code{character} vector with the name to use
#'    when labelling the additional survival data. Defaults to 'additional'
#'    if not specified. This label will identify the additional survival data
#'    in the returned \code{data.frame}
#'
#' @return A \code{data.frame} with columns cohorts, sampleNames, OS and OSstatus
#'   containing the per sample survival information.
#'
#' @export
extractSurvivalData <- function(expressionData, survivalData, survivalDataLabel="additonal") {
  sampleLengths <- c(vapply(expressionData, nrow, FUN.VALUE=numeric(1)),
                     nrow(survivalData))
  names(sampleLengths)[length(expressionData) + 1 ] <- survivalDataLabel
  cohorts <- unlist(mapply(rep, x=names(sampleLengths), times=sampleLengths, SIMPLIFY=FALSE))
  samples <- c(unlist(lapply(expressionData, rownames)),
               as.character(survivalData$ID))
  os <- c(unlist(lapply(expressionData, `[[`, "OS")),
          as.numeric(survivalData$OS))
  osStatus <- c(unlist(lapply(expressionData, `[[`, "OS_Status")),
                as.numeric(survivalData$OS_Status))
  return(data.frame("cohorts"=cohorts,
                    "sampleNames"=samples,
                    "OS"=as.numeric(os),
                    "OSstatus"=as.numeric(osStatus)))
}

#' Get sample metaclusters
#'
#' @param clusters A \code{list} of per cohort consensus clustering results,,
#'   as returned by the `conClustAllCohorts` function in this package.
#'
#' @export
extractMetaClusters <- function(clusters) {
  metaClusters <- lapply(clusters, `[[`, "metaClasses")
  sampleNames <- lapply(metaClusters, names)
  sampleLengths <- lapply(sampleNames, length)
  cohorts <- unlist(mapply(rep, x=names(clusters), times=sampleLengths, SIMPLIFY=FALSE))
  mClust <- data.frame("cohorts"=cohorts,
                       "sampleNames"=unlist(sampleNames),
                       "metaClusters"=unlist(metaClusters))
  mClust[!duplicated(mClust$sampleNames), ]
}

#' Merge two data.frames or data.tables on cohorts and some other shared column
#'
#' @data1 A \code{data.frame} or \code{data.table} with the columns cohorts
#'   as well as `sharedColumn`
#' @data1 A \code{data.frame} or \code{data.table} with the columns cohorts
#'   as well as `sharedColumn`
#' @param sharedColumn A \code{character} vector containing the name
#'   of a shared column to merge on, in addition to cohorts.
#'
#' @export
mergeSurvOn <- function(data1, data2, sharedColumn) {
  if (is(data2, 'data.table') && is(data2, 'data.table')) {
    df <- as.data.frame(data.table::merge(data1, data2,
                                          by=c("cohorts", sharedColumn)))
  } else {
    df <- base::merge(data1, data2, by=c("cohorts", sharedColumn))
  }
  return(df)
}

#' Add cluster labels to the metaClusters column of a survival data.frame
#'
#' @param metaClusterSurvival A \code{data.frame} containing per sample
#'   metaclusters and survival data, as returned by the `mergeSurvOn` function
#'   in this package.
#' @param clusterLabels A \code{character} vector with labels for each metacluster
#'   in `metaClusterSurvival`. Matches by postion, with metacluster 1 being labelled
#'   with the first item in `clusterLabels` and so on.
#'
#' @return The `metaClusterSurvival` \code{data.frame} with the 'metaClusters'
#'   column labelled.
#'
#' @export
annotateMetaClusters <- function(metaClusterSurvival, clusterLabels) {
  DT <- as.data.table(metaClusterSurvival, keep.rownames=TRUE)
  DT[, metaClusters := as.character(metaClusters)]
  i <- 1
  for (val in na.omit(unique(DT$metaClusters))) {
    DT[metaClusters == val, metaClusters := clusterLabels[i]]
    i <- i + 1
  }
  return(as.data.frame(DT[, -'rn'], row.names=DT$rn))
}

#' Fit a proportional hazards model to survival data using the `survial::coxph`
#'    function.
#'
#'
#' @param metaClusterSurvAnnot A \code{data.frame} containg per sample survival data
#'   and labelled metaclusters, as returned by the `annotateMetaClusters` function
#'   in this package.
#'
#' @return A \code{coxph} object containing the fitted proportional hazards
#'    model.
#'
#' @import survival
#' @export
fitProportionalHazardsModel <- function(metaClusterSurvAnnot) {
  coxph(Surv(OS, OSstatus) ~ metaClusters + strata(cohorts),
        data=metaClusterSurvAnnot)
}

#' Fit Kaplan-Meyer survival curves by metacluster
#'
#' @param metaClusterSurvAnnot A \code{data.frame} or other matrix like object
#'    containing the metacluster labels in the 'metaClusters' column,
#'    overall survival time in the 'OS' column and survival status in the
#'    'OSstatus' column.
#'
#' @import survival
#' @export
fitSurvivalCurves <- function(metaClusterSurvAnnot) {
  fit <- survfit(Surv(OS, OSstatus) ~ metaClusters,
          data=metaClusterSurvAnnot)
  survDiff <- survdiff(Surv(OS, OSstatus) ~ metaClusters,
                 data=metaClusterSurvAnnot)
  pval <- 1 - pchisq(survDiff$chisq, length(survDiff$n) - 1)
  return(list("fit"=fit, "pval"=pval))
}

#' Plot the survival curves for each predicted metacluster
#'
#' @param survivalCurves A \code{list} whose first item is the `survfit` object
#'    from curve fitting and second item is the associated pvalue for the
#'    difference between survival curves.
#' @param title A \code{character} vector containing the plot title. Optional
#'    and defaults to ''.
#' @param saveDir An optional \code{character} vector specifying the path
#'    to the directory where the plot should be saved. If excluded, fileName
#'    will not work.
#' @param palette An optional RColorBrewer palette to colour the survival curves,
#'   defaults to 'Set1'.
#' @param inset An optional parameter specifying the how far to inset the legend
#'   in the plot. Useful for adjusting legend position in plot grids. Defualts
#'   to zero.
#' @param showPlot A \code{boolean} indicating wether or not to draw the plot.
#'   For internal use in `plotCohortwiseSurvivalCurves` to prevent plotting
#'   twice.
#' @param fileName An optional \code{character} vector specifying the
#'    name and extension of the file to save the plot it. This is passed to
#'    the `ggplot2::ggsave`.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggsave
#' @importFrom grid grid.draw
#' @importFrom ggplotify base2grob as.ggplot
#' @export
plotSurvivalCurves <- function(survivalCurves, title="", showPlot=TRUE,
                               palette="Set1", inset=0, saveDir, fileName) {
  plotFunction <- as.expression(call(".plotSurvivalCurve", survivalCurves, title,
                                     palette, inset))

  grob <- base2grob(plotFunction)

  if(!missing(saveDir) && !missing(fileName))
    ggsave(grob, file=file.path(saveDir, fileName))

  if (showPlot) {
    grid.draw(grob)
  } else {
    return(grob)
  }
}

#' Draw a survival curve; for use in plotSurvivalCurves to label plots by
#'
#' @param survivalCurves A \code{survfit} object as returned by the
#'     `survival::survfit` function.
#' @param title An optional title for the plot.
#' @param palette An optional RColorBrewer palette to colour the survival curves,
#'   defaults to 'Set1'.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales scientific
#' @keywords internal
#' @export
.plotSurvivalCurve <- function(survivalCurves, title, palette, inset=0) {
  plot(survivalCurves$fit, col=brewer.pal(n=8, palette), lwd=2,
       xlab="Days", ylab="Survival Probability")
  legend("topright", title=paste0("P = ", scientific(survivalCurves$pval, 2)), bty="n",
         fill=brewer.pal(n=8, palette),
         legend=sort(gsub('.*=', '', names(survivalCurves$fit$strata))),
         inset=inset)
  title(title, line=0.2)
}

#' Plot a grid of per cohort survival curves
#'
#' @param metaClusterSurvAnnot A \code{data.frame} containing the per cohort
#'    survival data
#' @param palette An optional RColorBrewer palette to colour the survival curves,
#'   defaults to 'Set1'.
#' @param inset An optional parameter specifying the how far to inset the legend
#'   in the plot. Useful for adjusting legend position in plot grids. Defualts
#'   to zero.
#' @param saveDir An optional \code{character} vector specifying the path
#'    to the directory where the plot should be saved. If excluded, fileName
#'    will not work.
#' @param fileName An optional \code{character} vector specifying the
#'    name and extension of the file to save the plot it. This is passed to
#'    the `ggplot2::ggsave`.
#'
#' @return A \code{ggplot} object
#'
#' @import data.table
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggsave
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplotify base2grob as.ggplot
#' @export
plotCohortwiseSurvCurves <- function(metaClusterSurvAnnot, plot=TRUE, inset=0.2,
                                     saveDir, fileName) {
  DT <- as.data.table(metaClusterSurvAnnot, keep.rownames=TRUE)
  perCohort <- split(DT, by='cohorts')
  cohortCurves <- lapply(perCohort, fitSurvivalCurves)
  plotGrobs <- mapply(plotSurvivalCurves, cohortCurves, title=names(perCohort),
                      MoreArgs=list(showPlot=FALSE, inset=inset), SIMPLIFY=FALSE)

  grob <- grid.arrange(grobs=plotGrobs, ncol=ceiling(sqrt(length(plotGrobs))))

  if(!missing(saveDir) && !missing(fileName))
    ggsave(grob, file=file.path(saveDir, fileName))

  if (plot) {
    grid.draw(grob)
  } else {
    return(grob)
  }
}
