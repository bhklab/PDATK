#' Draw a plot showing the network graph of the significant cluster edges
#' 
#'
#' @param clusterEdges A \code{data.table} containing the statistics for significant
#'    inter-cohort cluster comparisons, as returned by `compareClusters`.
#' @param seed
#' @param savePath
#' @param fileName
#'
#' @import igraph
#' @export
plotClusterNetwork <- function(clusterEdges, seed=NULL, savePath, fileName, ...) {
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
  colours <-  c("blue", brewer.pal(n=8, name="Dark2")[c(4, 5)])
  
  plotNetwork <- call('.plotNetwork', ugraph, metaClusters, coords, colours)
  
  # Plot the network graph
  plot <- base2grob(as.expression(plotNetwork))
  
  if (!missing(savePath) && !missing(fileName))
    ggsave(file.path(savePath, fileName), plot)
  
  grid.draw(plot)
}


#'
#'
#'
#'
#'
.plotNetwork <- function(ugraph, metaClusters, coords, colours) {
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
         legend=c("Basal","Exocrine","Classical"), 
         pt.cex=1.8, pch=21, pt.bg=colours, bty='n', col=colours)
}


#' 
#'
#' @param expressionData
#' @param survivalData
#'
#' @importFrom 
#' @export
extractSurvivalData <- function(expressionData, survivalData) {
  sampleLengths <- c(vapply(expressionData, nrow, FUN.VALUE=numeric(1)),
                     'yang'=nrow(survivalData))
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

#'
#'
#'
#'
#'
#'
extractMetaClusters <- function(clusters, survivalData) {
  metaClusters <- lapply(clusters, `[[`, "metaClasses")
  sampleNames <- lapply(metaClusters, names)
  sampleLengths <- lapply(sampleNames, length)
  cohorts <- unlist(mapply(rep, x=names(clusters), times=sampleLengths, SIMPLIFY=FALSE))
  mClust <- data.frame("cohorts"=cohorts, 
                       "sampleNames"=unlist(sampleNames), 
                       "metaClusters"=unlist(metaClusters))
  mClust[!duplicated(mClust$sampleNames), ]
}

#'
#'
#' @data1 A \code{data.frame} 
#' @data1 A \code{data.frame}
#' @param sharedColumn
#'
#' @export
mergeOn <- function(data1, data2, sharedColumn) {
  df <- merge(data1, data2, by=c("cohorts", sharedColumn))
  return(df)
}

#'
#'
#'
#'
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

#'
#'
#'
#' 
#'
fitProportionalHazardsModel <- function(metaClusterSurvAnnot) {
  coxph(Surv(OS, OSstatus == 1) ~ metaClusters + strata(cohorts), 
        data=metaClusterSurvAnnot)
}

#'
#'
#'
#'
#'
fitSurvivalCurves <- function(metaClusterSurvAnnot) {
  fit <- survfit(Surv(OS, OSstatus == 1) ~ metaClusters, 
          data=metaClusterSurvAnnot)
  survDiff <- survdiff(Surv(OS, OSstatus == 1) ~ metaClusters, 
                 data=metaClusterSurvAnnot)
  pval <- 1 - pchisq(survDiff$chisq, length(survDiff$n) - 1)
  return(list("fit"=fit, "pval"=pval))
}

#'
#'
#'
#'
#'
plotSurvivalCurves <- function(survivalCurves, title="", plot=TRUE, saveDir, fileName) {
  plotFunction <- as.expression(call(".plotSurvivalCurve", survivalCurves, title))
  
  grob <- base2grob(plotFunction)

  if(!missing(saveDir) && !missing(fileName))
    ggsave(grob, file=file.path(saveDir, fileName))
  
  if (plot) {
    grid.draw(grob)
  } else {
    return(grob)
  }
}

#'
#'
#'
#'
.plotSurvivalCurve <- function(survivalCurves, title) {
  plot(survivalCurves$fit, col=c("tomato", "green", "purple"), lwd=2,
       xlab="Days", ylab="Survival Probability")
  legend("bottomleft", paste0("P = ", round(survivalCurves$pval, 2)), bty="n")
  title(title, line=0.2)
}

#'
#'
#'
#'
#'
plotCohortwiseSurvivalCurves <- function(metaClusterSurvAnnot, plot=TRUE, saveDir, fileName) {
  DT <- as.data.table(metaClusterSurvAnnot, keep.rownames=TRUE)
  perCohort <- split(DT, by='cohorts')
  cohortCurves <- lapply(perCohort, fitSurvivalCurves)
  plotGrobs <- mapply(plotSurvivalCurves, cohortCurves, title=names(perCohort), 
                      MoreArgs=list(plot=FALSE), SIMPLIFY=FALSE)
  
  grob <- grid.arrange(grobs=plotGrobs, ncol=ceiling(sqrt(length(plotGrobs))))
  
  if(!missing(saveDir) && !missing(fileName))
    ggsave(grob, file=file.path(saveDir, fileName))
  
  if (plot) {
    grid.draw(grob)
  } else {
    return(grob)
  }
}