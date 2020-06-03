#' Calculate the weighted MAD across all cohorts per gene
#'
#' @param cohortMADrankings A \code{list} of \code{data.frames} containing
#'     the genes, madValuesa nd ranking for each cohort. As returned by
#'     the `rowMADdfs` function.
#'
#' @return A \code{numeric} vector
#'
#' @importFrom matrixStats weightedMad
#' @export
calcGenewiseWeightedMADs <- function(cohortMADrankings) {

  cohortMADvals <- lapply(cohortMADrankings, `[[`, "madValues")
  cohortNcols <- vapply(cohortMADrankings,
                        function(cohort) attributes(cohort)$datasetCols,
                        FUN.VALUE=numeric(1))

  perGeneMADs <- lapply(seq_along(rownames(cohortMADrankings[[1]])),
                        function(i, MADs) vapply(MADs, `[`, i=i,
                                                 FUN.VALUE=numeric(1)),
                        MADs=cohortMADvals)

  vapply(perGeneMADs,
         function(x, cohortNcols, na.rm) weightedMad(x, cohortNcols, na.rm=na.rm),
         cohortNcols, na.rm=TRUE,
         FUN.VALUE=numeric(1))
}

#' Calculate the weighted MAD across all cohorts per gene and
#'     return a `data.frame` of the results.
#'
#' @param cohortMADrankings A \code{list} of \code{data.frames} containing
#'     the genes, madValuesa nd ranking for each cohort. As returned by
#'     the `rankAllCohortGenesByMAD` function.`
#'
#' @return A \code{data.frame} with columns genes, weightedMADs and rankings
#'     ordered by decreasing weightedMAD values.
#'
#' @importFrom dplyr dense_rank
#' @export
calcGenewiseWeightedMADdf <- function(cohortMADrankings) {
  weightedMADs <- calcGenewiseWeightedMADs(cohortMADrankings)
  data <- data.frame(
    "genes"=cohortMADrankings[[1]]$genes,
    "weightedMADs"=weightedMADs,
    "rankings"=dense_rank(-weightedMADs)
  )
  data[order(data$weightedMADs, decreasing=TRUE), ]
}

#' Extract the top n genes from each cohort
#'
#' @param cohortMADrankings A \code{list} of \code{data.frames} containing
#'     the genes, madValuesa nd ranking for each cohort. As returned by
#'     the `rowMADdfs` function.
#' @param n A \code{numeric} vector with the integer number of top genes
#'     to return
#'
#' @return The `cohortMADrankings` \code{data.frame} for the top n genes
#'     in each cohort.
#'
#' @export
getTopGenes <- function(cohortMADrankings, n) {
  topN <- seq_len(n)
  MADs <- lapply(cohortMADrankings, `[[`, "madValues")
  topRankIdxs <- lapply(MADs,
                        function(MAD, topN) order(MAD, decreasing=TRUE)[topN],
                        topN=topN)
  topGeneDFs <- lapply(seq_along(cohortMADrankings),
                       function(i, topRankIdx, cohortMADrankings)
                         cohortMADrankings[[i]][topRankIdx[[i]], ],
                       topRankIdxs,
                       cohortMADrankings)
  structure(topGeneDFs,
            .Names=names(cohortMADrankings))
}

#' Get a vector of the unique top n genes for each cohort appended to the
#'     top m genes from the
#'
#' @param cohortMADrankings A \code{list} of \code{data.frames} containing
#'     the genes, madValuesa nd ranking for each cohort. As returned by
#'     the `rowMADdfs` function.
#' @param geneWiseMADdf A \code{data.frame} containing the ordered gene
#'     rankings across all cohorts, as returned by `calcGeneWiseWeightedMADdf`.
#' @param n A \code{numeric} vector with the integer number of top genes
#'     to return from each cohort.
#' @param m A \code{numeric} vector with the integer number of top genes
#'     to return from the meta-cohort rankings.
#'
#' @return A \code{character} vector of the unique top n ranked genes
#'     per cohort and top m genes in the meta-cohort ranking.
#'
#' @export
getMetaGenes <- function(cohortMADrankings, geneWiseWeigthedMADdf, n, m) {
    top20genesPerCohort <- getTopGenes(cohortMADrankings, n)
    topNgenesPerCohort <- Reduce(c, lapply(top20genesPerCohort, `[`, i=TRUE, j="genes"))
    topMgenesMeta <- geneWiseWeigthedMADdf$genes[seq_len(m)]
    unique(c(topNgenesPerCohort, topMgenesMeta))
}

#' Consensus cluster on gene expression for a cohort
#'
#' @param cohort A \code{matrix} of expression data for a cohort, with rows
#'     as genes and columns as samples. As in a list item from
#'     `preprocessCohorts`.
#' @param maxK A \code{numeric} vector containing the integer number max
#'     number of clusters to use. See
#'    `?ConsensusClusterPlus` for more information.
#' @param distance A \code{character} vector containing the name of
#'    the distance metric to use in `ConsensusClusterPlus`. See
#'    `?ConsensusClusterPlus` for more information.
#' @param method A \code{character} vector containing the name of
#'    the clustering method to use in `ConsensusClusterPlus`. See
#'    `?ConsensusClusterPlus` for more information.
#' @param reps A \code{numeric} vector containing the integer number of times
#'    to repeat clustering. Defaults to 1000.
#' @param suppressPlots A \code{logical} vector indicating whether to suppress
#'    the plots fomr ConsensusClusterPlus. Defaults to TRUE.
#' @param plotTo A \code{character} vector passed to the plot argument of
#'    `ConsensusClusterPlus`. Defaults to NULL.
#' @param seed A \code{numeric} vector containing an integer random seed to
#'    set for sampling in the `ConsensusClustPlus` function.
#' @param ... Fallthrough arguments to `ConsensusClusterPlus`. If specficied
#'     the function defaults are overridden.
#'
#' @return A \code{list} of consensus clustering results with items optimumK,
#'     clusterTable, centroidClusters and classes.
#'
#' @section Warning: This function uses random numbers, remember to set.seed()
#'     before running it to ensure reproducible results!
#'
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @export
conClustCohort <- function(cohort, maxK, distance, method, reps=1000,
                                   suppressPlots=TRUE, plotTo=NULL, seed=NULL, ...) {

  if (suppressPlots) pdf(NULL)

  # Compute cluster ensemble
  if (!missing(...)) {
    clusters <- ConsensusClusterPlus(cohort, maxK=maxK,
                                     distance=distance, clusterAlg=method,
                                     ...)
  } else {
    clusters <- ConsensusClusterPlus(cohort, maxK=maxK,
                                     distance=distance, clusterAlg=method,
                                     reps=reps,
                                     pItem=0.8,
                                     pFeature=1,
                                     innerLinkage="complete",
                                     corUse="pairwise.complete.obs",
                                     plot=plotTo,
                                     seed=seed)
  }

  if (suppressPlots) dev.off()

  # Find optimal K value
  Kvals <- seq(2, maxK)
  consensusMatrixes <- lapply(clusters[Kvals], `[[`, "consensusMatrix")
  lowerTris <- lapply(consensusMatrixes, lower.tri)
  lowerTriConMatrix <- mapply(`[`, consensusMatrixes, lowerTris, SIMPLIFY=FALSE)
  ecdfs <- lapply(lowerTriConMatrix, ecdf)
  subinterval <- c(0.1, 0.9)
  propAmbigClust <- vapply(ecdfs,
                           function(ecdf, subinterval)
                               ecdf(subinterval[2]) - ecdf(subinterval[1]),
                           subinterval=subinterval,
                           FUN.VALUE=numeric(1))
  optimalK <- Kvals[which.min(propAmbigClust)]

  # Generate return results
  clusterTable <- table(clusters[[optimalK]]$consensusClass)

  clusterClasses <- clusters[[optimalK]]$consensusClass
  centroidClusters <- sapply(unique(clusterClasses),
                            .calcClusterCentroid,
                            t(cohort), clusterClasses)

  list(
    "optimalK"=optimalK,
    "clusterTable"=clusterTable,
    "centroidClusters"=centroidClusters,
    "classes"=clusterClasses
  )
}

#' Consensus cluster a list of normalized expression matixes
#'
#' @param cohortL A \code{list} of expression data matrixes, with rows
#'     as genes and columns as samples. As retruned by the `preprocCohorts`
#'     function.
#' @param maxK A \code{numeric} vector containing the integer number max
#'     number of clusters to use. See
#'    `?ConsensusClusterPlus` for more information.
#' @param distance A \code{character} vector containing the name of
#'    the distance metric to use in `ConsensusClusterPlus`. See
#'    `?ConsensusClusterPlus` for more information.
#' @param method A \code{character} vector containing the name of
#'    the clustering method to use in `ConsensusClusterPlus`. See
#'    `?ConsensusClusterPlus` for more information.
#' @param reps A \code{numeric} vector containing the integer number of times
#'    to repeat clustering. Defaults to 1000.
#' @param seed A \code{numeric} vector containing an integer random seed to
#'    set for sampling in the `ConsensusClustPlus` function. Defaults to NULL,
#'    as in no seed.
#' @param nthread An optional argument specifying the number of threads to
#'    parallelize over. If excluded, the default number of threads from the
#'    `BiocParallel` package will be used.
#'
#' @importFrom BiocParallel bplapply
#' @export
conClustAllCohorts <- function(preprocCohorts, maxK=5, distance="pearson",
                               method="hc", nthread, reps=1000, seed=NULL) {
  # Temporily change number of cores to parallelize over
  opts <- options()
  options("mc.cores"=nthread)
  on.exit(options(opts))

  bplapply(preprocCohorts,
           conClustCohort,
           maxK=maxK, distance=distance,
           method=method, seed=seed,
           reps=reps)
}

#' Calculate the centroid from a
#'
#' @param i The cluster index
#' @param cohort The cohort data with genes as columns, rows as samples
#' @param clusterClasses The
#'
#'
#' @importFrom matrixStats colMeans
#' @keywords internal
.calcClusterCentroid <- function(i, cohort, clusterClasses) {
  clusterIdxs <- which(clusterClasses == i)

  if (length(clusterIdxs) == 1) {
    cohort[clusterIdxs, ]
  } else {
    colMeans(cohort[clusterIdxs, ], na.rm=TRUE)
  }
}

#' Subset all cohorts to meta genes, transform
#'
#' @param cohorts A \code{list} of cohort \code{data.frame}s or \code{matrix}es
#' @param metaGenes A \code{character} vector of meta genes to subset on, as
#'    returned by `getMetaGenes`.
#' @param center A \code{logical} vector passed to the `center` argument of
#'    the `base::scale`` function. Defaults to TRUE if excluded.
#' @param scale A \cpde{logical} vector passed to the `scale` argument of
#'    the `base::scale` function. Defaults to FALSE if excluded.
#'
#' @export
preprocCohorts <- function(cohorts, metaGenes, center=TRUE, scale=FALSE) {
  cohortSubsets <- .subsetCohorts(cohorts, genes=metaGenes)
  structure(lapply(cohortSubsets, .scaleGenewise,
                   center=center, scale=scale),
            .Names=names(cohorts))
}

#' Rank genes in a list of expression matrices by MAD
#'
#' @param cohortDataL A \code{list} of cohort expression matrixes.
#' @param nthread A \code{numeric} vector containing the integer number of
#'    threads to parallelize over.
#'
#' @return A \code{list} of `data.frame`s containng the gene names, MAD score
#'    and rankingsfor each cohort in `cohortDataL`
#'
#' @importFrom BiocParallel bplapply
#' @export
rankAllCohortGenesByMAD <- function(cohortDataL, nthread) {
  # Temporily change number of cores to parallelize over
  opts <- options()
  options("mc.cores"=nthread)
  on.exit(options(opts))

  bplapply(cohortsDataL, function(cohort) calcRowMADdf(cohort))
}

#' Convert to matrix, transform, scale row-wise, transform, return
#'
#' @param cohortSubet A \code{data.frame} containing the cohort
#'    data subset to the genes of interest. Rows are features,
#'    columns are samples.
#'
#' @keywords internal
.scaleGenewise <- function(cohortSubset, center, scale) {
  t(scale(t(data.matrix(cohortSubset)),
          center=center, scale=scale)
    )
}

#' Subset all cohorts in a list of
#'
#' @param cohorts A \code{list} of cohorts to subset
#' @param genes A \code{character} vector of genes to subset to, if
#'    excluded all genes are returned.
#' @param columns A \code{character} vector of samples to subset to,
#'    if excluded all samples are returned.
#'
#' @keywords internal
.subsetCohorts <- function(cohorts, genes=TRUE, columns=TRUE) {
  lapply(cohorts, `[`, i=genes, j=columns)
}

#' Find all non-self pair-wise combinations of cohorts
#'
#' @param clusterNames A \code{character} vector of cohort names
#'
#' @return A \code{data.frame} with the index of all non-self pair-wise cohort combinations. Rownames
#'     are the names of the two clusters being compared.
#'
#' @export
findAllCohortPairs <- function(clusterNames) {
    pairs <- expand.grid(x=seq_along(clusterNames), y=seq_along(clusterNames))
    namePairs <- expand.grid(x=clusterNames, y=clusterNames)

    # Remove self comparisons
    allPairs <- pairs[-which(pairs[, 1] == pairs[, 2]), ]
    allNames <- namePairs[-which(namePairs[, 1] == namePairs[, 2]), ]

    # Paste together pair names
    pairNames <- mapply(paste, allNames[, 1], allNames[, 2], MoreArgs=list(sep="-"))

    # Assign names as rownames to pair data.frame
    rownames(allPairs) <- pairNames

    return(allPairs)
}


#' Calculate the MSM threshold between
#'
#' @param cohortPair
#' @param allConClusters
#' @param allProcCohorts
#'
#' @return A \code{data.table} 
#'
#' @import data.table
#' @export
calcMSMthresholds <- function(cohortPair, allConClusters, allProcCohorts) {

  # Get the name of the comparison and a vector of the cohort list indexes
  comparison <- rownames(cohortPair)
  cohortPair <- as.numeric(cohortPair)

  # Find all pair-wise K comparisons between the two consensus clusters
  clust1Ks <- seq_len(allConClusters[[cohortPair[1]]]$optimalK)
  clust2Ks <- seq_len(allConClusters[[cohortPair[2]]]$optimalK)

  # Extract data
  clust1Centroid <- na.omit(allConClusters[[cohortPair[1]]]$centroidClusters)
  clust2Data <- na.omit(allProcCohorts[[cohortPair[2]]])
  clust2Classes <- na.omit(allConClusters[[cohortPair[2]]]$classes)

  # Find common genes
  sharedGenes <- intersect(rownames(clust1Centroid), rownames(clust2Data))

  # Subset to common genes
  clust1Centroid <- clust1Centroid[sharedGenes, ]
  clust2Data <- clust2Data[sharedGenes, ]

  # Calculate the threshold values for each K comparison
  corList <- vector("list", max(clust1Ks) * max(clust2Ks))
  k <- 1
  for (i in clust1Ks) {
    for (j in clust2Ks) {
      classIdxs <- which(clust2Classes == j)
      centroid <- clust1Centroid[, i]
      meanThresh <- mean(unlist(lapply(classIdxs,
                            function(idx, centroid, data) {
                              as.numeric(
                                cor.test(
                                  centroid,
                                  data[, idx],
                                  method="pearson")$estimate
                                )},
                            centroid=centroid,
                            data=clust2Data)))
      corList[[k]] <- data.table('comparison'=comparison,
                                 'cohort1'=cohortPair[1],
                                 'cohort2'=cohortPair[2],
                                 'c1Clust'=i,
                                 'c2Clust'=j,
                                 "threshold"=meanThresh)
      k <- k + 1
    }
  }
  rbindlist(corList, fill=TRUE)
}

#' Calculate the MSM thresholds for each non-self pair-wise comparison of
#'    cohort consensus clustering results.
#'
#' @param cohortPairs
#' @param allConClusters
#' @param allProcCohorts
#' @param nthread
#'
#' @importFrom BiocParallel bplapply
#' @import data.table
#' @export
calcAllMSMthresholds <- function(cohortPairs, allConClusters,
                                 allProcCohorts, nthread) {
  if (!missing(nthread)) {
    opts <- options()
    options(mc.cores=nthread)
    on.exit(options(opts))
  }

  DTlist <- bplapply(seq_len(nrow(cohortPairs)),
           function(idx, pairs, cluster, data) {
             calcMSMthresholds(pairs[idx, ], cluster, data)
           },
           pairs=cohortPairs,
           cluster=allConClusters,
           data=allProcCohorts)
  rbindlist(DTlist, fill=TRUE)
}


#' Calculate reproduction statistics for consensus clustering between two
#'    cohorts using `reps` random samples.
#'
#' @param conClusters A \code{list} of concensus clustering results, as returned
#'    by the `conClustCohort` function or as a list item from `conClustAllCohorts`.
#' @param procCohort A \code{matrix} of normalized gene expression data with genes
#'    as rows and samples as columns.
#' @param reps The number of permutations to use in `clusterRepro` function. Defaults
#'    to 100.
#' @param seed The integer seed to set when reproducing clusters.
#'
#' @return A \code{data.table} containing the p.value, Number, Actual.IGP,
#'   Actual.Size values for each cluster pair-wise cluster comparison.
#'
#' @importFrom clusterRepro clusterRepro
#' @import data.table
#' @export
calcClusterRepro <- function(conClusters, procCohort, reps=100,
                             seed=NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Extract appropriate data
  # Find all pair-wise K comparisons between the two consensus clusters
  clust1Ks <- seq_len(conClusters[[1]]$optimalK)
  clust2Ks <- seq_len(conClusters[[2]]$optimalK)

  # Extract data
  clust1Centroid <- na.omit(conClusters[[1]]$centroidClusters)
  clust2Data <- na.omit(procCohort[[2]])
  clust2Classes <- na.omit(conClusters[[2]]$classes)

  # Find common genes
  sharedGenes <- intersect(rownames(clust1Centroid), rownames(clust2Data))

  # Subset to common genes
  clust1Centroid <- clust1Centroid[sharedGenes, ]
  clust2Data <- clust2Data[sharedGenes, ]

  # Reproduce clusters
  repClusters <- clusterRepro(clust1Centroid, clust2Data,
                              Number.of.permutations=reps)
  return(as.data.table(repClusters))
}


#' Permutations test for null distribution of of each consensus clustering results
#'
#' Calculate reproducibility statistics for each consensus clustering result
#'    using random samples from the the second cohort in a comparison.
#'
#' @param MSMthresholds A \code{data.table} where columns 1 through 5 represent
#' @param allConClusters A \code{list} of consensus clustering results for
#'     all non-self pair-wise comparisons of cohorts, as returned by the
#'     `conClustAllCohorts` function.
#' @param allProcCohorts A \code{list} of preprocessed gene expression data
#'     for each cohort. As returned by the `preprocCohorts` function.
#' @param reps A \code{numeric} vector containing the integer number of
#'     random samples to reproduce clustering over.
#' @param nthread A \code{numeric} vector containing the integer number of
#'     threads to parallelize over.
#' @param seed A \code{numeric} vector containing the desired seed to be used
#'     for sampling.
#'
#' @importFrom BiocParallel bpplappy DoparParam
#' @importFrom doParallel registerDoParallel
#' @export
calcAllClusterRepro <-  function(MSMthresholds, allConClusters, allProcCohorts,
                                 reps, nthread, seed=NULL) {
  if (!missing(nthread)) {
    opts <- options()
    options(mc.cores=nthread)
    on.exit(options(opts))
  }

  uniqueThresholds <- unique(MSMthresholds[, .(comparison, cohort1, cohort2)])

  registerDoParallel(nthread)
  p <- DoparParam()

  reproList <-
      bplapply(seq_len(nrow(uniqueThresholds)),
               function(idx, thresholds, cluster, data, n, seed) {
                 comp <- thresholds[idx, ]$comparison
                 print(comp)
                 cohortIdxs <- c(thresholds[idx, ]$cohort1,
                               thresholds[idx, ]$cohort2)
                 cl <- calcClusterRepro(cluster[cohortIdxs],
                                        data[cohortIdxs],
                                        reps=n,
                                        seed=seed)
                 cl$comparison <- rep(comp, nrow(cl))
                 cl
               },
               thresholds=uniqueThresholds,
               cluster=allConClusters,
               data=allProcCohorts,
               n=reps,
               seed=seed,
               BPPARAM=p)

  rbindlist(reproList)
}


#' Compare the consensus clustering results for cohort one of a pair
#'    against the null distribution for cohort2
#'
#' @param MSMthresholds A \code{data.table} containing the thresholds between
#'    each cluster for each comparison, as returned by `calcAllMSMthresholds`.
#' @param allClusterRepro A \code{data.table} containing the cluster reproduction
#'    statistics, as calculated with `calcAllClusterRepro`.
#'
#' @importFrom ClusterRepro ClusterRepro
#' @import data.table
#' @export
compareClusters <- function(MSMthresholds, allClusterRepro, pValue, actualIGP, minThresh) {
  # Subset to best clustering match per cluster in cohort 1
  MSMmaxThresholds <- MSMthresholds[, .('c2Clust'=which.max(threshold),
                                        'threshold'=max(threshold)),
                                    by=.(comparison, c1Clust)]
  # Get indexes of significant rows
  sigClustIdx <- allClusterRepro[, na.omit(.I[p.value < pValue & Actual.IGP > actualIGP])]

  # Subset to significant rows and add old idxs to join on
  sigClusters <- allClusterRepro[sigClustIdx, ][, idx := sigClustIdx]
  sigMaxThresholds <- MSMmaxThresholds[sigClustIdx, ][, idx := sigClustIdx]

  # Join on idx, filter below minThresh and drop idx
  clustStats <- merge(sigMaxThresholds, sigClusters, on='idx', sort=FALSE)[threshold > minThresh, -'idx']
  return(clustStats)
}


#' Retrieve the network graph information based on a set of significant cluster
#'    edges.
#'
#' @param clusterEdges A \code{list} of significant cluster edges
#' @param seed An integer seed to set for for the community search.
#'
#' @importFrom igraph graph_from_edgelist layout_with_fr as.undirected
#'     fastgreedy.community
#' @importFrom RColorBrewer brewer.pal
#' @export
getClusterNetworkData <- function(clusterEdges, seed=NULL) {
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

  return(list("graph"=graph, "coords"=coords,
              "ugraph"=ugraph, "metaClusters"=metaClusters))
}

#' Find all possible clusters in all cohorts based on the optimal K value.
#'
#' @param allConClusters A \code{list} of consensus clustering results, as
#'    returned by the `conClustAllCohorts` function.
#'
#' @return A \code{character} vector containing the cluster number appended to
#'    the cohort name.
#'
#' @export
findAllPossibleClusters <- function(allConClusters) {
  optKs <- lapply(allConClusters, `[[`, "optimalK")
  possibleClusters <- lapply(optKs, seq_len)
  clusterList <- lapply(names(allConClusters),
                        function(name, possibleClusters)
                          paste(name, possibleClusters[[name]], sep="-"),
                        possibleClusters=possibleClusters)
  return(unlist(clusterList))
}

#' Return which cohorts were not successfully clustered in `conClustAllCohorts`
#'
#' @param allPossibleClusters A \code{character} vector of cohort clusters,
#'    as returned by the `findAllPossibleClusters` function.
#' @param metaClusters A \code{communities} object from the `igraph` package,
#'    as in the `metaClusters` item in the list returned by the
#'    `getClusterNetworkData` function.
#'
#' @return A \code{character} vector of the cohort clusters which were not
#'    successfully metaclustered.
#'
#' @export
findCohortsNotClustered <- function(allPossibleClusters, metaClusters) {
  setdiff(allPossibleClusters, metaClusters$name)
}

#' Extract the per cohort cluster and meta cluster indexes from metaclustering
#'    results.
#'
#' @param metaClusters A \code{communities} object from the `igraph` package,
#'    as in the `metaClusters` item in the list returned by the
#'    `getClusterNetworkData` function.
#'
#' @return A \code{data.table} with the columns cohort, cohortCluster and
#'    metaCluster columns, matching the clusters from each individual cohort
#'    clustering with the calculated metaclusters.
#'
#' @import data.table
#' @export
getCohortwiseClasses <- function(metaClusters, notClustered) {
  cohortClusters <- strsplit(c(metaClusters$names, notClustered), '-')
  cohortClasses <- data.table(do.call(rbind, cohortClusters))
  cohortClasses$metaCluster <- c(metaClusters$membership, rep("NA", length(notClustered)))
  colnames(cohortClasses) <- c("cohort", "cohortCluster", "metaCluster")
  cohortClasses[, `:=`(cohortCluster=as.numeric(cohortCluster),
                       metaCluster=c(as.numeric(metaCluster)))]
  return(cohortClasses[order(cohort), ])
}

#' Replace the metacluster indexes with name of the metaclusters
#'
#' @param allConClusters A \code{list} of consensus clustering results, as
#'   returned byt he `conClustAllCohorts` function.
#' @param clusterwiseSamples A \code{data.table} with the samples per cluster
#'   in a list column, as returne dby the `getClusterwiseSamples` function.
#'
#' @import data.table
#' @export
annotateSampleMetaClasses <- function(allConClusters, clusterwiseSamples) {
  setorderv(clusterwiseSamples, cols=c('cohort', 'cohortCluster'))
  cohorts <- split(clusterwiseSamples, by='cohort')
  sampleMetaClusters <-
      lapply(cohorts,
             function(cohort, conClusters) {
               cohortClusters <- cohort$cohortCluster
               cohortClusters <- na.omit(cohortClusters)
               samples <- numeric()
               for (cl in cohortClusters) {
                 samplesLength <- length(unlist(cohort[cohortClusters==cl, ]$samples))
                 namedMetaclusters <- rep(cohort[cohortCluster==cl, ]$metaCluster, samplesLength)
                 names(namedMetaclusters) <- unlist(cohort[cohortClusters == cl, ]$samples)
                 samples <- append(samples,
                                   namedMetaclusters)
               }
               return(samples)
             }, conClusters=allConClusters)[names(allConClusters)]
  annotatedClusters <- mapply(.annotateConClusters,
                              cohort=sampleMetaClusters,
                              cluster=allConClusters,
                              SIMPLIFY=FALSE)
}

#' Helper function for `annotatedSampleMetaClasses` to convert metacluster
#'     indexes to names.
#'
#' @param cohort A \code{data.table} of non-NA per cohort metacluster identities,
#'    as in the `data.table` returned from `getClusterwiseSamples` if subset
#'    to a single cohort.
#' @param cluster A \code{list} of consensus clustering results, as retured
#'    by the `conClustCohort` function or an item in list returned by
#'    `conClustAllCohorts`.
#'
#' @keywords internal
.annotateConClusters <- function(cohort, cluster) {
  cluster$metaClasses <- cohort[names(cluster$classes)]
  names(cluster$metaClasses) <- names(cluster$classes)
  return(cluster)
}

#' Append a list column with the vector of sample names in each cohort cluster
#'   to the `data.table` returned by the `getCohortwiseClasses` function.
#'
#' @param cohortwiseClasses A \code{data.table} containing the meta cluster
#'    index for each cohort cluster, as returned by the `getCohortwiseClasses`
#'    function in this package.
#' @param allConClusters A \code{list} of per cohort conensus clustering results,
#'    as returned by the `conClustAllCohorts` function in this package.
#'
#' @return A \code{data.table} equivalent to `cohortwiseClasses`, except with
#'   a \code{list} column called samples, containing the samples for each
#'   cohort per cluster per metacluster.
#'
#' @import data.table
#' @export
getClusterwiseSamples <- function(cohortwiseClasses, allConClusters) {
  metaclusterDTs <- split(cohortwiseClasses, by="metaCluster")
  sampleNames <- lapply(metaclusterDTs,
         function(DT, conClust) {
           cohort <- DT$cohort
           cluster <- DT$cohortCluster
           m <- mapply(.extractSamples,
                       cohort=cohort,
                       cluster=cluster,
                       MoreArgs=list(allConCluster=conClust),
                       SIMPLIFY=FALSE)
           names(m) <- paste0(DT$cohort, DT$cohortCluster)
           m
         },
         conClust=allConClusters)
  for (cluster in unique(cohortwiseClasses$metaCluster)) {
    cohortwiseClasses[metaCluster==cluster, samples := sampleNames[[as.character(cluster)]]]
  }
  return(cohortwiseClasses)
}

#' Get the sample names in a cohort cluster
#'
#' @param allConClusters A \code{list} of per cohort conensus clustering results,
#'    as returned by the `conClustAllCohorts` function in this package.
#' @param cohort A \code{data.table} of non-NA per cohort metacluster identities,
#'    as in the `data.table` returned from `getCohortwiseClasses` if subset
#'    to a single cohort.
#' @param cluster A \code{list} of consensus clustering results, as retured
#'    by the `conClustCohort` function or an item in list returned by
#'    `conClustAllCohorts`.
#'
#' @return a \code{character} vector of sample names for that cohort cluster
#'    combination.
#'
#' @keywords internal
.extractSamples <- function(allConClusters, cohort, cluster) {
  classes <- allConClusters[[cohort]]$classes
  names(classes)[classes == cluster]
}