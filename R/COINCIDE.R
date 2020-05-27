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
calcGeneWiseWeightedMADs <- function(cohortMADrankings) {
  
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
#'     the `rowMADdfs` function.
#'
#' @return A \code{data.frame} with columns genes, weightedMADs and rankings
#'     ordered by decreasing weightedMAD values.
#'
#' @importFrom dplyr dense_rank
#' @export
calcGeneWiseWeightedMADdf <- function(cohortMADrankings) {
  weightedMADs <- calcGeneWiseWeightedMADs(cohortMADrankings)
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

##################### CLUSTER DATA ##########################################################################

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
#' @param suppressPlots A \code{logical} vector indicating whether to suppress
#'    the plots fomr ConsensusClusterPlus. Defaults to TRUE.
#' @param plotTo A \code{character} vector passed to the plot argument of 
#'    `ConsensusClusterPlus`. Defaults to NULL.
#' @param ... Fallthrough arguments to `ConsensusClusterPlus`. If specficied
#'     the function defaults are overridden.
#' 
#' @return A \code{list} of consensus clustering results with items optimumK,
#'     clusterTable, centroidClusters and classes.
#' 
#' @section Warning: This function uses random numbers, remember to set.seed()
#'     before running it to ensure reproducible results!
#' 
#' @importFROM IDPmisc NaRV.omit
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @export
consensusClusterCohort <- function(cohort, maxK, distance, method, 
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
                                     reps=1000,
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

#' 
#' @param i The cluster index
#' @param cohort The cohort data with genes as columns, rows as samples
#' @param clusterClasses The 
#'
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
preprocessCohorts <- function(cohorts, metaGenes, center=TRUE, scale=FALSE) {
  cohortSubsets <- .subsetCohorts(cohorts, genes=metaGenes)
  structure(lapply(cohortSubsets, .scaleGenewise,
                   center=center, scale=scale),
            .Names=names(cohorts))
}

#' Convert to matrix, transform, scale row-wise, transform, return
#'
#' @param cohortSubet A \code{data.frame} containing the cohort
#'    data subset to the genes of interest. Rows are features,
#'    columns are samples.
#'
#' @export
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
#' @export
#' @keywords internal
.subsetCohorts <- function(cohorts, genes=TRUE, columns=TRUE) {
  lapply(cohorts, `[`, i=genes, j=columns)
}

#'
#'
#'
#'
#'
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

#'
#'
#'
#'
#'
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
#'    cohorts using `numReps` random samples.
#'
#' @param conClusters
#' @param procCohorts
#' @param numReps
#' @param seed
#'
#' @export
calcClusterRepro <- function(conClusters, procCohorts, numReps=100, 
                             seed=NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Extract appropriate data
  # Find all pair-wise K comparisons between the two consensus clusters
  clust1Ks <- seq_len(conClusters[[1]]$optimalK)
  clust2Ks <- seq_len(conClusters[[2]]$optimalK)
  
  # Extract data
  clust1Centroid <- na.omit(conClusters[[1]]$centroidClusters)
  clust2Data <- na.omit(procCohorts[[2]])
  clust2Classes <- na.omit(conClusters[[2]]$classes)
  
  # Find common genes
  sharedGenes <- intersect(rownames(clust1Centroid), rownames(clust2Data))
  
  # Subset to common genes
  clust1Centroid <- clust1Centroid[sharedGenes, ]
  clust2Data <- clust2Data[sharedGenes, ]
  
  # Reproduce clusters
  repClusters <- clusterRepro(clust1Centroid, clust2Data, 
                              Number.of.permutations=numReps)
  return(as.data.table(repClusters))
}


#' Permutations test for null distribution of of each consensus clustering results 
#' 
#' Calculate reproducibility statistics for each consensus clustering result 
#'    using random samples from the the second cohort in a comparison.
#' 
#' @param MSMthresholds A \code{data.table} where columns 1 through 5 represent
#' @param allConClusters A \code{list} of consensus clustering results for
#'     all non-self pair-wise comparisons of cohorts.
#' @param allProcCohorts A \code{list} of preprocessed expression cohorts
#' @param numReps A \code{numeric} vector containing the integer number of 
#'     random samples to reproduce clustering over.     
#' @param nthread A \code{numeric} vector containing the integer number of 
#'     threads to parallelize over.
#' @param seed A \code{numeric} vector containing the desired seed to be
#'     for sampling.
#'
#' @export
calcAllClusterRepro <-  function(MSMthresholds, allConClusters, allProcCohorts, 
                                 numReps, nthread, seed=NULL) {
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
                                        numReps=n,
                                        seed=seed)
                 cl$comparison <- rep(comp, nrow(cl))
                 cl
               },
               thresholds=uniqueThresholds,
               cluster=allConClusters,
               data=allProcCohorts,
               n=numReps,
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
#' @export
compareClusters <- function(MSMthresholds, allClusterRepro, pValue, actualIGP, minThresh) {
  # Subset to best clustering match per cluster in cohort 1
  MSMmaxThresholds <- MSMthresholds[, .('c2Clust'=which.max(threshold), 'threshold'=max(threshold)),
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
#'
#' @param clusterEdges
#' @param seed
#'
#' @
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

#'
#'
#'
#'
findAllPossibleClusters <- function(allConClusters) {
  optKs <- lapply(allConClusters, `[[`, "optimalK")
  possibleClusters <- lapply(optKs, seq_len)
  clusterList <- lapply(names(allConClusters), 
                        function(name, possibleClusters)
                          paste(name, possibleClusters[[name]], sep="-"),
                        possibleClusters=possibleClusters)
  return(unlist(clusterList))
}

#'
#'
#'
#'
findCohortsNotClustered <- function(allPossibleClusters, metaClusters) {
  setdiff(allPossibleClusters, metaClusters$name)
}

#'
#'
#'
#'
getCohortwiseClasses <- function(metaClusters, notClustered) {
  cohortClusters <- strsplit(c(metaClusters$names, notClustered), '-')
  cohortClasses <- data.table(do.call(rbind, cohortClusters))
  cohortClasses$metaCluster <- c(metaClusters$membership, rep("NA", length(notClustered)))
  colnames(cohortClasses) <- c("cohort", "cohortCluster", "metaCluster")
  cohortClasses[, `:=`(cohortCluster=as.numeric(cohortCluster),
                       metaCluster=c(as.numeric(metaCluster)))]
  return(cohortClasses[order(cohort), ])
}

#'
#'
#' @param
#' @param 
#'
#'
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

.annotateConClusters <- function(cohort, cluster) {
  cluster$metaClasses <- cohort[names(cluster$classes)]
  names(cluster$metaClasses) <- names(cluster$classes)
  return(cluster)
}

#' 
#'
#' @param cohortwiseClasses
#' @param allConClusters
#'
#' @return A \code{data.table} equivalent to `cohortwiseClasses`, except with
#'   a \code{list} column called samples, containing the samples for each
#'   cohort per cluster per metacluster.
#'
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
#' @param allConClusters
#' @param cohort
#' @param cluster
#'
#' @return a \code{character} vector of sample names for that cohort cluster
#'    combination.
#' 
.extractSamples <- function(allConClusters, cohort, cluster) {
  classes <- allConClusters[[cohort]]$classes
  names(classes)[classes == cluster]
}