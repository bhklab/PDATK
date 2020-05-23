########################################################################################################
### Required libraries
#########################################################################################################

require(IDPmisc)
library(dplyr)
library(matrixStats)
library(factoextra)
library(NbClust)
library(cluster)
library(ConsensusClusterPlus)
library(clusterRepro)
library(igraph)
library(gtools)
library(plyr)
library(survival)
library(KMsurv)
library(survminer)
library(limma)
library(piano)
require(ggplot2)
library(ggpubr)
library(doParallel)
library(foreach)
library(effsize)
library(vcdExtra)
library(survcomp)
library(RColorBrewer)
library(lattice)

#########################################################################################################
### Loading Normal, GTEx and other cohorts
#########################################################################################################

# load('../data/common_genes_cohorts_new.RData')       
# load('../data/normal_GTEx.RData')
# load('../data/adjacent_normal.RData')


#########################################################################################################
### Sourcing helping functions
#########################################################################################################
source('../R/cluster_me_cont.R')
source('../R/comare_interclusters_cont_pearson.R')
source('../R/mgsub_function.R')
source('../R/msm_threshold.r')


########################################################################################################
### Ranking the genes based on MAD
#########################################################################################################


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
  set.seed("1987", sample.kind="Rounding")
  
  # Get the name of the comparison and a vector of the cohort list indexes
  comparison <- names(cohortPair)
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
  meanCors <- vector("numeric", max(clust2Ks))
  for (i in clust1Ks) {
    for (j in clust2Ks) {
      classIdxs <- which(clust2Classes == j)
      centroid <- clust1Centroid[, i]
      meanCors[j] <- mean(unlist(lapply(classIdxs,
                            function(idx, centroid, data) {
                              as.numeric(
                                cor.test(
                                  centroid,
                                  data[, idx],
                                  method="pearson")$estimate
                                )}, 
                            centroid=centroid,
                            data=clust2Data)))
    }
  }
  return(meanCors)
}

#' 
#' 
#' 
#' 
#' 
#' 
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
  DTlist
  # rbindlist(DTlist, fill=TRUE)[, `:=`(comparison=rownames(cohortPairs), 
  #                                     x=cohortPairs$x, y=cohortPairs$y)]
}



calcClusterRepro <- function(MSMthreshold, allConClusters, allProcCohorts, 
                             seed, sampleKind) {
  
  if (!missing(seed)) set.seed(seed, sample.kind=sampleKind)
  
  # Get cluster indexes
  cohort1 <- MSMthreshold$x
  cohort2 <- MSMthreshold$y
  
  # Get cluster names
  cohNames <- unlist(strsplit(MSMthreshold$comparison, '-'))
  
  # Extract appropriate data
  # Find all pair-wise K comparisons between the two consensus clusters
  clust1Ks <- seq_len(allConClusters[[cohort1]]$optimalK)
  clust2Ks <- seq_len(allConClusters[[cohort2]]$optimalK)
  
  # Extract data
  clust1Centroid <- na.omit(allConClusters[[cohort1]]$centroidClusters)
  clust2Data <- na.omit(allProcCohorts[[cohort2]])
  clust2Classes <- na.omit(allConClusters[[cohort2]]$classes)
  
  # Find common genes
  sharedGenes <- intersect(rownames(clust1Centroid), rownames(clust2Data))
  
  # Subset to common genes
  clust1Centroid <- clust1Centroid[sharedGenes, ]
  clust2Data <- clust2Data[sharedGenes, ]
  
  # Reproduce clusters
  repClusters <- clusterRepro(clust1Centroid, clust2Data, 
                              Number.of.permutations=500)
  return(repClusters)
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
#'     
#' @param nthread 
#' @param seed
#' @param sampleKind
#' @export
calcAllClusterRepro <-  function(MSMthresholds, allConClusters, allProcCohorts, 
                                 nthread, seed, sampleKind=NULL) {
  if (!missing(nthread)) {
    opts <- options()
    options(mc.cores=nthread)
    on.exit(options(opts))
  }
  
  if (missing(seed)) {
    reproList <- bplapply(seq_len(nrow(MSMthresholds)),
                          function(idx, thresholds, cluster, data) {
                            print(thresholds[idx, ]$comparison)
                            calcClusterRepro(thresholds[idx, ], cluster, data)
                          },
                          thresholds=MSMthresholds,
                          cluster=allConClusters,
                          data=allProcCohorts)
  } else {
    reproList <- 
      bplapply(seq_len(nrow(MSMthresholds)),
               function(idx, thresholds, cluster, data, seed, kind) {
                        print(thresholds[idx, ]$comparison)
                        calcClusterRepro(thresholds[idx, ], cluster, data, 
                                         seed=seed, sampleKind=sampleKind)
                        },
               thresholds=MSMthresholds,
               cluster=allConClusters,
               data=allProcCohorts,
               seed=seed,
               kind=sampleKind)
  }
  reproList
}


#' Compare the consensus clustering results for cohort one of a pair 
#'    against the null distribution for cohort2
#' 
#' 
#' 
#' @param MSMthresholds
#' 
#' @importFrom ClusterRepro ClusterRepro
#' @export
compareClusters <- function(MSMthreshold, conClusterResult, 
                            seed, sampleKind=NULL) {
  
  # Get cluster indexes
  cohort1 <- MSMthreshold$x
  cohort2 <- MSMthreshold$y
  
  # Get cluster names
  cohNames <- unlist(strsplit(MSMthreshold$comparison, '-'))
  
  # Extract appropriate data
  # Find all pair-wise K comparisons between the two consensus clusters
  clust1Ks <- seq_len(allConClusters[[cohort1]]$optimalK)
  clust2Ks <- seq_len(allConClusters[[cohort2]]$optimalK)
  
  # Extract data
  clust1Centroid <- na.omit(allConClusters[[cohort1]]$centroidClusters)
  clust2Data <- na.omit(allProcCohorts[[cohort2]])
  clust2Classes <- na.omit(allConClusters[[cohort2]]$classes)
  
  # Find common genes
  sharedGenes <- intersect(rownames(clust1Centroid), rownames(clust2Data))
  
  # Subset to common genes
  clust1Centroid <- clust1Centroid[sharedGenes, ]
  clust2Data <- clust2Data[sharedGenes, ]
  
  # Extract non-null threshold values
  thresholds <- as.numeric(MSMthreshold[, .SD, .SDcols=clust2Ks])
  thresholds <- replace(thresholds, is.na(thresholds), 0)
  # Fill missing thresholds if cluster 1 has more ks than cluster 2
  # This is to make the for loop work
  if (length(thresholds) < length(clust1Ks)) {
    .fillVector(thresholds, length(clust1Ks), 0)
  }
  
  results <- vector("list", length(clust1Ks))
  cummThresh <- numeric()
  for (i in clust1Ks) {
    cummThresh <- append(cummThresh, thresholds[i])
    if (!is.na(repClusters$p.value[i])) {
      if (repClusters$Actual.IGP[i] > 0.5 &&
          repClusters$p.value[i] < 0.05 &&
          max(cummThresh) > 0) {
        results[[i]] <-
          c(paste(cohNames[1], i, sep="-"),
            paste(cohNames[2], which.max(cummThresh), sep="-"),
            max(cummThresh),
            repClusters$Actual.IGP[i],
            repClusters$p.value[i]
          )
      }
    }
  }
  notNull <- results[!vapply(results, is.null, logical(1))]
  if (length(notNull) == 0) {
    # Create an empty data.table with 8 named columns and return
    return(
      data.table("coh1Cluster"=numeric(), 
                 "coh2Cluster"=numeric(), 
                 "maxThreshold"=numeric(), 
                 "actualIGP"=numeric(), 
                 "pValue"=numeric(),
                 "comparison"=numeric(),
                 "cohort1"=numeric(),
                 "cohort2"=numeric())
    )
  }
  # Assemble the data.table for this cluster comparison
  DTlist <- lapply(notNull, function(vec) transpose(data.table(vec)))
  resultDT <- rbindlist(DTlist)
  colnames(resultDT) <- c("coh1Cluster", "coh2Cluster", "maxThreshold", "actualIGP", "pValue")
  resultDT[, `:=`(comparison=rep(MSMthreshold$comparison, length(DTlist)),
                  cohort1=rep(cohort1, length(DTlist)),
                  cohort2=rep(cohort2, length(notNull)))]
  return(resultDT)
}







pdf("../results/densityplot.pdf")

densityplot(unlist(threshold_msm ))

cores=detectCores()
cl <- makeCluster(cores[1]-1, outfile="") #not to overload your computer
registerDoParallel(cl)

###################################3

edges=list()
edges= foreach(i = 1:dim(allPairs)[1], 
               .packages = c( "doParallel",
                              "ConsensusClusterPlus", 
                              "clusterRepro",
                              "foreach")) %dopar% {                                                        
  compare_interclusters(datasets[allPairs[i,1]], datasets[allPairs[i,2]], 
                        clusters[[allPairs[i,1]]],clusters[[allPairs[i,2]]],
                        data.frame(data[allPairs[i,2]]),i)
  
}
edges_define = matrix(unlist(edges), ncol = 5, byrow = TRUE)
save(edges_define, file='../results/edges.RData')

######### FIND EDGES ###############################
edges_defined = edges_define[,1:2]

######### FIND METACLUSTERS ########################
g = graph_from_edgelist(edges_defined)
coords = layout_with_fr(g)


#### Fast greedy

g1=as.undirected(g)
c5=fastgreedy.community(g1,weights = as.numeric(edges_define[,3]))
table(membership(c5))
plot(c5, g1, layout=coords,vertex.size=6)

l <- layout_with_fr(g1)
#plot(g1, vertex.color=membership(c5),  vertex.size=6,edge.arrow.size=0.5,vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.4, vertex.color="gray",edge.arrow.size=.4,layout=l)


colrs <-  c( "blue", brewer.pal(n = 8, name = "Dark2")[c(4,5)])
shapes <- c("circle", "square", "sphere")
plot(g1, vertex.color=colrs[membership(c5)],  vertex.shape= "sphere", vertex.size=6,edge.arrow.size=0.5,vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.1, vertex.color=colrs[membership(c5)],
     edge.arrow.size=.4, layout=layout.davidson.harel(g1))


zz=sapply(1:length(membership(c5)), function(x) strsplit(names(membership(c5))[x], "-")[[1]][1])
x1= match(zz, unique(zz))

cc= rainbow(22)

pdf("../results/Meta_clusters_Network_community.pdf")

plot.igraph(g1, vertex.color=colrs[membership(c5)],  vertex.shape= "sphere", vertex.size=6,edge.arrow.size=0.5,vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.1, vertex.color=colrs[membership(c5)],
     edge.arrow.size=.4, layout=layout_with_dh(g1), layout=coords)
legend('topleft',legend=c("Basal","Exocrine","Classical"), pt.cex=1.8,pch=21, pt.bg=colrs, bty='n', col=colrs)

################## Samples not clustered
datasets = c("ICGC_seq","PCSI","TCGA","Kirby","OUH","Winter","Collisson","Zhang","Chen","UNC","ICGC_arr","Balagurunathan","Pei","Grutzmann","Badea", "haider","lunardi","yang","hamidi","janky","bauer", "Normals","GTEX_Normals")
clusters = list(z_ICGC, z_PCSI, z_TCGA, z_Kirby, z_OUH, z_winter,z_Collisson,z_zhang,z_chen,z_unc,z_icgc_arr,z_balagurunathan, z_pei, z_grutzmann, z_badea, z_haider, z_lunardi, z_yang, z_hamidi, z_janky, z_bauer, z_normal,z_gtex)

total_clusters=list()
for(i in 1:length(datasets)){
  
  total_clusters[[i]]=paste(datasets[i], 1:clusters[[i]]$optimumK, sep="-" )
  
}

not_found=setdiff(unlist(total_clusters), c5$names)

not_found

########### Samples in classes ###################
datasets = c("ICGC_seq","PCSI","TCGA","Kirby","OUH","Winter","Collisson","Zhang","Chen","UNC","ICGC_arr","Balagurunathan","Pei","Grutzmann","Badea", "haider","lunardi","yang","hamidi","janky","bauer", "Normals","GTEX_Normals")
clusters = list(z_ICGC, z_PCSI, z_TCGA, z_Kirby, z_OUH, z_winter,z_Collisson,z_zhang,z_chen,z_unc,z_icgc_arr,z_balagurunathan, z_pei, z_grutzmann, z_badea, z_haider, z_lunardi, z_yang, z_hamidi, z_janky, z_bauer, z_normal,z_gtex)

samples=list()
for(i in 1: length(unique(membership(c5)))){
  
  clust = names(membership(c5))[which(membership(c5)==i)]  
  
  dataset_no = sapply(1:length(clust), function(x) which(strsplit(clust[x],"-")[[1]][1] == datasets))
  cluster_no = as.numeric(sapply(1:length(clust), function(x) strsplit( names(membership(c5))[which(membership(c5)==i)],"-")[[x]][2]))
  
  samples[[i]]=sapply(1:length(dataset_no),function(x) names(clusters[[dataset_no[x]]]$classes)[which(clusters[[dataset_no[x]]]$classes == cluster_no[x])])
  
}

sample_groups<- list()
for(i in 1:length(samples)){
  
  sample_groups[[i]]=unlist(samples[[i]])
  
}

#save(sample_groups, file="/Users/vandanasandhu/Desktop/Subtyping_PDACs/PDAC_meta_Classes.RData")

########### Samples classes cohort wise #################################################
clusters = list(PCSI= z_PCSI, TCGA=z_TCGA, ICGC_seq=z_ICGC, 
                Kirby = z_Kirby, OUH= z_OUH, Winter= z_winter,
                Collisson=z_Collisson,Zhang=z_zhang,Chen=z_chen,
                UNC=z_unc,ICGC_arr=z_icgc_arr,Balagurunathan=z_balagurunathan, 
                Grutzmann=z_grutzmann,Badea= z_badea,Pei=z_pei,
                hamidi=z_hamidi,yang=z_yang, lunardi=z_lunardi,
                janky=z_janky, bauer=z_bauer, haider=z_haider)
#, Normals=z_normal, GTEX_Normals= z_gtex
#)


x=unlist(strsplit(names(membership(c5)),"-"))
mm=matrix(unlist(x), ncol=2, byrow = TRUE)
mm=data.frame(mm)
colnames(mm)= c("dataset","cohort_clusters")

mm$meta_classes=as.numeric( membership(c5))
mm=mm[order(mm$dataset),]

### ADDING NOT FOUND ROWS
dim(mm)

j=1
for(i in (dim(mm)[1] +1 ):(dim(mm)[1]+length(not_found))){
  
  pp = unlist(strsplit(not_found[j],"-"))
  mm[i,]=c(pp[1],pp[2],NA)
  j=j+1
}

for(i in 1:length(clusters)){
  
  z= which(mm[,1] %in% names(clusters)[i])
  clusters[[i]]$meta_classes=mapvalues(clusters[[i]]$classes, mm[,2][z], mm[,3][z]) 
  
}

save(clusters,file= '../results/clusters.RData')

