## ----message=FALSE----------------------------------------------------------------------------------------------------------------------------------------
library(PDATK)
library(Biobase)
library(data.table)
library(ggplot2)
library(randomForest)
library(survival)
library(survminer)
library(gridExtra)
library(SummarizedExperiment)

## ----load_cohort_data-------------------------------------------------------------------------------------------------------------------------------------
metaGxPancreas <- readRDS('../data/metaGxPancreas.rds')

## ----get_common_genes---------------------------------------------------------------------------------------------------------------------------
commonGenes <- getCommonGenes(metaGxPancreas)
length(commonGenes)

## ----extract_gene_expression_profile---------------------------------------------------------------------------------------------------------------------------
cohortsDataL <- subsetCohortExprs(metaGxPancreas, commonGenes)
dim(cohortsDataL$TCGA)
mode(cohortsDataL$TCGA)

## ----ranking_genes_based_on_mad---------------------------------------------------------------------------------------------------------------------------
cohortMADrankings <- rankAllCohortGenesByMAD(cohortsDataL, nthread=15)
head(cohortMADrankings$TCGA[order(cohortMADrankings$TCGA$rank),])

## ----rank_cohort_genes_by_MAD-----------------------------------------------------------------------------------------------------------------------------
genewiseWeightedMADdf <- calcGenewiseWeightedMADdf(cohortMADrankings)
head(genewiseWeightedMADdf)

## ----get_top_ranked_genes---------------------------------------------------------------------------------------------------------------------------------
top20genesPerCohort <- getTopGenes(cohortMADrankings, n=20)
head(top20genesPerCohort$TCGA)

## ----get_top_genes_and_meta_genes-------------------------------------------------------------------------------------------------------------------------
metaGenes <- getMetaGenes(cohortMADrankings, genewiseWeightedMADdf, n=20, pct=20)
length(metaGenes)
resultsDir <- file.path("..", "results")
saveRDS(metaGenes, file.path(resultsDir, "metaGenes.rds"))


## ----preprocess_cohort_expression_matrixes----------------------------------------------------------------------------------------------------------------
preprocCohortL <- preprocCohorts(cohortsDataL, metaGenes)
dim(preprocCohortL$TCGA)


## ----compute_consensus_cluster----------------------------------------------------------------------------------------------------------------------------
consensusClusters <- conClustAllCohorts(preprocCohortL, maxK=5, reps=1000,
                                        distance="pearson", method="hc",
                                        nthread=4, seed=1987)

saveRDS(consensusClusters, file.path("..", "results", "consensusClusters.rds"))


## ----consensus_cluster_normals_and_GTEx_normals-----------------------------------------------------------------------------------------------------------
normalGTEx <- readRDS('../data/normalGTEx.rds')
adjNormal <- readRDS('../data/adjacentNormal.rds')
normalCohortL <- list("normals"=adjNormal,
                      "normalsGTEx"=normalGTEx)
preprocNormalL <- preprocCohorts(normalCohortL, metaGenes)
conClustNormals <- conClustAllCohorts(preprocNormalL, maxK=5, reps=1000,
                                      distance="pearson", method="hc",
                                      nthread=4, seed=1987)


## ----append_normals_to_cohorts----------------------------------------------------------------------------------------------------------------------------
allConClusters <- c(consensusClusters, conClustNormals)
allProcCohorts <- c(preprocCohortL, preprocNormalL)


## ----get_all_pairs----------------------------------------------------------------------------------------------------------------------------------------
cohortPairs <- findAllCohortPairs(names(allConClusters))


## ----calculate_thresholds---------------------------------------------------------------------------------------------------------------------------------
MSMthresholds <- calcAllMSMthresholds(cohortPairs, allConClusters, allProcCohorts, nthread=16)


## ----compare_all_clusters---------------------------------------------------------------------------------------------------------------------------------
if (file.exists(file.path('..', 'results', 'allClusterRepro.rds'))) {
    allClusterRepro <- readRDS(file.path('..', 'results', 'allClusterRepro.rds'))
} else {
    a <- Sys.time()
    allClusterRepro <- calcAllClusterRepro(MSMthresholds, allConClusters,
                                           allProcCohorts, nthread=14, reps=500,
                                           seed=1987)
    b <- Sys.time()
    b - a
}


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(allClusterRepro, file.path('..', 'results', 'allClusterRepro.rds'))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
clusterEdges <- compareClusters(MSMthresholds, allClusterRepro, pValue=0.05, actualIGP=0.5, minThresh=0)
saveRDS(clusterEdges, file.path('..', 'results', 'clusterEdges.rds'))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
graphData <- getClusterNetworkData(clusterEdges, seed=1987)
graphData$metaClusters


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
possibleClusters <- findAllPossibleClusters(allConClusters)
notClustered <- findCohortsNotClustered(possibleClusters, graphData$metaClusters)
notClustered


## ----cohortwise_classes-----------------------------------------------------------------------------------------------------------------------------------
cohortwiseClasses <- getCohortwiseClasses(graphData$metaClusters, notClustered)
cohortwiseClasses


## ----cohortCluster_samples--------------------------------------------------------------------------------------------------------------------------------
clusterwiseSamples <- getClusterwiseSamples(cohortwiseClasses, allConClusters)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
annotatedClusters <- annotateSampleMetaClasses(allConClusters, clusterwiseSamples)
saveRDS(annotatedClusters, file.path("..", "results", "annotatedClusters.rds"))

## ----cluster_network, fig.height=10, fig.width=10---------------------------------------------------------------------------------------------------------
plotClusterNetwork(clusterEdges, seed=1987,
                   clusterLabels=c("Cluster 1", "Cluster 2", "Cluster 3"),
                   saveDir=file.path('..', 'results'),
                   fileName='clusterNetworkGraph.pdf')

## ----load_expression_data---------------------------------------------------------------------------------------------------------------------------------
load(file=file.path('..', 'data', 'PDACexpressionData.rda'))
load(file=file.path('..', 'data', 'yangSurvival.rda'))

## ----extract_survival_data--------------------------------------------------------------------------------------------------------------------------------
survivalData <- extractSurvivalData(metaGxPancreas)

## ----extract_meta_clusters--------------------------------------------------------------------------------------------------------------------------------
metaClusters <- extractMetaClusters(annotatedClusters)
head(metaClusters)

## ----merge_to_samples_with_survival_data------------------------------------------------------------------------------------------------------------------
metaClusterSurvival <- mergeSurvOn(survivalData, metaClusters, "sampleNames")
head(metaClusterSurvival)

## ----annotate_metaClusters--------------------------------------------------------------------------------------------------------------------------------
metaClusterSurvAnnot <- annotateMetaClusters(metaClusterSurvival, c("Cluster 1", "Cluster 2", "Cluster 3"))
head(metaClusterSurvAnnot)

## ----fit_cox_model----------------------------------------------------------------------------------------------------------------------------------------
proportionalHarzardsModel <- fitProportionalHazardsModel(metaClusterSurvAnnot)
summary(proportionalHarzardsModel)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
metaClusterSurvAnnot <- as.data.table(metaClusterSurvAnnot)[metaClusters %in% c("Cluster 1", "Cluster 2", "Cluster 3")]
#metaClusterSurvAnnot <- as.data.table(metaClusterSurvAnnot)[metaClusters %in% c("Cluster 1", "Cluster 2")]
survivalCurves <- fitSurvivalCurves(metaClusterSurvAnnot)

## ----fig.height=6, fig.width=8----------------------------------------------------------------------------------------------------------------------------
plotSurvivalCurves(survivalCurves, saveDir=file.path("..", "results"), fileName="survivalMetaclusters.pdf")

## ---- fig.height=12, fig.width=12-------------------------------------------------------------------------------------------------------------------------
plotCohortwiseSurvCurves(metaClusterSurvAnnot,  saveDir=file.path("..", "results"), fileName="cohortSurvivalMetaclusters.pdf")

## ----get_metaclass_by_cohort_and_sample-------------------------------------------------------------------------------------------------------------------
sampleMetaClassDT <- na.omit(extractSampleMetaClasses(annotatedClusters)[metaClasses %in% c(1, 2, 3)]) # Remove this to include other cluster
sampleMetaClassDT <- na.omit(extractSampleMetaClasses(annotatedClusters)[metaClasses %in% c(1, 2)]) # Remove this


annotSampMetaClassDT <- annotateSampleMetaClassDT(sampleMetaClassDT)
head(annotSampMetaClassDT)


## ----read_in_gene_signature_data--------------------------------------------------------------------------------------------------------------------------
geneSigL <- readRDS(file=file.path("..", "data", "geneSigList.rds"))


computeSigScoreDT_2 <- function(cohortsDataL, sampleMetaClassDT, `signatureGenes`)
{
    cohortsSigGenes <- lapply(cohortsDataL, `[`, i = signatureGenes,
                              j = TRUE)
    cohortsSigGenes <- lapply(cohortsSigGenes, na.omit)
    normalizedCohorts <- normalizeCohortsList(cohortsSigGenes)
    geneScoreList <- lapply(normalizedCohorts, colMeans)
    sampleNames <- unlist(lapply(geneScoreList, names))
    geneScoreDT <- data.table(samples = sampleNames, sigScores = unlist(geneScoreList))
    sampleClassDT <- annotateSampleMetaClassDT(sampleMetaClassDT,
                                               c("Cluster 1", "Cluster 2", "Cluster 3"))
    sigScoreDT <- merge(geneScoreDT[!duplicated(samples)], sampleClassDT[!duplicated(samples)],
                        on = "samples")
}


calcAllGeneSetSigScores_2 <- function(geneSigL, cohortsDataL, sampleMetaClassDT)
{
    mapply(computeSigScoreDT_2, geneSigL, MoreArgs = list(cohortsDataL = cohortsDataL,
                                                          sampleMetaClassDT = sampleMetaClassDT), SIMPLIFY = FALSE)
}


sigScoreL <- calcAllGeneSetSigScores_2(geneSigL, cohortsDataL, annotSampMetaClassDT)
