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

rm(list=ls(all=TRUE))

# -------------------------------------------------------------------------
# 1. COINCIDE Implementation ----------------------------------------------
# -------------------------------------------------------------------------

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
genewiseWeightedMADdf <- calcPerGeneWeigthedMAD(cohortMADrankings)
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


# -------------------------------------------------------------------------
# 2. Network Community Search ---------------------------------------------
# -------------------------------------------------------------------------


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



# -------------------------------------------------------------------------
# 3. Gene Signature Comparison -----------------------------------------------
# -------------------------------------------------------------------------


## ----get_metaclass_by_cohort_and_sample-------------------------------------------------------------------------------------------------------------------
sampleMetaClassDT <- na.omit(extractSampleMetaClasses(annotatedClusters)[metaClasses %in% c(1, 2, 3)]) # Remove this to include other cluster
annotSampMetaClassDT <- annotateSampleMetaClassDT(sampleMetaClassDT, c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"))
annotSampMetaClassDT <- annotSampMetaClassDT[metaClasses %in% c("Cluster 1", "Cluster 2", "Cluster 3")] # Remove this to include other cluster
head(annotSampMetaClassDT)


## ----read_in_gene_signature_data--------------------------------------------------------------------------------------------------------------------------
geneSigL <- readRDS(file=file.path("..", "data", "geneSigList.rds"))
sigScoreL <- calcAllGeneSetSigScores(geneSigL, cohortsDataL, annotSampMetaClassDT)


## ----plot_gene_signatures_between_metaclasses, fig.height=10, fig.width=12--------------------------------------------------------------------------------
comparisons <- list(c("Cluster 1", "Cluster 2"),
                    c("Cluster 2", "Cluster 3"),
                    c("Cluster 3", "Cluster 1"))
plotSigScoreL(sigScoreL, comparisons, saveDir=file.path("..", "data"), fileName="sigScorePlots.pdf")




# -------------------------------------------------------------------------
# 4. Biomarker Comparison -------------------------------------------------
# -------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
biomarkers <- c("EGFR", "SNAI2", "KRT81", "GATA6", "CYP3A5", "HNF1A",
                         "REG1A", "REG1B", "REG3A", "CEL")
names(biomarkers) <- c(rep("Cluster 1", 3), rep("Cluster 2", 3), rep("Cluster 3", 4)) # Maybe this should be Basal,
                                                                                      # Classical and Exocrine since
                                                                                      # the genes are known biomarkers
                                                                                      # for those subtypes?
geneBiomarkerScoreL <- lapply(biomarkers, computeGeneBiomarkerScores,
                              cohortsDataL=cohortsDataL,
                              annotSampMetaClassDT=annotSampMetaClassDT)
names(geneBiomarkerScoreL) <- biomarkers


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
biomarkerScoreBoxplotL <- boxplotBiomarkerScoresL(geneBiomarkerScoreL,
                                                  biomarkers,
                                                  comparisons)


## ----basal_class_boxplots, fig.height=10, fig.width=10----------------------------------------------------------------------------------------------------
## TODO:: Put this inside a function to reduce libraries needed for analysis
scoreClasses <- names(biomarkerScoreBoxplotL)
scoreIdxs <- which(scoreClasses %in% "Cluster 1")
cluster1Plots <- grid.arrange(grobs=biomarkerScoreBoxplotL[scoreIdxs], ncol=2)
ggsave(file.path("..", "results", "basalBiomarkerBoxPlots.pdf"))
grid::grid.draw(cluster1Plots)


## ----classical_class_boxplots, fig.height=10, fig.width=10------------------------------------------------------------------------------------------------
scoreIdxs1 <- which(scoreClasses %in% "Cluster 2")
cluster2Plots <- grid.arrange(grobs=biomarkerScoreBoxplotL[scoreIdxs], ncol=2)
ggsave(file.path("..", "results", "classicalBiomarkerBoxPlots.pdf"))
cluster2Plots


## ----exocrine_class_boxplots, fig.height=10, fig.width=10-------------------------------------------------------------------------------------------------
scoreIdxs <- which(scoreClasses %in% "Cluster 3")
cluster3Plots <- grid.arrange(grobs=biomarkerScoreBoxplotL[scoreIdxs], ncol=2)
ggsave(file.path("..", "results", "exocrineBiomarkerBoxPlots.pdf"))
cluster3Plots



# -------------------------------------------------------------------------
# 5. Single Sample Classifier ---------------------------------------------
# -------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
cohortCommMat <- computeCohortCommMat(cohortsDataL, metaGenes)


## ----train_single_sample_classifiers----------------------------------------------------------------------------------------------------------------------
uniqueMetaClasses <- na.omit(unique(sampleMetaClassDT$metaClasses))
singleSampleClassifiers <- lapply(uniqueMetaClasses,
                                  trainSingleSampleClassifer,
                                  cohortComMat=cohortCommMat,
                                  sampleMetaClassDT=sampleMetaClassDT,
                                  maxK=20)
names(singleSampleClassifiers) <- uniqueMetaClasses

## ----get_predicted_genes----------------------------------------------------------------------------------------------------------------------------------
classTSPs <- lapply(singleSampleClassifiers, `[[`, "TSPs")
topGenes1 <- lapply(classTSPs, `[`, i=TRUE, j=1)
topGenes2 <- lapply(classTSPs, `[`, i=TRUE, j=2)


## ----get_top_features-------------------------------------------------------------------------------------------------------------------------------------
lengthTopGenes <- lapply(topGenes1, length)
topGenesClasses <- mapply(rep, names(topGenes1), times=lengthTopGenes, SIMPLIFY=FALSE)
featuresDT <- data.table("class"=unlist(topGenesClasses),
                         "genes"=unlist(lapply(c(topGenes1, topGenes2), unlist)))
genePairs <- paste(unlist(topGenes1), unlist(topGenes2), sep=">")


## ----subset_matrix_to_top_genes---------------------------------------------------------------------------------------------------------------------------
keepSamples <- intersect(colnames(cohortCommMat), unique(sampleMetaClassDT$samples))
cohortTopGenesMat <- cohortCommMat[featuresDT$genes, keepSamples]
sampleMetaClassDT <- sampleMetaClassDT[!duplicated(samples) & samples %in% keepSamples, ]
annotTopGeneMetaClassDT <- annotateSampleMetaClassDT(sampleMetaClassDT, c("Basal", "Classical", "Exocrine"))
metaClassesFactor <- as.factor(annotTopGeneMetaClassDT$metaClasses)


## ----compare_top_genes------------------------------------------------------------------------------------------------------------------------------------
g1Gt2Matrix <- calcTopGenes1Gt2Matrix(cohortTopGenesMat, topGenes1, topGenes2, genePairs)
g1Gt2Matrix[1:5, 1:5]


## ----make_model_predictions-------------------------------------------------------------------------------------------------------------------------------
if (file.exists(file.path("..", "results", "sampleClassPreds.rds"))) {
    sampleClassPreds <- readRDS(file.path("..", "results", "sampleClassPreds.rds"))
} else {
    a <- Sys.time()
    sampleClassPreds <- predictSingleSampleClasses(g1Gt2Matrix, metaClassesFactor, nthread=12)
    saveRDS(sampleClassPreds, file.path('..', 'results', 'sampleClassPreds.rds'))
    b <- Sys.time()
    b - a
}


## ----predict_binary_rf_model------------------------------------------------------------------------------------------------------------------------------
binaryRFmodel <- randomForest(t(g1Gt2Matrix), metaClassesFactor)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
topGenesDT <- data.table("topGenes1"=unlist(topGenes1), "topGenes2"=unlist(topGenes2), "genePairs"=genePairs, "classes"=unlist(topGenesClasses))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
classPredDT <- data.table(sampleClassPreds)
classPredDT$samples <- names(sampleClassPreds)
singleSampClassifDT <- merge(classPredDT, sampleMetaClassDT, on=samples, sort=FALSE)




# -------------------------------------------------------------------------
# 6. Compass Subtyping Waterfall ------------------------------------------
# -------------------------------------------------------------------------


## ----compass_load_data------------------------------------------------------------------------------------------------------------------------------------
compassSurvDT <- fread(file.path("..", "data", "compassSurvDT.csv"))
compassLogExprMat <- read.csv(file.path("..", "data", "compassLogExprMat.csv"),
                              row.names="X")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
sampleMetaClassPredDT <- predictSampleMetaClass(compassLogExprMat, binaryRFmodel,
                                                topGenesDT, g1Gt2Matrix,
                                                metaClassesFactor)
table(sampleMetaClassPredDT$predClass)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
sampClassPredwSurvivalDT <- merge(sampleMetaClassPredDT, compassSurvDT, by="sample")
setorderv(sampClassPredwSurvivalDT, cols="tumorResponse", order=-1)
head(sampClassPredwSurvivalDT)


## ----tumor_response_all_drug_all_classes, fig.height=8, fig.width=12--------------------------------------------------------------------------------------
waterfallPlotTumorResponse(sampClassPredwSurvivalDT, saveDir=resultsDir,
                           fileName="compassTumorResponseWaterfallPlot.pdf")


## ----tumor_response_basalvclassical_ffx, fig.height=8, fig.width=12---------------------------------------------------------------------------------------
## TODO:: Make a function to compare the survival
survBasalvClassicalDT <-
  sampClassPredwSurvivalDT[predClass %in% c("Basal", "Classical") &
                             drug %in% c("FFx", "FFx x 1")]

classSurvComparison <- survfit(Surv(OS, OSstatus) ~ predClass,
                               data=survBasalvClassicalDT)

pValue <- surv_pvalue(classSurvComparison)$pval.txt

waterfallPlotTumorResponse(survBasalvClassicalDT, noXaxis=TRUE, pVal=pValue,
                           saveDir=resultsDir,
                           fileName="basalClassicalTumorResponseWaterfallPlot.pdf")


## ----prop_hazards_model-----------------------------------------------------------------------------------------------------------------------------------
classSurvProHarzardsModel <- coxph(Surv(OS, OSstatus) ~ predClass,
                                  data=survBasalvClassicalDT)
summary(classSurvProHarzardsModel)


## ----compass_NLR------------------------------------------------------------------------------------------------------------------------------------------
classNLRcompDT <- sampClassPredwSurvivalDT[order(NLR, decreasing=TRUE), ]
classNLRcompDT <-classNLRcompDT[!is.na(NLR), ]
boxplotClassNLR(classNLRcompDT)


## ----fig.height=8, fig.width=12---------------------------------------------------------------------------------------------------------------------------
compassPlotBiomarkers(compassLogExprMat, sampClassPredwSurvivalDT,
                      biomarkers=biomarkers[names(biomarkers) %in% c("Cluster 1", "Cluster 2")])




# -------------------------------------------------------------------------
# 7. PharmacoGx Cellline Subtyping ----------------------------------------
# -------------------------------------------------------------------------


## ----load_and_format_data---------------------------------------------------------------------------------------------------------------------------------
PGxCelllinesL <- readRDS(file.path('..', 'data', "PharmacoGxCelllines.rds"))
PGxCLLAvgRepsL <- lapply(PGxCelllinesL, limma::avereps)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
PGxCCLpreds <- lapply(PGxCLLAvgRepsL, predictSampleMetaClass,
                      classifModel=binaryRFmodel, topGenesDT=topGenesDT,
                      trainData=g1Gt2Matrix, trainLabels=metaClassesFactor)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(PGxCCLpreds, file.path("..", "results", "PGxSubtypePredictions"))




# -------------------------------------------------------------------------
# 8. PharmacoGx Cellline Drug Sensitivity ---------------------------------
# -------------------------------------------------------------------------


## ----load_data_and_process--------------------------------------------------------------------------------------------------------------------------------
CCLsensitivityL <- readRDS(file.path("..", "data", "CCLdrugSensitivityAUCs.rds"))


## ----format_data------------------------------------------------------------------------------------------------------------------------------------------
CCLsensitivityDtL <- formatCCLsensitivityLtoDT(CCLsensitivityL)


## ----merge_drugs_and_classes------------------------------------------------------------------------------------------------------------------------------
mergedCelllineDT <- mergeClassAndDrugs(PGxCCLpreds, CCLsensitivityDtL)


## ----get_wilcox_test_pvals--------------------------------------------------------------------------------------------------------------------------------
wilcoxPvals <- computeWilcoxOnSharedDrugs(mergedCelllineDT,
                                          class1="Basal", class0="Classical")
wilcoxPvals$CTRPv2


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
concordanceInds <- computeConcIndOnSharedDrugs(mergedCelllineDT, class="Basal")
concordanceInds


## ----meta_stats-------------------------------------------------------------------------------------------------------------------------------------------
metaStatsDT <- calcMetaStats(concordanceInds)


## ----AUC_boxplots-----------------------------------------------------------------------------------------------------------------------------------------
# Exclude GDSC due to NA stats
plot <- boxplotAUCperSubtypePerDataset(mergedCelllineDT[dataset != "GDSC"],
                                       concordanceInds[names(concordanceInds) != "GDSC"])


## ----fig.height=12, fig.width=12--------------------------------------------------------------------------------------------------------------------------
## TODO:: Put this into a function
dims <- ceiling(sqrt(length(plot)))
ggarrange(plotlist=plot, ncol=dims, nrow=dims)




# -------------------------------------------------------------------------
# 9. PDXE Drug Subtyping --------------------------------------------------
# -------------------------------------------------------------------------


## ----load_PDX_data----------------------------------------------------------------------------------------------------------------------------------------
# Why is first dataset here? Nothing is done with it?
PDXgeneExprDT <- fread(file.path("..", "data", "PDXgeneExpr.csv"))
PDXgeneExprSet <- readRDS(file.path("..", "data", 'PDXexprSet.rds'))
PDXdrugResponse <- fread(file.path('..', 'data', 'PDXdrugResponse.csv'))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
PDXgeneExprClassDT1 <- preprocPDXdata(PDXgeneExprDT, binaryRFmodel, topGenesDT,
                                      g1Gt2Matrix, metaClassesFactor)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
PDXgeneExprClassDT2 <- preprocPDXdata(exprs(PDXgeneExprSet), binaryRFmodel,
                                      topGenesDT, g1Gt2Matrix,
                                      metaClassesFactor)
table(PDXgeneExprClassDT2$predClass)


## ----format_drug_response---------------------------------------------------------------------------------------------------------------------------------
PDXdrugResp <- PDXdrugResponse[, .(patient.id, drug, AAC)]
colnames(PDXdrugResp) <- c('sample', 'drug', 'AAC')


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
PDXmergedDT <- merge(PDXgeneExprClassDT2, PDXdrugResp, by="sample")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
plots <- boxplotPDXsubtypePerDrug(PDXmergedDT)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
## TODO:: Add plot saving arguments to boxplotDPXsubtypePerDrug
plotGrids <- ggarrangePlotL(plots, 5) # This function is buggy right now, need to fix
for (i in seq_along(plotGrids)) {
  ggsave(plotGrids[[i]], file=file.path(resultsDir, paste0('PDXdrugRespPerDrug', i, ".pdf")))
}

## ----fig.height=12, fig.width=8---------------------------------------------------------------------------------------------------------------------------
plotGrids




# -------------------------------------------------------------------------
# 10. Cellularity Comparison ----------------------------------------------
# -------------------------------------------------------------------------


## ----load_cellularity_data--------------------------------------------------------------------------------------------------------------------------------
cellularityFileL <- list.files(file.path("..", "data"),
                               pattern="cellularity.*tsv",
                               full.names=TRUE)
cohortCellularity <- lapply(cellularityFileL, fread, sep="\t")
names(cohortCellularity) <- gsub('^.*/|_.*$', '', cellularityFileL)
cohortCellularityDT <- rbindlist(cohortCellularity)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
mergedCohortCellDT <- merge(cohortCellularityDT, annotSampMetaClassDT, by.x="sample", by.y="samples")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
boxplotCellarityByCohort(mergedCohortCellDT,
                         saveDir=file.path("..", "results"),
                         fileName="cohortCellularityBoxplot.pdf")




# -------------------------------------------------------------------------
# 11. Copy Number Abberation ----------------------------------------------
# -------------------------------------------------------------------------


## ----load_CNA_data----------------------------------------------------------------------------------------------------------------------------------------
copyNumSNPposDT <- fread(file.path("..", "data", "copyNumSNPpos.csv"))
osloTCGAcnScores <- fread(file.path('..', 'data', 'osloTCGAcnScores.csv'))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
if (file.exists(file.path(resultsDir, "CNAscoreByProbe.csv"))) {
  CNAscoreByProbeDT <- fread(file.path(resultsDir, 'CNAscoreByProbe.csv'))
} else {
  setDTthreads(14)
  CNAscoreByProbeDT <- countCNAscoresByProbe(osloTCGAcnScores, copyNumSNPposDT, annotSampMetaClassDT)
  fwrite(CNAscoreByProbeDT, file.path(resultsDir, "CNAscoreByProbe.csv"))
}

## ---------------------------------------------------------------------------------------------------------------------------------------------------------

## TODO:: Add plot saving functionality to plotPerChromosomeCNAscores
chrScorePlots <- plotPerChromosomeCNAscores(CNAscoreByProbeDT)

## ----fig.height=14, fig.width=14--------------------------------------------------------------------------------------------------------------------------
chrScorePlotGridL <- ggarrangePlotL(chrScorePlots, 3)
for (i in seq_along(chrScorePlotGridL)) {
  ggsave(chrScorePlotGridL[[i]], file=file.path(resultsDir, paste0('CNAplotGrid', i, ".pdf")))
}




# -------------------------------------------------------------------------
# 12. Calculate Cluster Meta-Effect Size ----------------------------------
# -------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
genewiseCohortsDT <- makeGenewiseCohortsDT(allProcCohorts, annotSampMetaClassDT)


## ----calculate_the_genewise_metaeffect_size---------------------------------------------------------------------------------------------------------------
# TODO:: Rewrite this with group by in data.table
if (file.exists(file.path(resultsDir, 'classByGnEffSize.csv'))) {
  classByGnEffSizeDT <- fread(file.path(resultsDir, 'classByGnEffSize.csv'))
} else {
  classByGnEffSizeDT <- calcClustMetaEstStats(genewiseCohortsDT)
  fwrite(classByGnEffSizeDT, file.path(resultsDir, 'classByGnEffSize.csv'))
}




# -------------------------------------------------------------------------
# 13. Pathway Analysis ----------------------------------------------------
# -------------------------------------------------------------------------


## ----load_oathway_data------------------------------------------------------------------------------------------------------------------------------------
geneSetL <- loadGnSetGMTs(file.path("..", "data"))
names(geneSetL)

## ----rank_each_metaclass----------------------------------------------------------------------------------------------------------------------------------
rankedMetaClassGenes <- rankMetaClassGenesByEffSize(classByGnEffSizeDT)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
referenceGenes <- rownames(cohortsDataL[[1]])
pathwayStatsDT <- computePathwayScores(rankedMetaClassGenes, geneSetL, referenceGenes)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
pathwayHMPs <- heatmapPathwayScores(pathwayStatsDT, significance=0.05)
for (i in seq_along(pathwayHMPs)) {
  ggsave(pathwayHMPs[[i]], file=file.path(resultsDir, paste0('pathwayHMP', names(geneSetL)[i], ".pdf")))
}

## ----plot_grid_heatmap_of_gene_signature_scores, fig.height=12, fig.width=10------------------------------------------------------------------------------
##TODO:: Add plot saving functionality to heatmapPathwaysScores
pathwayHMPgridL <- ggarrangePlotL(pathwayHMPs, 3)
pathwayHMPgridL




# -------------------------------------------------------------------------
# 14. Published Classifier Subtyping --------------------------------------
# -------------------------------------------------------------------------


## ---------------------------------------------------------------------------------------------------------------------------------------------------------
classifFileL <- list.files(file.path('..', 'data'), pattern="centroid",
                          ignore.case=TRUE, full.names=TRUE)

classifL <- lapply(classifFileL, readRDS)
# Format file path into the name of the classifier
names(classifL) <- gsub('^.*\\/|[Cc]entroid*|\\.[^\\.]*$', '', classifFileL)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
publishedClassifSubtypeDT <- subtypeDataLwClassifCentroidL(PGxCelllinesL, classifL,
                                                          PGxCCLpreds, seed=1987)

## ---------------------------------------------------------------------------------------------------------------------------------------------------------
classifComparison <- plotClassifierComparisons(publishedClassifSubtypeDT, allClassif=TRUE)
ggsave(classifComparison, file=file.path(resultsDir, "publishedClassifSubtypeComparison.pdf"))

## ----caclualte_association_stats--------------------------------------------------------------------------------------------------------------------------
assocStatsDT <- calcAssocStats(publishedClassifSubtypeDT)


## ----heatmap_classif_assoc_stats--------------------------------------------------------------------------------------------------------------------------
classifCorHMP <- heatmapClassifCors(assocStatsDT)
ggsave(classifCorHMP, file=file.path(resultsDir, "publishedClassifCorrHMP.pdf"))

