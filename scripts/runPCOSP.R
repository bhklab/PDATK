library(PDATK)
library(qs)

# Initialize a directory for result data
resultsDir <- file.path("results", "PCOSP")

if (!file.exists(resultsDir)) {
    dir.create(resultsDir)
    dir.create(file.path(resultsDir, 'figures'))
}

# load the package data
cohortList <- qread('data/cohort_list.qs')

# subset to only shared genes in all cohorts
commonGenes <- findCommonGenes(cohortList)
cohortList <- subset(cohortList, subset=commonGenes)

# label cohort molecular datatypes

# extract only the ICGC cohorts for initial model training/validation
ICGCcohortList <- cohortList[grepl('icgc', names(cohortList),
    ignore.case=TRUE)]

# extract the evaluation cohorts (i.e., not used to train/validate the model)
validationCohortList <- cohortList[!grepl('icgc', names(cohortList),
    ignore.case=TRUE)]

# find samples with both types of molecular data
## TODO:: Is it confusing to use subset to subset list items? Maybe make
##   generic subsetItems or subsubset?
commonSamples <- findCommonSamples(ICGCcohortList)
ICGCtrainCohorts <- subset(ICGCcohortList, select=commonSamples)
ICGCtestCohorts <- subset(ICGCcohortList, select=commonSamples, invert=TRUE)

# Drop early deaths, which we attribute to disease severity and thus aren't
#   usefulf or classifying survival
ICGCtrainCohorts <- dropNotCensored(ICGCtrainCohorts)

# Merge the training cohorts into a single SurvivalExperiment, with each
#   mDataType as an assay.
ICGCtrain <- merge(ICGCtrainCohorts[[1]], ICGCtrainCohorts[[2]],
    cohortNames=names(ICGCcohortList))

# construct a PCOSP model object
PCOSPmodel <- PCOSP(ICGCtrain, randomSeed=1987)

PCOSPmodel <- trainModel(PCOSPmodel, numModels=10, minAccuracy=0.6)

saveRDS(PCOSPmodel, file=file.path(resultsDir, "1_PCOSPmodel.rds"))


# -------------------------------------------------------------------------
# 2. Calculate and Validate Meta-Estimate Statistics ----------------------
# -------------------------------------------------------------------------

# add the ICGC testing data to the model
validationCohortList <- c(ICGCtestCohorts, validationCohortList)

predictionCohortList <- predictClasses(validationCohortList, model=PCOSPmodel)

validatedPCOSPmodel <- validateModel(PCOSPmodel, predictionCohortList)







## ----calculate_validation_stats------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
validationStats <- calculateValidationStats(validationCohorts,
                                            selectedModels,
                                            seqCohorts=seqCohorts,
                                            nthread=14)

saveRDS(validationStats, file=file.path(resultsDir, "2a_validationStats.rds"))

### forestPlotMetaEstimates
# -------------------------------------------------------------------------


## ----forest_plots--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
forestPlotMetaEstimate(validationStats,
                       stat="dIndex",
                       isSummary=c(rep(FALSE, 10), rep(TRUE, 3)),
                       filePath=file.path(resultsDir, "figures"),
                       fileName="2b_forestPlotDindexPCOSP.pdf")


## ----forest_plot_cIndex--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
forestPlotMetaEstimate(validationStats,
                       stat="cIndex",
                       isSummary=c(rep(FALSE, 10), rep(TRUE, 3)),
                       filePath=file.path(resultsDir, "figures"),
                       fileName="2c_forestPlotCindexPCOSP.pdf")

### Calculate AUC Distribution and Plot ROC Curves
# -------------------------------------------------------------------------


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
aucStats <- aucMetaEstimates(validationCohorts, validationStats,
                            seqCohorts=seqCohorts)
saveRDS(aucStats, file=file.path(resultsDir, "2c_aucStatsPCOSP.rds"))


## ----format_val_cohorts--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
formattedValCohorts <- formatValidationCohorts(validationCohorts)

saveRDS(formattedValCohorts, file=file.path(resultsDir, "2e_formattedValiCohorts.rds"))


## ----get_PCOSP_scores----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
PCOSPscores <- validationStats$probabilities


## ----plot_PCOSP_roc_curve------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
colours <- c("chartreuse3", "magenta", "#fb9a99", "turquoise3", "darkgoldenrod1",
             "wheat4", "green", "red", "cornflowerblue", "mediumorchid2")

## FIXME:: Find out where warnings/messages are coming from and suppress them
plotROCcurves(formattedValCohorts,
              PCOSPscores,
              colours,
              filePath=file.path(resultsDir, "figures"),
              fileName="2f_PCOSProcCurvesPlot")


# -------------------------------------------------------------------------
# 3. Model Random Label Shuffling -----------------------------------------
# -------------------------------------------------------------------------


## ----3_model_random_label_shuffling--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

system.time({
    randomLabelModels <- buildRandomLabelShufflingModel(trainingCohorts,
                                                        seqCohort="ICGCSEQ",
                                                    numModels=1000, nthread=14)
})
saveRDS(randomLabelModels, file=file.path(resultsDir, "3a_randomLabelModels.rds"))


## ----density_plot--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
densityPlotModel(formattedValCohorts,
                 selectedModels=randomLabelModels,
                 seqCohorts,
                 vlines=c(0.69, 0.72, 0.70),
                 nthread=15,
                 title="Random reshuffling of labels",
                 filePath=file.path(resultsDir, "figures"),
                 fileName="3b_densityPlotRandomLabelShuffling.pdf")



# -------------------------------------------------------------------------
# 4. Random Gene Assignment Model -----------------------------------------
# -------------------------------------------------------------------------

## ----4_random_gene_assignment_models-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## These results will differ from the Code Ocean capsule due to sampling.
## For the original results use: original = TRUE; this will be about 5-10x slower
##     than the new method if multiple CPU cores are available.
system.time({
    randomGeneModels <- buildRandomGeneAssignmentModels(trainingCohorts,
                                                        seqCohort="ICGCSEQ",
                                                        numModels=1000,
                                                        nthread=15,
                                                        original=FALSE)
})
saveRDS(randomGeneModels, file=file.path(resultsDir, "4a_randomGeneAssignModel.rds"))


## ----random_gene_assignment_density--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
densityPlotModel(formattedValCohorts,
                 selectedModels=randomGeneModels,
                 seqCohorts,
                 vlines=c(0.69, 0.72, 0.70),
                 nthread=15,
                 title="Random gene assignment",
                 filePath=file.path(resultsDir, "figures"),
                 fileName="4b_densityPlotRandomGeneAssignment.pdf")


# -------------------------------------------------------------------------
# 5. HyperGSA Pathway Analysis --------------------------------------------
# -------------------------------------------------------------------------


## ----hyper_gsea_models---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

dataDir <- file.path("..", "data", "PCOSP")
pathwayDataDir <- file.path(dataDir, "PCOSP_Pathway_db_files")
referenceGenes <- read.table(file.path(pathwayDataDir, "reference_genes.txt"))[, 1]
PCOSPgenes <- read.table(file.path(pathwayDataDir, "pcosp_genes.txt"))[, 1]

GSCfiles <- c("hallmark.h.all.v6.1.symbols.gmt",
              "canonical.c2.cp.v6.1.symbols.gmt",
              "GOmolecular.c5.mf.v6.1.symbols.gmt",
              "GOcellular.c5.cc.v6.1.symbols.gmt")

GSCpaths <- vapply(GSCfiles,
                   function(file, pathwayDataDir) file.path(pathwayDataDir, file),
                   pathwayDataDir,
                   FUN.VALUE=character(1))

GSAmodels <- lapply(GSCpaths,
                  function(filePath, PCOSPgenes, referenceGenes, adjMethod)
                      fitGSAtoGeneSetCollection(filePath=filePath,
                                                 PCOSPgenes=PCOSPgenes,
                                                 referenceGenes=referenceGenes,
                                                 adjMethod=adjMethod),
                  PCOSPgenes=PCOSPgenes,
                  referenceGenes=referenceGenes,
                  adjMethod="fdr"
                  )
names(GSAmodels) <- gsub("\\..*", "", GSCfiles)
saveRDS(GSAmodels, file=file.path(resultsDir, "5_GSAmodels.rds"))


# -------------------------------------------------------------------------
# 6. Clinical Model Comparison --------------------------------------------
# -------------------------------------------------------------------------


## ----clinical_model_comparison-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load(file.path(dataDir, "PCOSP_clinicalFeatures.rda"))
load(file.path(dataDir, "PCOSP_cohortClasses.rda"))

# Rename clinical features for plots and reorder them
names(clinicalFeatures) <- c("PCSI", "ICGCclinical", "TCGA", "ICGCMICRO", "OUH")
clinicalFeatures <- clinicalFeatures[c("TCGA", "PCSI", "ICGCMICRO", "OUH", "ICGCclinical")]

fitClinicalModels <- summarizeClinicalModels(clinicalFeatures, cohorts=c("ICGCclinical"))

ICGCclinicalModel <- fitClinicalModels$ICGCclinical

saveRDS(fitClinicalModels, file=file.path(resultsDir, "6a_fitClinicalModels.rds"))
saveRDS(ICGCclinicalModel, file=file.path(resultsDir, "6b_ICGCclinicalModel.rds"))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Rename cohortClasses to match cohots
rnCohortClasses <- cohortClasses
names(rnCohortClasses) <- c("PCSI", "TCGA", "ICGC", "KIRBY", "ICGCMICRO", "UNC",
                            "CHEN", "OUH", "WINTER", "ZHANG", "COLLISON",
                            "ICGCMICRO_ALL")
rnCohortClasses <- rnCohortClasses[c("TCGA", "PCSI", "KIRBY", "ICGCMICRO",
                                     "UNC", "CHEN", "OUH", "ZHANG",
                                     "WINTER", "COLLISON")]

cohorts <- c("TCGA", "PCSI", "ICGCMICRO", "OUH")

modelComparisonStats <- compareClinicalModels(clinicalFeatures, rnCohortClasses,
                                   models="ICGCclinical",
                      cohorts=cohorts)
saveRDS(modelComparisonStats, file=file.path(resultsDir, "6c_clinicalModelComparisonStats.rds"))


## ----barplot_model_comparison--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
barplotModelComparison(modelComparisonStats,
                       model="ICGCclinical",
                       names=c("TCGA","PCSI","ICGCMICRO","OUH"),
                       colours=c("palevioletred1", "darkgrey"),
                       filePath=file.path(resultsDir, "figures"),
                       fileName="6d_modelComparisonBarPlot.pdf")


## ----calculate_model_comparison_pvals------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
modelComparisonPvals <- metaEstimateComparisonAUCs(modelComparisonStats,
                                                   model="ICGCclinical",
                                                   seqCohorts=c("PCSI", "TCGA"))
modelComparisonPvals
saveRDS(modelComparisonPvals, file=file.path(resultsDir, "6e_modelComparisonPvals.rds"))


# -------------------------------------------------------------------------
# 7. Clinical Model Dindex and Cindex -------------------------------------
# -------------------------------------------------------------------------


## ----clinical_model_dindex_cindex----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

clinicalFeatures <- clinicalFeatures[grepl(paste(cohorts, collapse="|"), names(clinicalFeatures))]
names(clinicalFeatures) <- cohorts

modelProbabilities <- calculateCohortProbabilties(fitClinicalModels, clinicalFeatures)

clinicalModelStats <- calculateModelDandCindex(modelProbabilities,
                                               clinicalFeatures,
                                               seqCohorts=c("PCSI", "TCGA"),
                                               model="ICGCclinical")
saveRDS(clinicalModelStats, file=file.path(resultsDir, "7a_clinicalModelStats.rds"))


## ----model_comparison_forest_plot_dindex, fig.height=11, fig.width=8.5---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
forestPlotModelComparison(clinicalModelStats,
                           stat="dIndex",
                           isSummary=c(rep(FALSE, 4), rep(TRUE, 3)),
                           filePath=file.path(resultsDir, "figures"),
                           fileName="7b_modelComparisonDindexForestPlot.pdf")


## ----model_comparison_forest_plot_cindex, fig.height=11, fig.width=8.5---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
forestPlotModelComparison(clinicalModelStats,
                           stat="cIndex",
                           isSummary=c(rep(FALSE, 4), rep(TRUE, 3)),
                           filePath=file.path(resultsDir, "figures"),
                           fileName="7c_modelComparisonCindexForestPlot.pdf")



# -------------------------------------------------------------------------
# 8. Existing Classifier Score Calculations -------------------------------
# -------------------------------------------------------------------------


## ----existing_classifier_forest_plots------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dataDir <- file.path("..", "data", "PCOSP")
classifierDir <- file.path(dataDir, "PCOSP_Classifier_Comparison")

# Haider data (provided by author)
load(file.path(dataDir, "PCOSP_haiderSigScores.rda"))

# Harmonize case/format of names to use for subsetting the list
names(haiderSigScores) <- gsub("-.*|_.*", "", tolower(names(haiderSigScores)))


# Birnbaum data
fileNames1 <- list.files(file.path(classifierDir, "Birnbaum"))[-1]
geneCoefFile1 <-  list.files(file.path(classifierDir, "Birnbaum"))[1] # bmc_genes.txt

birnbaumSigScores <- readGenefuFilesToSigScores(file.path(classifierDir, "Birnbaum"),
                                     fileNames1,
                                     geneCoefFile1)
saveRDS(birnbaumSigScores, file=file.path(resultsDir, "8a_birnbaumSigScores.rds"))

# Chen data
fileNames2 <- list.files(file.path(classifierDir, "Chen"))[-1]
geneCoefFile2 <- list.files(file.path(classifierDir, "Chen"))[1] # chen_genes.txt

chenSigScores <- readGenefuFilesToSigScores(file.path(classifierDir, "Chen"),
                                            fileNames2,
                                            geneCoefFile2,
                                            signed=TRUE)

saveRDS(chenSigScores, file=file.path(resultsDir, "8b_chenSigScores.rds"))

# PCOSP data
PCOSPscores <- validationStats$probabilities

# Harmonize case/format of names to use for subsetting the list
names(PCOSPscores) <- gsub("-.*|_.*", "", tolower(names(PCOSPscores)))
saveRDS(PCOSPscores, file=file.path(resultsDir, "8c_PCOSPscores.rds"))


# -------------------------------------------------------------------------
# 9. Existing Classifier Forestplots --------------------------------------
# -------------------------------------------------------------------------



## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
validationCohorts <- lapply(validationCohorts,
                            function(cohort)
                                as.data.frame(cohort,
                                              row.names=rownames(cohort)))

as.numeric.factor <- function(factor) {
    if (is.factor(factor)) as.numeric(levels(factor)[factor]) else as.numeric(factor)
}

# Get cohort survival data
cohortSurvivalData <- lapply(validationCohorts,
                             function(cohort) {
                                 OS <- as.numeric.factor(cohort$OS)
                                 OS_Status <- as.numeric.factor(cohort$OS_Status)
                                 data.frame("OS"=OS, "OS_Status"=OS_Status,
                                     row.names=rownames(cohort))
                                 })

# Fix TCGA data to match study files?
PCOSPscores$tcga <- PCOSPscores$tcga[cohortSurvivalData$TCGA$OS > 180]
cohortSurvivalData$TCGA <- cohortSurvivalData$TCGA[cohortSurvivalData$TCGA$OS > 180, ]

# Harmonize case/format of names to use for subsetting the list
names(cohortSurvivalData) <- gsub("-.*|_.*", "", tolower(names(cohortSurvivalData)))
names(cohortSurvivalData)[names(cohortSurvivalData) == "icgcmicro"] <- 'icgc'
names(cohortSurvivalData)[names(cohortSurvivalData) == "collison"] <- 'collisson'

cohortSurvivalData <- cohortSurvivalData[names(haiderSigScores)]

# Put our classifier data into a list
metaestimateData <- list("haider"=haiderSigScores, "chen"=chenSigScores, "birnbaum"=birnbaumSigScores,
                         "pcosp"=PCOSPscores, "survival"=cohortSurvivalData)

# Subset cohorts/samples, reorder to match
metaestimateData <- subsetSharedCohortsAndSamples(metaestimateData)
saveRDS(metaestimateData, file=file.path(resultsDir, "9a_metaestimateDataClassifier.rds"))

seqCohorts <- c("tcga", "pcsi", "kirby")

classifierStats <- lapply(metaestimateData[names(metaestimateData) != "survival"],
                          constructMetaEstimatesDF,
                          cohortData=metaestimateData[["survival"]],
                          seqCohorts=seqCohorts,
                          hetero=rep(TRUE, 3))
saveRDS(classifierStats, file=file.path(resultsDir, "9b_classifierStats.rds"))


## ----classifier_forest_plot_c_index, fig.height=11, fig.width=8.5--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
forestPlotClassifierModelComparision(classifierStats, "cIndex",
                        names=c("Haider, 2014", "Chen, 2015", "Birnbaum, 2017", "PCOSP"),
                        filePath=file.path(resultsDir, "figures"),
                        fileName="9c_classifierComparisonForestPlot.pdf")


## ----classifier_forest_plot_d_index, fig.height=11, fig.width=8.5--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
forestPlotClassifierModelComparision(classifierStats, "dIndex",
                                names=c("Haider, 2014", "Chen, 2015", "Birnbaum, 2017", "PCOSP"),
                                filePath=file.path(resultsDir, "figures"),
                                fileName="9c_classifierComparisonForestPlot.pdf")


# -------------------------------------------------------------------------
# 10. Subtype Specific Meta-Estimates -------------------------------------
# -------------------------------------------------------------------------


## ----subtype_statistics--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Ensure all cohorts are data.frames
allValidationCohorts <- lapply(allValidationCohorts,
                               function(cohort)
                                   as.data.frame(cohort,
                                                 row.names=rownames(cohort)))

# Match names with cohortClasses and reorder cohorts to match between class and cohort data
names(allValidationCohorts) <- tolower(names(allValidationCohorts))
names(allValidationCohorts)[which(names(allValidationCohorts) == "icgcmicro")] <- "icgc_arr"
names(allValidationCohorts)[which(names(allValidationCohorts) == "collison")] <- "collisson"

# Remove ICGC_array_all
allValidationCohorts <- allValidationCohorts[!(names(allValidationCohorts) %in% c("icgc_array_all", "icgc"))]
cohortClasses <- cohortClasses[!(names(cohortClasses) %in% c("icgc_array_all", "icgc"))]

# Reorder validation cohorts to same order as the forest plot names
plotNames <- c("tcga", "pcsi", "kirby", "icgc_arr", "unc", "chen", "ouh", "zhang",  "winter", "collisson")
cohortClasses <- cohortClasses[plotNames]
allValidationCohorts <- allValidationCohorts[plotNames]

saveRDS(allValidationCohorts, file=file.path(resultsDir, "10a_allValidationCohorts.rds"))

# Get cohort survival data
allCohortSurvivalData <- lapply(allValidationCohorts,
                             function(cohort) {
                                 OS <- as.numeric.factor(cohort$OS)
                                 OS_Status <- as.numeric.factor(cohort$OS_Status)
                                 data.frame("OS"=OS, "OS_Status"=OS_Status,
                                 row.names=rownames(cohort))
                                 })
saveRDS(allCohortSurvivalData, file=file.path(resultsDir, "10b_allCohortSurvivalData.rds"))

# Get probabilties for all samples
allCohortProbabilites <- calculatePCOSPscores(allValidationCohorts, selectedModels, nthread=15)
saveRDS(allCohortProbabilites, file=file.path(resultsDir, "10c_allCohortProbabilities.rds"))

# Set sequencing cohorts
seqCohorts <- c("tcga", "pcsi", "kirby")

# Calculate meta-estimates between subtypes per cohort
cohortSubtypeStats <- calculateCohortSubtypeStats(allCohortProbabilites, allCohortSurvivalData, cohortClasses, seqCohorts=seqCohorts)
saveRDS(cohortSubtypeStats, file=file.path(resultsDir, "10d_cohortSubtypeStats.rds"))


## ----forest_plot_subtype_comparison_dindex, fig.height=11, fig.width=8.5-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
forestPlotCohortSubtypeComparison(cohortSubtypeStats, "dIndex",
                                  cohortNames=c(NA, "TCGA (n=146)", "Basal-subtype (n=29)", "Classical-subtype (n=117)",
                                                NA, "PCSI (n=118)", "Basal-subtype (n=25)", "Classical-subtype (n=93)",
                                                NA, "Kirby (n=51)", "Basal-subtype (n=13)", "Classical-subtype (n=38)",
                                                NA, "ICGC-array (n=178)", "Basal-subtype (n=156)", "Classical-subtype (n=22)",
                                                NA, "UNC (n=125)", "Basal-subtype (n=27)", "Classical-subtype (n=98)",
                                                NA, "Chen (n=63)", "Basal-subtype (n=20)", "Classical-subtype (n=43)",
                                                NA, "OUH (n=48)", "Basal-subtype (n=08)", "Classical-subtype (n=40)",
                                                NA, "Zhang (n=42)", "Basal-subtype (n=13)", "Classical-subtype (n=29)",
                                                NA, "Winter (n=30)", "Basal-subtype (n=8)", "Classical-subtype (n=22)",
                                                NA, "Collisson (n=27)", "Basal-subtype (n=1)", "Classical-subtype (n=26)",
                                                NA),
                                  summaryNames=c(NA,"Sequencing", "Basal-subtype", "Classical-subtype",
                                                 NA, "Microarray", "Basal-subtype", "Classical-subtype",
                                                 NA, "Overall", "Basal-subtype", "Classical-subtype"),
                                  filePath=file.path(resultsDir, "figures"),
                                  fileName="10e_subtypeDindexComparisonForestPlot.pdf"
                                  )


## ----forest_plot_subtype_comparison_cindex, fig.height=11, fig.width=8.5-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
forestPlotCohortSubtypeComparison(cohortSubtypeStats, "cIndex",
                                  cohortNames=c(NA, "TCGA (n=146)", "Basal-subtype (n=29)", "Classical-subtype (n=117)",
                                                NA, "PCSI (n=118)", "Basal-subtype (n=25)", "Classical-subtype (n=93)",
                                                NA, "Kirby (n=51)", "Basal-subtype (n=13)", "Classical-subtype (n=38)",
                                                NA, "ICGC-array (n=178)", "Basal-subtype (n=156)", "Classical-subtype (n=22)",
                                                NA, "UNC (n=125)", "Basal-subtype (n=27)", "Classical-subtype (n=98)",
                                                NA, "Chen (n=63)", "Basal-subtype (n=20)", "Classical-subtype (n=43)",
                                                NA, "OUH (n=48)", "Basal-subtype (n=08)", "Classical-subtype (n=40)",
                                                NA, "Zhang (n=42)", "Basal-subtype (n=13)", "Classical-subtype (n=29)",
                                                NA, "Winter (n=30)", "Basal-subtype (n=8)", "Classical-subtype (n=22)",
                                                NA, "Collisson (n=27)", "Basal-subtype (n=1)", "Classical-subtype (n=26)",
                                                NA),
                                  summaryNames=c(NA,"Sequencing", "Basal-subtype", "Classical-subtype",
                                                 NA, "Microarray", "Basal-subtype", "Classical-subtype",
                                                 NA, "Overall", "Basal-subtype", "Classical-subtype"),
                                  filePath=file.path(resultsDir, "figures"),
                                  fileName="10f_subtypeCindexComparisonForestPlot.pdf"
                                  )


# -------------------------------------------------------------------------
# END ---------------------------------------------------------------------
# -------------------------------------------------------------------------


