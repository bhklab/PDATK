# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Run PCOSP Model Script
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


## ----build_and_select_models---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(PCOSP)
library(MetaGxPancreas)

# -------------------------------------------------------------------------
# 1. PCOSP Model Building -------------------------------------------------
# -------------------------------------------------------------------------

# NOTE: The train/test splits were extracted from the original PCOSP .rda files and are assumed to be random
#       It may be worth checking for correlations between training/testing data to ensure this assumption is valid

## -- ICGC sequencing
ICGCseqTrain <- c("ICGC_0006", "ICGC_0007", "ICGC_0009", "ICGC_0020", "ICGC_0021",
                  "ICGC_0149", "ICGC_0150", "ICGC_0153", "ICGC_0169", "ICGC_0185",
                  "ICGC_0025", "ICGC_0026", "ICGC_0031", "ICGC_0033", "ICGC_0037",
                  "ICGC_0048", "ICGC_0051", "ICGC_0052", "ICGC_0053", "ICGC_0054",
                  "ICGC_0055", "ICGC_0059", "ICGC_0061", "ICGC_0063", "ICGC_0066",
                  "ICGC_0067", "ICGC_0075", "ICGC_0087", "ICGC_0088", "ICGC_0099",
                  "ICGC_0115", "ICGC_0124", "ICGC_0134", "ICGC_0135", "ICGC_0139",
                  "ICGC_0103", "ICGC_0105", "ICGC_0108", "ICGC_0109", "ICGC_0114",
                  "ICGC_0188", "ICGC_0192", "ICGC_0199", "ICGC_0201", "ICGC_0205",
                  "ICGC_0206", "ICGC_0207", "ICGC_0214", "ICGC_0223", "ICGC_0224",
                  "ICGC_0227", "ICGC_0230", "ICGC_0235", "ICGC_0295", "ICGC_0296",
                  "ICGC_0300", "ICGC_0301", "ICGC_0303", "ICGC_0304", "ICGC_0309",
                  "ICGC_0312", "ICGC_0313", "ICGC_0315", "ICGC_0321", "ICGC_0326",
                  "ICGC_0338", "ICGC_0365", "ICGC_0391", "ICGC_0392", "ICGC_0393",
                  "ICGC_0395", "ICGC_0406", "ICGC_0412", "ICGC_0415", "ICGC_0417",
                  "ICGC_0419", "ICGC_0420", "ICGC_0486", "ICGC_0518", "ICGC_0521",
                  "ICGC_0140", "ICGC_0141", "ICGC_0143", "ICGC_0144", "ICGC_0146",
                  "ICGC_0526", "ICGC_0535", "ICGC_0536", "ICGC_0543")

ICGCseqValidation <- c("ICGC_0006", "ICGC_0007", "ICGC_0009", "ICGC_0020", "ICGC_0021",
                       "ICGC_0048", "ICGC_0051", "ICGC_0052", "ICGC_0053", "ICGC_0054",
                       "ICGC_0025", "ICGC_0026", "ICGC_0031", "ICGC_0033", "ICGC_0037",
                       "ICGC_0055", "ICGC_0059", "ICGC_0061", "ICGC_0063", "ICGC_0066",
                       "ICGC_0067", "ICGC_0075", "ICGC_0087", "ICGC_0088", "ICGC_0099",
                       "ICGC_0103", "ICGC_0105", "ICGC_0108", "ICGC_0109", "ICGC_0114",
                       "ICGC_0115", "ICGC_0124", "ICGC_0134", "ICGC_0135", "ICGC_0139",
                       "ICGC_0140", "ICGC_0141", "ICGC_0143", "ICGC_0144", "ICGC_0146",
                       "ICGC_0149", "ICGC_0150", "ICGC_0153", "ICGC_0169", "ICGC_0185",
                       "ICGC_0188", "ICGC_0192", "ICGC_0199", "ICGC_0201", "ICGC_0205",
                       "ICGC_0206", "ICGC_0207", "ICGC_0212", "ICGC_0214", "ICGC_0215",
                       "ICGC_0223", "ICGC_0224", "ICGC_0227", "ICGC_0230", "ICGC_0235",
                       "ICGC_0295", "ICGC_0296", "ICGC_0300", "ICGC_0301", "ICGC_0303",
                       "ICGC_0304", "ICGC_0309", "ICGC_0312", "ICGC_0313", "ICGC_0315",
                       "ICGC_0321", "ICGC_0326", "ICGC_0338", "ICGC_0354", "ICGC_0365",
                       "ICGC_0391", "ICGC_0392", "ICGC_0393", "ICGC_0395", "ICGC_0406",
                       "ICGC_0412", "ICGC_0415", "ICGC_0417", "ICGC_0419", "ICGC_0420",
                       "ICGC_0486", "ICGC_0502", "ICGC_0507", "ICGC_0518", "ICGC_0521",
                       "ICGC_0522", "ICGC_0526", "ICGC_0535", "ICGC_0536", "ICGC_0543")

ICGCseqTrainTestSplit <- length(ICGCseqTrain) / (length(ICGCseqTrain) + length(ICGCseqValidation))
print(ICGCseqTrainTestSplit)

# -- ICGC microarray
ICGCmicroTrain <- c("ICGC_0006", "ICGC_0007", "ICGC_0009", "ICGC_0020", "ICGC_0021",
                    "ICGC_0055", "ICGC_0059", "ICGC_0061", "ICGC_0063", "ICGC_0066",
                    "ICGC_0025", "ICGC_0026", "ICGC_0031", "ICGC_0033", "ICGC_0037",
                    "ICGC_0067", "ICGC_0075", "ICGC_0087", "ICGC_0088", "ICGC_0099",
                    "ICGC_0103", "ICGC_0105", "ICGC_0108", "ICGC_0109", "ICGC_0114",
                    "ICGC_0115", "ICGC_0124", "ICGC_0134", "ICGC_0135", "ICGC_0139",
                    "ICGC_0140", "ICGC_0141", "ICGC_0143", "ICGC_0144", "ICGC_0146",
                    "ICGC_0149", "ICGC_0150", "ICGC_0153", "ICGC_0169", "ICGC_0185",
                    "ICGC_0188", "ICGC_0192", "ICGC_0199", "ICGC_0201", "ICGC_0205",
                    "ICGC_0206", "ICGC_0207", "ICGC_0214", "ICGC_0223", "ICGC_0224",
                    "ICGC_0227", "ICGC_0230", "ICGC_0235", "ICGC_0295", "ICGC_0296",
                    "ICGC_0300", "ICGC_0301", "ICGC_0303", "ICGC_0304", "ICGC_0309",
                    "ICGC_0312", "ICGC_0313", "ICGC_0315", "ICGC_0321", "ICGC_0326",
                    "ICGC_0338", "ICGC_0365", "ICGC_0391", "ICGC_0392", "ICGC_0393",
                    "ICGC_0395", "ICGC_0406", "ICGC_0412", "ICGC_0415", "ICGC_0417",
                    "ICGC_0419", "ICGC_0420", "ICGC_0486", "ICGC_0518", "ICGC_0521",
                    "ICGC_0048", "ICGC_0051", "ICGC_0052", "ICGC_0053", "ICGC_0054",
                    "ICGC_0526", "ICGC_0535", "ICGC_0536", "ICGC_0543")

ICGCmicroValidation <- c("ICGC_0001", "ICGC_0002", "ICGC_0003", "ICGC_0004", "ICGC_0005",
                         "ICGC_0019", "ICGC_0023", "ICGC_0024", "ICGC_0027", "ICGC_0028",
                         "ICGC_0029", "ICGC_0030", "ICGC_0032", "ICGC_0034", "ICGC_0035",
                         "ICGC_0036", "ICGC_0038", "ICGC_0039", "ICGC_0040", "ICGC_0041",
                         "ICGC_0042", "ICGC_0043", "ICGC_0044", "ICGC_0046", "ICGC_0047",
                         "ICGC_0049", "ICGC_0050", "ICGC_0056", "ICGC_0057", "ICGC_0058",
                         "ICGC_0060", "ICGC_0062", "ICGC_0064", "ICGC_0065", "ICGC_0069",
                         "ICGC_0071", "ICGC_0077", "ICGC_0078", "ICGC_0079", "ICGC_0080",
                         "ICGC_0089", "ICGC_0090", "ICGC_0091", "ICGC_0092", "ICGC_0093",
                         "ICGC_0094", "ICGC_0095", "ICGC_0096", "ICGC_0097", "ICGC_0101",
                         "ICGC_0102", "ICGC_0104", "ICGC_0106", "ICGC_0107", "ICGC_0110",
                         "ICGC_0111", "ICGC_0112", "ICGC_0116", "ICGC_0118", "ICGC_0119",
                         "ICGC_0122", "ICGC_0125", "ICGC_0128", "ICGC_0136", "ICGC_0145",
                         "ICGC_0147", "ICGC_0152", "ICGC_0154", "ICGC_0157", "ICGC_0158",
                         "ICGC_0159", "ICGC_0161", "ICGC_0162", "ICGC_0164", "ICGC_0166",
                         "ICGC_0167", "ICGC_0170", "ICGC_0173", "ICGC_0175", "ICGC_0177",
                         "ICGC_0178", "ICGC_0180", "ICGC_0187", "ICGC_0189", "ICGC_0190",
                         "ICGC_0194", "ICGC_0197", "ICGC_0202", "ICGC_0203", "ICGC_0208",
                         "ICGC_0209", "ICGC_0216", "ICGC_0217", "ICGC_0218", "ICGC_0219",
                         "ICGC_0225", "ICGC_0226", "ICGC_0229", "ICGC_0232", "ICGC_0233",
                         "ICGC_0234", "ICGC_0236", "ICGC_0237", "ICGC_0238", "ICGC_0298",
                         "ICGC_0299", "ICGC_0302", "ICGC_0305", "ICGC_0307", "ICGC_0308",
                         "ICGC_0310", "ICGC_0311", "ICGC_0314", "ICGC_0316", "ICGC_0317",
                         "ICGC_0319", "ICGC_0323", "ICGC_0325", "ICGC_0336", "ICGC_0337",
                         "ICGC_0339", "ICGC_0340", "ICGC_0341", "ICGC_0342", "ICGC_0343",
                         "ICGC_0346", "ICGC_0349", "ICGC_0350", "ICGC_0352", "ICGC_0356",
                         "ICGC_0359", "ICGC_0360", "ICGC_0361", "ICGC_0362", "ICGC_0363",
                         "ICGC_0364", "ICGC_0394", "ICGC_0396", "ICGC_0397", "ICGC_0398",
                         "ICGC_0400", "ICGC_0402", "ICGC_0403", "ICGC_0404", "ICGC_0405",
                         "ICGC_0407", "ICGC_0408", "ICGC_0410", "ICGC_0414", "ICGC_0416",
                         "ICGC_0418", "ICGC_0421", "ICGC_0422", "ICGC_0423", "ICGC_0487",
                         "ICGC_0488", "ICGC_0490", "ICGC_0493", "ICGC_0499", "ICGC_0514",
                         "ICGC_0523", "ICGC_0530", "ICGC_0531", "ICGC_0537", "ICGC_0538",
                         "ICGC_0014", "ICGC_0015", "ICGC_0016", "ICGC_0017", "ICGC_0018",
                         "ICGC_0008", "ICGC_0010", "ICGC_0011", "ICGC_0012", "ICGC_0013",
                         "ICGC_0540", "ICGC_0542", "ICGC_0544")

ICGCmicroTrainTestSplit <- length(ICGCmicroTrain) / (length(ICGCmicroTrain) + length(ICGCmicroValidation))
print(ICGCmicroTrainTestSplit)

load(file.path('..', 'data', 'PCOSP_trainingCohorts.rda'), verbose=TRUE)
load(file.path('..', 'data', 'PCOSP_validationCohorts.rda'), verbose=TRUE)


oldICGCdata <- mapply(rbind, trainingCohorts, validationCohorts[c('ICGC_arr', 'ICGC_seq')], SIMPLIFY=FALSE)
oldData <- c(oldICGCdata, validationCohorts[!(names(validationCohorts) %in% c("ICGC_arr", "ICGC_seq", "ICGC_array_all"))])
names(oldData)[1:2] <- c('icgcmicro', 'icgcseq')
names(oldData)[8] <- c("collison")

.multiGrepL <- function(patterns, x, ...) vapply(patterns, grepl, x=x, ..., FUN.VALUE=logical(length(x)))
nameIdxs <- unlist(apply(.multiGrepL(names(oldData), x=names(cohortSEs), ignore.case=TRUE), 2, which))
names(oldData) <- names(cohortSEs)[nameIdxs]

oldSurivalData <- lapply(oldData, function(cohort)
    data.table(cohort, keep.rownames="sample_id")[, .(sample_id, OS, OS_Status)])

pancreasData <- loadPancreasDatasets()
cohortSEs <- pancreasData$SummarizedExperiments

oldSurvivalData <- lapply(oldData, function(cohort) data.table(cohort, keep.rownames='sample_id'))
cohortSEs <- cohortSEs[names(oldSurivalData)]

updateSumExpSurvival <- function(SE, survivalDT) {
    colDataDT <- data.table(colData(SE), keep.rownames='sample_id')

    colDataDT[is.na(days_to_death)]
}






























# ----split_data-------------------------------------------------------------------------------------------------------------------


# Find common genes
commonGenes <- Reduce(intersect, lapply(cohortSEs, rownames))
cohortSEs <- lapply(cohortSEs, `[`, i=commonGenes, j=TRUE)

# -- Extract ICGC sequencing and microarray
ICGCmicro_SE <- cohortSEs$ICGCMICRO_SumExp
ICGCseq_SE <- cohortSEs$ICGCSEQ_SumExp

ICGCmicro_exprs <- data.table(t(assay(ICGCmicro_SE, 'exprs')), keep.rownames='sample_id')
ICGCseq_exprs <- data.table(t(assay(ICGCseq_SE, 'exprs')), keep.rownames='sample_id')

ICGCmicro_survival <- as.data.table(as.data.frame(colData(ICGCmicro_SE)[, c("unique_patient_ID", "days_to_death", "vital_status")]))
ICGCseq_survival <- as.data.table(as.data.frame(colData(ICGCseq_SE)[, c("unique_patient_ID", "days_to_death", "vital_status")]))

# Death to 1, Living to 0
ICGCmicro_survival[vital_status == 'deceased', vital_status := 1]
ICGCmicro_survival[vital_status == 'living', vital_status := 0]
ICGCseq_survival[vital_status == 'deceased', vital_status := 1]
ICGCseq_survival[vital_status == 'living', vital_status := 0]

# -- Updating missing survival data from Vandanas old data

# Read in old data and make into data.table
load(file.path('..', 'data', 'PCOSP_trainingCohorts.rda'), verbose=TRUE)
load(file.path('..', 'data', 'PCOSP_validationCohorts.rda'), verbose=TRUE)

oldMicroSurvival <- rbind(data.table(trainingCohorts$icgc_array_cohort, keep.rownames=sample_id)[, .(rn, OS, OS_Status)],
                  data.table(validationCohorts$ICGC_arr, keep.rownames=sample_id)[, .(rn, OS, OS_Status)])

.factorToNumeric <- function(col) as.numeric(levels(col)[col])

oldMicroSurvival[, `:=`(OS = .factorToNumeric(OS), OS_Status = .factorToNumeric(OS_Status))]

oldSeqSurvival <- rbind(data.table(trainingCohorts$icgc_seq_cohort, keep.rownames=TRUE)[, .(rn, OS, OS_Status)],
                  data.table(validationCohorts$ICGC_seq, keep.rownames=TRUE)[, .(rn, OS, OS_Status)])
oldSeqSurvival[, `:=`(OS = .factorToNumeric(OS), OS_Status = .factorToNumeric(OS_Status))]

# Rename columns
colnames(ICGCmicro_survival) <- c("sample_id", "OS", "OS_Status")
colnames(ICGCseq_survival) <- c("sample_id", "OS", "OS_Status")

# Fill missing values in OS survival data
ICGCmicro_survival[is.na(OS),
                   `:=`(OS = oldMicroSurvival[rn %in% sample_id][order(rn)]$OS,
                        OS_Status = as.numeric(oldMicroSurvival[rn %in% sample_id][order(rn)]$OS_Status))]

ICGCseq_survival[is.na(OS),
                 `:=`(OS = oldSeqSurvival[rn %in% sample_id][order(rn)]$OS,
                      OS_Status = as.numeric(oldSeqSurvival[rn %in% sample_id][order(rn)]$OS_Status))]

# Update original SEs
colData(ICGCmicro_SE)$days_to_death <- ICGCmicro_survival$OS
colData(ICGCmicro_SE)$vital_status <- ICGCmicro_survival$OS_Status
colData(ICGCseq_SE)$days_to_death <- ICGCseq_survival$OS
colData(ICGCseq_SE)$vital_status <- ICGCseq_survival$OS_Status

## TODO:: Heewon - save these as the new version of ICGCSEQ and ICGCMICRO

# Attach expr and survival data as a data.frame
ICGCmicro <- merge(ICGCmicro_exprs, ICGCmicro_survival, by='sample_id')
ICGCmicro <- data.frame(ICGCmicro[, -'sample_id'], row.names=ICGCmicro$sample_id)
ICGCseq <- merge(ICGCseq_exprs, ICGCseq_survival, by='sample_id')
ICGCseq <- data.frame(ICGCseq[, -'sample_id'], row.names=ICGCseq$sample_id)

# Remove extraneous data to save RAM
rm(ICGCmicro_survival, ICGCseq_survival, ICGCseq_SE, ICGCmicro_SE, ICGCseq_exprs, ICGCmicro_exprs); gc()

# Split the data
ICGCmicro_train <- ICGCmicro[ICGCmicroTrain, ]
ICGCmicro_validation <- ICGCmicro[ICGCmicroValidation, ]
ICGCseq_train <- ICGCseq[ICGCmicroTrain, ]
ICGCseq_validation <- ICGCseq[ICGCmicroValidation, ]

#
trainingCohorts <- list("icgc_array_cohort"=ICGCmicro_train, 'icgc_seq_cohort'=ICGCseq_train)

`%notin%` <- Negate(`%in%`)
valCohorts <- cohortSEs[names(cohortSEs) %notin% c("ICGCMICRO_SumExp", "ICGCSEQ_SumExp")]
cohortsWithSurvival <- paste0(c("TCGA", "PCSI", "KIRBY", "OUH", "UNC", "CHEN", "ZHANG", "COLLISSON"), '_SumExp')
valCohorts <- valCohorts[names(valCohorts) %in% cohortsWithSurvival]
validationCohorts <- lapply(calCohorts, )

## ----train_model---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Initialize a directory for result data
resultsDir <- file.path("..", "results")

# Build the models using the selected random seed
system.time({
selectedModels <- buildPCOSPmodels(trainingCohorts, numModels=1000, nthread=15)
})

saveRDS(selectedModels, file=file.path(resultsDir, "1_selectedModels.rds"))


# -------------------------------------------------------------------------
# 2. Calculate and Validate Meta-Estimate Statistics ----------------------
# -------------------------------------------------------------------------

# Keep a copy of all the cohorts for later
allValidationCohorts <- validationCohorts

# Select the validation cohorts we want to include
include <- c("TCGA", "PCSI", "Kirby", "ICGC_arr", "UNC", "Chen", "OUH",
             "Zhang", "Winter", "Collisson")

# Subset to selected cohorts
validationCohorts <- validationCohorts[include]

# Fix names to look nice
names(validationCohorts) <- gsub("_arr", "-array", names(validationCohorts))

# Identify which cohorts contain sequencing data
seqCohorts <- c("PCSI", "TCGA", "Kirby")


## ----calculate_validation_stats------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
validationStats <- calculateValidationStats(validationCohorts,
                                            selectedModels,
                                            seqCohorts=seqCohorts,
                                            nthread=15)

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
set.seed(1987, sample.kind="Rounding")

system.time({
    randomLabelModels <- buildRandomLabelShufflingModel(trainingCohorts,
                                                    numModels=1000, nthread=15)
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
set.seed(1987, sample.kind="Rounding")

## These results will differ from the Code Ocean capsule due to sampling.
## For the original results use: original = TRUE; this will be about 5-10x slower
##     than the new method if multiple CPU cores are available.
system.time({
    randomGeneModels <- buildRandomGeneAssignmentModels(trainingCohorts,
                                                        numModels=1000,
                                                        nthread=15,
                                                        original=TRUE)
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

dataDir <- file.path("..", "data")
pathwayDataDir <- file.path(dataDir, "Pathway_db_files")
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
load(file.path(dataDir, "clinicalFeatures.rda"))
load(file.path(dataDir, "cohortClasses.rda"))

# Rename clinical features for plots and reorder them
names(clinicalFeatures) <- c("PCSI", "ICGCclinical", "TCGA", "ICGCarray", "OUH")
clinicalFeatures <- clinicalFeatures[c("TCGA", "PCSI", "ICGCarray", "OUH", "ICGCclinical")]

fitClinicalModels <- summarizeClinicalModels(clinicalFeatures, cohorts=c("ICGCclinical"))

ICGCclinicalModel <- fitClinicalModels$ICGCclinical

saveRDS(fitClinicalModels, file=file.path(resultsDir, "6a_fitClinicalModels.rds"))
saveRDS(ICGCclinicalModel, file=file.path(resultsDir, "6b_ICGCclinicalModel.rds"))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Rename cohortClasses to match cohots
rnCohortClasses <- cohortClasses
names(rnCohortClasses) <- c("PCSI", "TCGA", "ICGC", "Kirby", "ICGCarray", "UNC", "Chen", "OUH", "Winter", "Zhang", "Collison", "ICGCarrayAll")

cohorts <- c("TCGA", "PCSI", "ICGCarray", "OUH")

modelComparisonStats <- compareClinicalModels(clinicalFeatures, rnCohortClasses,
                                   models="ICGCclinical",
                      cohorts=cohorts)
saveRDS(modelComparisonStats, file=file.path(resultsDir, "6c_clinicalModelComparisonStats.rds"))


## ----barplot_model_comparison--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
barplotModelComparison(modelComparisonStats,
                       model="ICGCclinical",
                       names=c("TCGA","PCSI","ICGCarray","OUH"),
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
dataDir <- file.path("..", "data")
classifierDir <- file.path(dataDir, "Classifier_Comparison")

# Haider data (provided by author)
load(file.path(dataDir, "haiderSigScores.rda"))

# Harmonize case/format of names to use for subsetting the list
names(haiderSigScores) <- gsub("-.*|_.*", "", tolower(names(haiderSigScores)))


# Birnbaum data
fileNames1 <- list.files(file.path(classifierDir, "Birnbaum"))[-1]
geneCoefFile1 <-  list.files(file.path(classifierDir, "Birnbaum"))[1]

birnbaumSigScores <- readGenefuFilesToSigScores(file.path(classifierDir, "Birnbaum"),
                                     fileNames1,
                                     geneCoefFile1)
saveRDS(birnbaumSigScores, file=file.path(resultsDir, "8a_birnbaumSigScores.rds"))

# Chen data
fileNames2 <- list.files(file.path(classifierDir, "Chen"))[-1]
geneCoefFile2 <- list.files(file.path(classifierDir, "Chen"))[1]

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
validationCohorts <- lapply(validationCohorts, function(cohort) as.data.frame(cohort, row.names=rownames(cohort)))

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

cohortSurvivalData <- cohortSurvivalData[names(haiderSigScores)]

# Put our classifier data into a list
metaestimateData <- list("haider"=haiderSigScores, "chen"=chenSigScores, "birnbaum"=birnbaumSigScores,
                         "pcosp"=PCOSPscores, "survival"=cohortSurvivalData)

# Subset cohorts/samples, reorder to match
metaestimateData <- subsetSharedCohortsAndSamples(metaestimateData)
saveRDS(metaestimateData, file=file.path(resultDir, "9a_metaestimateDataClassifier.rds"))

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
allValidationCohorts <- lapply(allValidationCohorts, function(cohort) as.data.frame(cohort, row.names=rownames(cohort)))

# Match names with cohortClasses and reorder cohorts to match between class and cohort data
names(allValidationCohorts) <- tolower(names(allValidationCohorts))
names(allValidationCohorts)[which(names(allValidationCohorts) == "icgc_seq")] <- "icgc"

# Remove ICGC_array_all
allValidationCohorts <- allValidationCohorts[!(names(allValidationCohorts) %in% c("icgc_array_all", "icgc"))]
cohortClasses <- cohortClasses[!(names(cohortClasses) %in% c("icgc_array_all", "icgc"))]

# Reorder validation cohrots to same order as the forest plot names
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


