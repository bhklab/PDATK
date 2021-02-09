library(testthat)
library(PDATK)

library(testthat)
library(PDATK)

data(sampleICGCmicro)
data(sampleCohortList)

# suppressWarnings({
#     clinicalModel <- ClinicalModel(sampleICGCmicro,
#         formula='prognosis ~ sex + age + T + N + M + grade')
#     trainedClincialModel <- trainModel(clinicalModel, numModels=10)
#     clinicalPredCohortList <- predictClasses(sampleCohortList[seq_len(2)],
#         model=trainedClincialModel)
#     validatedClinicalModel <- validateModel(trainedClincialModel,
#         valData=clinicalPredCohortList)
# })
#
# test_that('trainModel method for PCOSP works correctly', {
#     expect_s4_class(trainModel(clinicalModel), 'ClinicalModel')
#     expect_true(length(models(trainedPCOSPmodel)) > 5)
#     # the model list has the correct items
#     expect_equal(names(models(trainedPCOSPmodel)[[1]]),
#         c('name', 'TSPs', 'score', 'labels'))
#     # added correct metadata during training
#     expect_equal(names(metadata(trainedPCOSPmodel)$modelParams),
#         c('randomSeed', 'RNGkind', 'minDaysSurvived', 'numModels',
#             'minAccuracy'))
# })
#
# test_that('predictClasses method for PCOSP works correctly', {
#     expect_warning(predictClasses(sampleCohortList[[1]], model=trainedPCOSPmodel),
#         '.*Dropped samples with NA survival data.*')
#     expect_true(all(mcols(PCOSPpredCohortList)$hasPredictions))
#     expect_true('PCOSP_prob_good' %in%
#         colnames(colData(PCOSPpredCohortList[[1]])))
#     expect_true(is.matrix(metadata(PCOSPpredCohortList[[1]])$PCOSPpredictions))
# })
#
# test_that('validateModel method for PCOSP  works correctly', {
#     expect_message(validateModel(trainedPCOSPmodel,
#         valData=PCOSPpredCohortList[[1]]), 'Setting levels: .*')
#     expect_true(metadata(validatedPCOSPModel)$isValidated)
#     expect_true(is.data.frame(validationStats(validatedPCOSPModel)))
#     expect_equal(colnames(validationStats(validatedPCOSPModel)),
#         c('statistic', 'estimate', 'se', 'lower', 'upper', 'p_value',
#             'n', 'isSummary', 'cohort', 'mDataType'))
#     expect_equal(validationData(validatedPCOSPModel), PCOSPpredCohortList)
# })