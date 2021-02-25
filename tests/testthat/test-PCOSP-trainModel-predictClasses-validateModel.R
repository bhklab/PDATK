library(testthat)
library(PDATK)
library(BiocParallel)

if (Sys.info()['sysname'] == 'Windows') {
    BiocParallel::register(BiocParallel::SerialParam())
}

data(sampleICGCmicro)
data(sampleCohortList)

suppressWarnings({
    PCOSPmodel <- PCOSP(sampleICGCmicro, randomSeed=1987)
    trainedPCOSPmodel <- trainModel(PCOSPmodel, numModels=5)
    PCOSPpredCohortList <- predictClasses(sampleCohortList[seq_len(2)],
        model=trainedPCOSPmodel)
    validatedPCOSPModel <- validateModel(trainedPCOSPmodel,
        valData=PCOSPpredCohortList)
})

test_that('trainModel method for PCOSP works correctly', {
    expect_s4_class(trainModel(PCOSPmodel, numModels=5), 'PCOSP')
    expect_true(length(models(trainedPCOSPmodel)) > 2)
    # the model list has the correct items
    expect_equal(names(models(trainedPCOSPmodel)[[1]]),
        c('name', 'TSPs', 'score', 'labels'))
    # added correct metadata during training
    expect_equal(names(metadata(trainedPCOSPmodel)$modelParams),
        c('randomSeed', 'RNGkind', 'minDaysSurvived', 'numModels',
            'minAccuracy'))
})

test_that('predictClasses method for PCOSP works correctly', {
    expect_warning(predictClasses(sampleCohortList[[1]], model=trainedPCOSPmodel),
        '.*Dropped samples with NA survival data.*')
    expect_true(all(mcols(PCOSPpredCohortList)$hasPredictions))
    expect_true('PCOSP_prob_good' %in%
        colnames(colData(PCOSPpredCohortList[[1]])))
    expect_true(is.matrix(metadata(PCOSPpredCohortList[[1]])$PCOSPpredictions))
})

test_that('validateModel method for PCOSP  works correctly', {
    expect_message(validateModel(trainedPCOSPmodel,
        valData=PCOSPpredCohortList[[1]]), 'Setting levels: .*')
    expect_true(metadata(validatedPCOSPModel)$isValidated)
    expect_true(is.data.frame(validationStats(validatedPCOSPModel)))
    expect_true(all(
        c('statistic', 'estimate', 'se', 'lower', 'upper', 'p_value',
            'n', 'cohort', 'mDataType') %in%
                colnames(validationStats(validatedPCOSPModel))))
    expect_equal(validationData(validatedPCOSPModel), PCOSPpredCohortList)
})
