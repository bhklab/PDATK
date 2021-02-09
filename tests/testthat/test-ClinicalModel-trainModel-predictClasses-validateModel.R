library(testthat)
library(PDATK)

data(sampleICGCmicro)
data(sampleCohortList)

suppressWarnings({
    clinicalModel <- ClinicalModel(sampleICGCmicro,
        formula='prognosis ~ sex + age + T + N + M + grade')
    trainedClinicalModel <- trainModel(clinicalModel)
    clinicalPredCohortList <- predictClasses(sampleCohortList[seq_len(2)],
        model=trainedClinicalModel)
    validatedClinicalModel <- validateModel(trainedClinicalModel,
        valData=clinicalPredCohortList)
})

test_that('trainModel method for ClinicalModel works correctly', {
    expect_s4_class(trainModel(clinicalModel), 'ClinicalModel')
    expect_true(length(models(trainedClinicalModel)) == 1)
    # the model list has the correct items
    expect_s3_class(models(trainedClinicalModel)[[1]], 'glm')
    # added correct metadata during training
    expect_equal(names(metadata(trainedClinicalModel)$modelParams),
        c('randomSeed', 'RNGkind', 'minDaysSurvived', 'formula'))
})

test_that('predictClasses method for ClinicalModel works correctly', {
    expect_warning(predictClasses(sampleCohortList[[1]], model=trainedClinicalModel),
        '.*Rows .* have levels that are not in the model, skipping.*')
    expect_true(all(mcols(PCOSPpredCohortList)$hasPredictions))
    expect_true('clinical_prob_good' %in%
        colnames(colData(PCOSPpredCohortList[[1]])))
})

test_that('validateModel method for ClinicalModel works correctly', {
    expect_message(validateModel(trainedPCOSPmodel,
        valData=PCOSPpredCohortList[[1]]), 'Setting levels: .*')
    expect_true(metadata(validatedPCOSPModel)$isValidated)
    expect_true(is.data.frame(validationStats(validatedPCOSPModel)))
    expect_equal(colnames(validationStats(validatedPCOSPModel)),
        c('statistic', 'estimate', 'se', 'lower', 'upper', 'p_value',
            'n', 'isSummary', 'cohort', 'mDataType'))
    expect_equal(validationData(validatedPCOSPModel), PCOSPpredCohortList)
})