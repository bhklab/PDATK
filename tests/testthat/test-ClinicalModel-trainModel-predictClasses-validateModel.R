library(testthat)
library(PDATK)
library(BiocParallel)


data(sampleICGCmicro)
data(sampleCohortList)

suppressWarnings({
    if (Sys.info()['sysname'] == 'windows') {
        BiocParallel::register(BiocParallel::SerialParam())
    }
    clinicalModel <- ClinicalModel(sampleICGCmicro,
        formula='prognosis ~ sex + age + T + N + M + grade', randomSeed=1987)
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

test_that('predictClasses for Clinical model warnings/errors as expected',{
    expect_warning(
        predictClasses(sampleCohortList[[1]], model=trainedClinicalModel),
            '.*Rows .* have levels that are not in the model, skipping.*')
    expect_error(predictClasses(sampleCohortList$UNC, model=trainedClinicalModel),
        '.*The columns .* are missing from the colData slot of the .*')
    models(trainedClinicalModel) <- c(models(trainedClinicalModel),
        SimpleList(anotherModel=NA))
    invalidClinicalModel <- trainedClinicalModel[, # exclude patients with NAs
        rownames(na.omit(colData(trainedClinicalModel)))]
    expect_warning(
        predictClasses(invalidClinicalModel, model=trainedClinicalModel),
        '.*There is more than one model in your ClinicalModel.*')
})

test_that('predictClasses method for ClinicalModel works correctly', {
    expect_true(all(mcols(clinicalPredCohortList)$hasPredictions))
    expect_true('clinical_prob_good' %in%
        colnames(colData(clinicalPredCohortList[[1]])))
})

## TODO:: implement test
# test_that('validateModel method for ClinicalModel warnings/errors correctly',
# {
#
# })

test_that('validateModel method for ClinicalModel works correctly', {
    expect_message(
        validateModel(trainedClinicalModel, valData=clinicalPredCohortList[[1]]),
        'Setting levels: .*|.*sample(s) with missing data.*')
    expect_true(metadata(validatedClinicalModel)$isValidated)
    expect_true(is.data.frame(validationStats(validatedClinicalModel)))
    expect_equal(colnames(validationStats(validatedClinicalModel)),
        c('statistic', 'estimate', 'se', 'lower', 'upper', 'p_value',
            'n', 'isSummary', 'cohort', 'mDataType'))
    expect_equal(validationData(validatedClinicalModel), clinicalPredCohortList)
})
