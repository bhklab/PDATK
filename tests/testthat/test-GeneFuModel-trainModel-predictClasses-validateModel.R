library(testthat)
library(PDATK)
library(S4Vectors)

data(existingClassifierData, verbose=TRUE)
data(sampleCohortList)

if (Sys.info()['sysname'] == 'Windows') {
    BiocParallel::register(BiocParallel::SerialParam())
}

set.seed(1987)
chenGeneFuModel <- GeneFuModel(randomSeed=1987)
models(chenGeneFuModel) <- SimpleList(list(chen=chen))
chenGeneFuPredCohortList <- predictClasses(sampleCohortList[1],
    model=chenGeneFuModel)
validatedChenModel <- validateModel(chenGeneFuModel,
    valData=chenGeneFuPredCohortList)

test_that('trainModel method for GeneFuModel errors as expected', {
    expect_error(trainModel(chenGeneFuModel),
        '.*Unfortunately we have not implemented model training for.*')
})

test_that('predictClasses method for GeneFuModel class errors as expected',
{
    models(chenGeneFuModel) <- SimpleList(chen=chen, chen_again=chen)
    expect_warning(predictClasses(sampleCohortList[1], model=chenGeneFuModel),
        '.*There is more than one model in the models slot of the.*')
})

test_that('predictClasses method for GeneFuModel class works correctly', {
    expect_s4_class(predictClasses(sampleCohortList[1], model=chenGeneFuModel),
        'CohortList')
    expect_true(all(mcols(chenGeneFuPredCohortList)$hasPrediction))
    expect_true('genefu_score' %in% colnames(colData(chenGeneFuPredCohortList[[1]])))
})

test_that('validateModel method for GeneFuModel class errors as expected',
{
    predSurvExp <- chenGeneFuPredCohortList[[1]]
    colData(predSurvExp)$prognosis <- NULL
    expect_warning(validateModel(chenGeneFuModel, predSurvExp),
        '.*The prognosis column is missing from the validation.*')
})

test_that('validateModel method for GeneFuModel class works correctly', {
    expect_s3_class(validationStats(validatedChenModel), 'data.frame')
    expect_true(
        all(c("statistic", "estimate", "se", "lower", "upper", "p_value",
            "n", "isSummary", "cohort", "mDataType") %in%
                colnames(validationStats(validatedChenModel))))
    expect_true(metadata(validatedChenModel)$isValidated)
})
