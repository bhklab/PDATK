library(testthat)
library(PDATK)

data(sampleICGCmicro)

test_that('SurvivalModel constructor errors with wrong trainCohorts type', {
    expect_error(SurvivalModel(trainCohorts=list()),
        '.*The trainCohorts argument is not a CohortList or SurvivalExperiment.*')
})

test_that('SurvivaModel constructor tells user about type coercion', {
    expect_message(
        SurvivalModel(
            CohortList(list(sampleICGCmicro, sampleICGCmicro),
                mDataTypes=c('rna_micro', 'rna_micro'))),
        '.* Merging trainCohorts `CohortList` into a `SurvivalExperiment` with shared genes and samples.*')
})

test_that('SurviavlModel constructor works with SurvivalExperiment trainData', {
    survModel <- SurvivalModel(sampleICGCmicro)
    expect_s4_class(survModel, 'SurvivalModel')
    expect_true(validObject(survModel))
})

test_that('SurviavlModel constructor works with CohortList, trainData', {
    suppressMessages({
        survModel <- SurvivalModel(
            CohortList(list(sampleICGCmicro, sampleICGCmicro),
                mDataTypes=c('rna_micro', 'rna_micro')))
    })
    expect_s4_class(survModel, 'SurvivalModel')
    expect_true(validObject(survModel))
})