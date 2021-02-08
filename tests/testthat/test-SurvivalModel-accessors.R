library(testthat)
library(PDATK)

data(sampleICGCmicro)
survModel <- SurvivalModel(trainCohorts=sampleICGCmicro)

test_that('SurvivalModel models accessors work correctly', {
    expect_s4_class(models(survModel), 'SimpleList')
    expect_silent(models(survModel) <- models(survModel))
})

test_that('SurvivalModel validationStats accessors work correctly', {
    expect_s3_class(validationStats(survModel), 'data.frame')
    expect_silent(validationStats(survModel) <- validationStats(survModel))
})

test_that('SurvivalModel validationData accessors work correctly', {
    expect_s4_class(validationData(survModel), 'CohortList')
    expect_silent(validationData(survModel) <- validationData(survModel))
})