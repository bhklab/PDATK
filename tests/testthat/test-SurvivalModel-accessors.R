library(testthat)
library(PDATK)

context("Testing SurvivalModel Class Accessor Methods")

test_that('SurivalModel constructor errors as expected', {
    data(sampleICGCmicro)
    sampleData
    expect_error(SurvivalExperiment(data=list()), 'unused argument')
    expect_error(SurvivalExperiment(sumExp=list()),
        '.*[PDATK::SurvivalExperiment] sumExp is not a `SummarizedExperiment`!.*')
    expect_error(SurvivalExperiment(sumExp=sampleICGCmicro, is_deceased='wrong_name'),
        '.*[PDATK::SurvivalExperiment] The columns .* are not present in the .*!')
})

test_that('SurvivalModel constructor can be made from raw data', {
    expect_is(SurvivalExperiment(assays=assays(sampleICGCmicro)),
        'SurvivalExperiment')
})

test_that('SurvivalModel constructor can be made from SummarizedExperiment', {

})

test_that(")