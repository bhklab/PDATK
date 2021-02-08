library(testthat)
library(PDATK)

data(sampleICGCmicro)

SE <- as(sampleICGCmicro, 'SummarizedExperiment')

test_that('SummarizedExperiment can be coerce to SurvivalExperiment', {
    survExp <- as(SE, 'SurvivalExperiment')
    expect_s4_class(survExp, 'SurvivalExperiment')
    expect_true(validObject(survExp))
})