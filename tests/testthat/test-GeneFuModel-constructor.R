library(testthat)
library(PDATK)

data(sampleICGCmicro)
data(sampleCohortList)

test_that('GeneFuModel constructor works correctly', {
    expect_s4_class({ gfModel <- GeneFuModel(); gfModel },
        'GeneFuModel')
    expect_true(validObject(gfModel))
})