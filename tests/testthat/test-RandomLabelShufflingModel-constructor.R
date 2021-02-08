library(testthat)
library(PDATK)

data(sampleICGCmicro)
data(sampleCohortList)

test_that('RandomLabelShufflingModel constructor works with SurvivalExperiment
    trainCohorts',
{
    expect_s4_class({
        RLSmodel <- RLSModel(sampleICGCmicro); RLSmodel
    },'RLSModel')
    expect_true(validObject(RLSmodel))
})

test_that('RandomLabelShufflingModel constuctor works with CohortList
    trainCohorts',
{
    suppressMessages({
        expect_s4_class({
            RLSmodel <- RLSModel(CohortList(list(cohort1=sampleICGCmicro,
                cohort2=sampleICGCmicro))); RLSmodel
        }, 'RLSModel')
        expect_true(validObject(RLSmodel))
    })
})