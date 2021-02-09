library(testthat)
library(PDATK)

data(sampleICGCmicro)
data(sampleCohortList)

test_that('RLSModel constructor works with SurvivalExperiment
    trainCohorts',
{
    expect_s4_class({
        RLSmodel <- RLSModel(sampleICGCmicro); RLSmodel
    },'RLSModel')
    expect_true(validObject(RLSmodel))
})

test_that('RLSModel constuctor works with CohortList
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