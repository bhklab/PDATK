library(testthat)
library(PDATK)

data(sampleICGCmicro)
data(sampleCohortList)

test_that('RGAModel constructor works with SurvivalExperiment
    trainCohorts',
{
    expect_s4_class({
        RGAmodel <- RGAModel(sampleICGCmicro); RGAmodel
    },'RGAModel')
    expect_true(validObject(RGAmodel))
})

test_that('RGAModel constuctor works with CohortList
    trainCohorts',
{
    suppressMessages({
        expect_s4_class({
            RGAmodel <- RGAModel(CohortList(list(cohort1=sampleICGCmicro,
                cohort2=sampleICGCmicro))); RGAmodel
        }, 'RGAModel')
        expect_true(validObject(RGAmodel))
    })
})