library(testthat)
library(PDATK)

data(sampleICGCmicro)
data(sampleCohortList)

test_that('ClinicalModel constructor errors as expected', {
    expect_error(ClinicalModel(sampleICGCmicro, formula=25),
        '.*The formula for a clinical model must either be a formula.*')
})

test_that('ClinicalModel constructor works as expected', {
    expect_s4_class({
        clinModel <- ClinicalModel(sampleICGCmicro,
            formula='prognosis ~ age + sex'); clinModel
    } , 'ClinicalModel')
    expect_true(validObject(clinModel))
})