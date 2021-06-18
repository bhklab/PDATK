library(testthat)
library(PDATK)

data(sampleICGCmicro)

test_that('SurivalExperiment constructor errors as expected', {
    expect_error(SurvivalExperiment(data=list()), 'unused argument')
})

test_that('SurivalExperiment constructor errors for missing days_surived or
    event_occurred column arguments',
    {
    expect_error(SurvivalExperiment(sampleICGCmicro,
        event_occurred='wrong_name'), '.*The columns .* are not present in the .*!')
})

test_that('SurvivalExperiment constructor successfully coerces event_occurred
    and survival_time when they are the wrong types',
{
    colData(sampleICGCmicro)$survival_time <-
        as.numeric(colData(sampleICGCmicro)$survival_time)
    colData(sampleICGCmicro)$event_occurred <-
        as.logical(colData(sampleICGCmicro)$event_occurred)
    expect_s4_class(SurvivalExperiment(sampleICGCmicro),
        'SurvivalExperiment')
    colData(sampleICGCmicro)$event_occurred <-
        ifelse(colData(sampleICGCmicro)$event_occurred, 'deceased', 'not_deceased')
    expect_s4_class(SurvivalExperiment(sampleICGCmicro),
        'SurvivalExperiment')
})

test_that('SurvivalExperiment constructor errors if survival_time or
    event_occurred columns are wrong types and coercion is not possible',
    {
    invalidExp1 <- invalidExp2 <- sampleICGCmicro
    colData(invalidExp1)$event_occurred <-
        as.factor(colData(invalidExp1)$event_occurred)
    colData(invalidExp2)$survival_time <-
        as.factor(colData(invalidExp2)$survival_time)
    expect_error(SurvivalExperiment(invalidExp1),
        '.*deceased is not in.*')
})

test_that('SurvivalExperiment constructor errors if survival_time or
    event_occurred columns are character but not coercibel to integers',
    {
    invalidExp1 <- invalidExp2 <- sampleICGCmicro
    # Character not coercible to integer
    colData(invalidExp1)$event_occurred <- 'A'
    colData(invalidExp2)$survival_time <- 'A'
    expect_error(SurvivalExperiment(invalidExp1),
        '.*The string deceased is not in the event_occurred column.*')
    expect_error(SurvivalExperiment(invalidExp2),
        '.*Tried to coerce survival_time from character to integer.*')
})

test_that('SurvivalExperiment constructor can be made from raw data', {
    survExp <- SurvivalExperiment(assays=assays(sampleICGCmicro),
        colData=colData(sampleICGCmicro), rowData=rowData(sampleICGCmicro),
        metadata=metadata(sampleICGCmicro))
    expect_s4_class(survExp, 'SurvivalExperiment')
    expect_true(validObject(survExp, complete=TRUE))
})

test_that('SurvivalExperimemt constructor can be made from
    SummarizedExperiment',
{
    survExp <- SurvivalExperiment(sampleICGCmicro)
    expect_s4_class(survExp, 'SurvivalExperiment')
    expect_true(validObject(survExp, complete=TRUE))
})

test_that('SurvivalExperiment constructor works to create empty
    SurvivalExperiment',
{
    expect_s4_class(SurvivalExperiment(), 'SurvivalExperiment')
})
