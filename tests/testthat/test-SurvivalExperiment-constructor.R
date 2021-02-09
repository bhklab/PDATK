library(testthat)
library(PDATK)

data(sampleICGCmicro)

test_that('SurivalExperiment constructor errors as expected', {
    expect_error(SurvivalExperiment(data=list()), 'unused argument')
    expect_error(SurvivalExperiment(sumExp=list()),
        '.*sumExp is not a `SummarizedExperiment`!.*')
})

test_that('SurivalExperiment constructor errors for missing days_surived or
    is_deceased column arguments',
    {
    expect_error(SurvivalExperiment(sumExp=sampleICGCmicro,
        is_deceased='wrong_name'), '.*The columns .* are not present in the .*!')
})

test_that('SurvivalExperiment constructor successfully coerces is_deceased
    and days_survived when they are the wrong types',
{
    colData(sampleICGCmicro)$days_survived <-
        as.float(colData(sampleICGCmicro)$days_survived)
    colData(sampleICGCmicro)$is_deceased <-
        as.logical(colData(sampleICGCmicro)$is_deceased)
    expect_s4_class(SurvialExpeiment(sampleICGCmicro))
})

test_that('SurvivalExperiment constructor errors if days_survived or
    is_deceased columns are wrong types and coercion is not possible',
    {
    invalidExp1 <- invalidExp2 <- sampleICGCmicro
    colData(invalidExp1)$is_deceased <-
        as.factor(colData(invalidExp1)$is_deceased)
    colData(invalidExp2)$days_survived <-
        as.factor(colData(invalidExp2)$days_survived)
    expect_error(SurvivalExperiment(sumExp=invalidExp1),
        '.*The is_deceased column is not logical or integer,.*')
    expect_error(SurvivalExperiment(sumExp=invalidExp2),
        '.*The days_survived column is not logical or integer,.*')
})

test_that('SurvivalExperiment constructor errors if days_survived or
    is_deceased columns are character but not coercibel to integers',
    {
    invalidExp1 <- invalidExp2 <- sampleICGCmicro
    # Character not coercible to integer
    colData(invalidExp1)$is_deceased <- 'A'
    colData(invalidExp2)$days_survived <- 'A'
    expect_error(SurvivalExperiment(sumExp=invalidExp1),
        '.*The string deceased is not in the is_deceased column.*')
    expect_error(SurvivalExperiment(sumExp=invalidExp2),
        '.*Tried to coerce days_survived from character to integer.*')
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
    survExp <- SurvivalExperiment(sumExp=sampleICGCmicro)
    expect_s4_class(survExp, 'SurvivalExperiment')
    expect_true(validObject(survExp, complete=TRUE))
})

test_that('SurvivalExperiment constructor works to create empty
    SurvivalExperiment',
{
    expect_s4_class(SurvivalExperiment(), 'SurvivalExperiment')
})
