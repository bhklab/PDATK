library(testthat)
library(PDATK)

data(sampleICGCmicro)
data(sampleCohortList)

suppressWarnings({
    PCOSPmodel <- PCOSP(sampleICGCmicro, randomSeed=1987)
    trainedPCOSPmodel <- trainModel(PCOSPmodel, numModels=10)
    PCOSPpredCohortList <- predictClasses(sampleCohortList[seq_len(2)],
        model=trainedPCOSPmodel)
    validatedPCOSPModel <- validateModel(trainedPCOSPmodel,
        valData=PCOSPpredCohortList)
})

test_that('ModelComparison constructor works with two SurvivalModel
    sub-classes',
{
    expect_s4_class(compareModels(validatedPCOSPModel, validatedPCOSPModel),
        'ModelComparison')
})

test_that('ModelComparison constructor works with one ModelComparison
    and one SurvivalModel sub-class',
{
    modComp <- compareModels(validatedPCOSPModel, validatedPCOSPModel)
    expect_s4_class(compareModels(modComp, validatedPCOSPModel),
        'ModelComparison')
})

## TODO:: Implement this method
# test_that('ModelComparison constructor works with two ModelComparison
#     objects',
# {
#     modComp <- compareModels(validatedPCOSPModel, validatedPCOSPModel)
#     expect_s4_class(compareModels(modComp, modComp), 'ModelComparsion')
# })