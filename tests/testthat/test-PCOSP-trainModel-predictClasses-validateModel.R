library(testthat)
library(PDATK)

data(sampleICGCmicro)
data(sampleCohortList)

PCOSPmodel <- PCOSP(sampleICGCmicro)
trainedPCOSPmodel <- trainModel(PCOSPmodel, numModels=10)
PCOSPpredCohortList <- predictClasses(sampleCohortList,
    model=trainedPCOSPmodel)
validatedPCOSPModel <- validateModel(trainedPCOSPmodel,
    valData=PCOSPpredCohortList)

test_that('PCOSP method for trainModel works correctly', {

})

test_that('PCOSP method for predictClasses works correctly', {

})

test_that('PCOSP method for validateModel works correctly')