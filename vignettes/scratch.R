#trainedPCOSPmodel <- trainModel(PCOSPmodel, numModels=10, minAccurary=0.6)

validatedPCOSPmodel <- validateModel(trainedPCOSPmodel,
    valData=PCOSPpredValCohorts)