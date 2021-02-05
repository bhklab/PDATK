#' ModelComparison Class Definition
#'
#' @md
#' @importClassesFrom S4Vectors DataFrame
#' @export
.ModelComparison <- setClass('ModelComparison', contains='DataFrame')


#' ModelComparison Constructor
#'
#' @param model1 An object with a `validationStats` method which returns a
#'   `data.table`. Probably this object should inherit from `SurvivalModel`.
#' @param model2 An object with a `validationStats` method which returns a
#'   `data.table`. Probably this object should inherit from `SurvivalModel`.
#' @param ... Not used.
#'
#' @return A `ModelComparison` object, which is a wrapper for `DataFrame`
#'   which is used for method dispatch.
#'
#' @examples
#' data(samplePCOSPmodel)
#' data(sampleClinicalModel)
#' data(samplePCSIsurvExp)
#'
#' # Train the models
#' trainedPCOSPmodel <-trainModel(samplePCOSPmodel, numModels=10, minAccuracy=0.6)
#' trainedClinicalModel <- trainModel(sampleClinicalModel)
#'
#' # Predict risk/risk-class
#' PCOSPpredPCSI <- predictClasses(samplePCSIsurvExp, model=trainedPCOSPmodel)
#' ClinicalPredPCSI <- predictClasses(samplePCSIsurvExp,
#'   model=trainedClinicalModel)
#'
#' # Validate the models
#' validatedPCOSPmodel <- validateModel(trainedPCSOPmodel,
#'   valData=PCOSPpredPCSI)
#' validatedClinicalModel <- validateModel(trainedClinicalModel,
#'   valData=ClinicalPredPCSI)
#'
#' # Compare the models
#' modelComp <- ModelComparison(validatedPCOSPmodel, validatedClinicalModel)
#' head(modelComp)
#'
#' @md
#' @importFrom methods is
#' @importFrom S4Vectors DataFrame
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder
#' @export
ModelComparison <- function(model1, model2, ...) {


    ## TODO:: Is it better to define a validationStats method for a
    ##   ModelComparsion? Then can't do class for model column.
    if (is(model1, 'ModelComparison')) {
        model1StatsDT <- as.data.table(model1)
    } else {
        model1StatsDT <- validationStats(model1)
        model1StatsDT[, model := class(model1)]
    }

    if (is(model2, 'ModelComparison')) {

    } else {
        model2StatsDT <- validationStats(model2)
        model2StatsDT[, model := class(model2)]
    }

    sharedCohorts <- intersect(model1StatsDT$cohort, model2StatsDT$cohort)
    modelComparisonDT <- rbind(model1StatsDT, model2StatsDT)
    modelComparisonDT <- modelComparisonDT[cohort %in% sharedCohorts, ]

    .ModelComparison(DataFrame(modelComparisonDT))
}
