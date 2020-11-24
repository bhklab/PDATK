#' Perform Validation on an `S4` Object Respresenting a Trained Model
#'
#' @param model An `S4` object
#'
#' @md
#' @export
setGeneric('validateModel', function(model, validationData, ...)
    standardGeneric('validateModel'))
#'
#' Evaluate the Performance of a List of Trained KTSP Models from a PCOSP
#'   Model
#'
#' @param model A `PCOSP` model which has been trained using `trainModel`.
#' @param validationData A `CohortList` containing one or more
#'   `SurvivalExperiment`s. The first assay in each `SurvivalExperiment` will
#'   be classified using all top scoring KTSP models in `models(model)`.
#' @param ... Fallthrough arguments to `BiocParallel::bplapply`, use this to
#'   configure the parallelization settings for this function. For example
#'   to specify BPARAM.
#'
#' @seealso BiocParallel::bplapply switchBox::SWAP.KTSP.Classify
#'
#' @return
#'
#' @importFrom BiocParallel bplapply
#' @importFrom switchBox SWAP.KTSP.Classify
#' @md
#' @export
setMethod('validateModel', signature(model='PCOSP',
    validationData='CohortList'), function(model, validationData, ...)
{

    # determine if the validation data already has predictions and if
    #   if the predictions were made using the same model
    if ('hasPredictions' %in% colnames(mcols(validationData))) {
        if (all(mcols(validationData)$hasPredictions)) {
            if (all.equal(model, metadata(validationData)$predictionModel)) {
                predCohortList <- validationData
            } else {
                warning(.warnMsg(.context(), 'The validationData argument ',
                    'has predictions, but the prediction model does not match',
                    'the model argument. Recalculation classes...'))
            }
        }
    } else {
        predCohortList <- predictClasses(validationData, model=model)
    }

    validationCohortList <- endoapply(predCohortList, validateModel, model=model)

    .getValidationStats <- function(x) metadata(x)$validationStats

    validationStats <- lapply(validationCohortList, .getValidationStats)

    for (i in seq_along(validationStats)) {
        validationStats[[i]]$cohort <- names(validationCohortList)[i]
        validationStats[[i]]$mDataType <- mcols(validationCohortList)[i]
    }

    validationDT <- rbindlist(validationStats)

})
#'
#' @param model A `PCOSP` model which has been trained using `trainModel`.
#' @param validationData A `SurvivalExperiment` to validate the model with.
#'
#' @return A `SurvivalExperiment` with the predicted classes and validation
#'   stats in the metadata slot.
#'
#' @seealso BiocParallel::bplapply switchBox::SWAP.KTSP.Classify
#'
#' @md
#' @importFrom survcomp D.index concordance.index combine.est
#' @export
setMethod('validateModel', signature(model='PCOSP',
    validationData='SurvivalExperiment'), function(model, validationData, ...)
{
    # determine if we need to rerun the classification model
    if (identical(metadata(model)$modelParams, metadata(validationData)$PCOSPparams))
    {
        survivalDF <- colData(validationData)[, c('sample_name', 'days_survived',
            'is_deceased', 'PCOSP_prob_good')]
        predSurvExp <- validationData
    } else {
        predSurvExp <- predictClasses(model, validationData)
        survivalDF <- colData(predSurvExp)[, c('sample_name', 'days_survived',
            'is_deceased', 'PCOSP_prob_good')]
    }

    # calculate the validation statistcs
    validationStats <- with(survivalDF,
        list(
            dIndex=D.index(x=1 - PCOSP_prob_good, surv.time=days_survived,
                surv.event=is_deceased, na.rm=TRUE, alpha=0.5,
                method.test='logrank'),
            cIndex=concordance.index(x=1 - PCOSP_prob_good, surv.time=days_survived,
                surv.event=is_deceased, method='noether', na.rm=TRUE)
        )
    )

    # assemble into a data.frame
    valStatsDF <- data.frame(
        statistic=c('D_index', 'concordance_index'),
        value=c(validationStats$dIndex$d.index, validationStats$cIndex$c.index),
        se=vapply(validationStats, `[[`, i='se', FUN.VALUE=numeric(1)),
        lower=vapply(validationStats, `[[`, i='lower', FUN.VALUE=numeric(1)),
        upper=vapply(validationStats, `[[`, i='upper', FUN.VALUE=numeric(1)),
        p_value=vapply(validationStats, `[[`, i='p.value', FUN.VALUE=numeric(1)),
        n=vapply(validationStats, `[[`, i='n', FUN.VALUE=numeric(1))
    )

    metadata(predSurvExp)$validationStats <- valStatsDF
    return(predSurvExp)
})