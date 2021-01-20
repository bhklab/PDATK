#' Plot ROC curves for an `S4` object
#'
#' @param object An `S4` object with a defined plotROC method.
#' @param ... Allow new parameters to be added to this generic.
#'
#' @return A `ggplot2` object containing the ROC curves.
#'
#' @export
setGeneric('plotROC', function(object, ...)
    standardGeneric('plotROC'))
#'
#'
#'
#'
#'
setMethod('plotROC', signature(object='PCOSP'),
    function(object, )
{
    # -- check for the correct columns
    ##TODO:: use mcols to specify if validation data has predictions
    valData <- validationData(object)
    if (length(valData) < 1)
        stop(.errorMsg(.context(), 'There is no validation data in this',
            'PCOSP object, please ensure you run validateModel before trying',
            ' to plot ROC curves!'))

    survivalDfList <- lapply(valData, colData)
    .hasPCOSPcol <- function(colData) 'PCOSP_prob_good' %in% colnames(colData)
    if (!all(unlist(lapply(survivalDfList, .hasPCOSPcol)))) {
        stop(.errorMsg(.context(), 'One or more of the SurvivalExperiment ',
            'objects in the validationSlot of the PCOSP object are missing ',
            'PCOSP_prob_good column. Please rerun predictClasses on your ',
            'validation data, then rerun validateModel with the updated ',
            'validation data for your model.'))
    }

    # -- calcualte the ROC sensitivities and specificities
    .calcROC <- function(df) {
        rocRes <- roc(df$is_deceased, 1 - df$PCOSP_prob_good)
        data.table(sensitivity=rocRes$sensitivities,
            specificity=rocRes$specificities)
    }

    rocDtList <- lapply(survivalDfList, .calcROC)
    for (i in seq_along(rocDtList)) {
        rocDtList[[i]][, cohort := names(rocDtList)[i]]
    }

    rocDT <- rbindlist(rocDtList)
})