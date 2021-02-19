#' Plot ROC curves for an `S4` object
#'
#' @param object An `S4` object with a defined plotROC method.
#' @param ... Allow new parameters to be added to this generic.
#'
#' @return A `ggplot` object containing the ROC curves.
#'
#' @md
#' @export
setGeneric('plotROC', function(object, ...)
    standardGeneric('plotROC'))
#'
#' Plot ROC curves for a `PCOSP` model object.
#'
#' @param object A `PCOSP` model which has been validated with with the
#'  `validateModel` method.
#' @param alpha A `float` specifying the significance level for the plot. Non-
#'   signficiant cohorts will have dotted lines.
#' @param ... Catch unnamed parameters. Not used.
#' @param xlabel A `character` vector specifying the x label.
#' @param ylabel A `character` vector specifying the y label.
#' @param title A `character` vector speciyfing the plot tile.
#'
#' @return A `ggplot` object containing the ROC curves.
#'
#' @examples
#' data(sampleTrainedPCOSPmodel)
#' data(sampleCohortList)
#'
#' # Make predictions
#' PCOSPpredCohortList <- predictClasses(sampleCohortList[1:4],
#'   model=sampleTrainedPCOSPmodel)
#'
#' # Validate the model
#' validatedPCOSPmodel <- validateModel(sampleTrainedPCOSPmodel,
#'   PCOSPpredCohortList)
#'
#' # Plot ROC curves
#' AUROCplot <- plotROC(validatedPCOSPmodel)
#'
#' @md
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder
#' @importFrom pROC roc coords
#' @importFrom ggplot2 ggplot geom_segment geom_step scale_x_reverse
#'   theme_classic scale_linetype_manual guides theme xlab ylab aes
#'   ggtitle labs
#' @export
setMethod('plotROC', signature(object='PCOSP'),
    function(object, alpha=0.05, ..., xlabel, ylabel, title)
{
    funContext <- .context(1)

    # -- check for the correct columns
    ##TODO:: use mcols to specify if validation data has predictions
    valData <- validationData(object)
    if (length(valData) < 1)
        stop(.errorMsg(funContext, 'There is no validation data in this',
            'PCOSP object, please ensure you run validateModel before trying',
            ' to plot ROC curves!'))

    survivalDfList <- lapply(valData, colData)
    .hasPCOSPcol <- function(colData) 'PCOSP_prob_good' %in% colnames(colData)
    if (!all(unlist(lapply(survivalDfList, .hasPCOSPcol)))) {
        stop(.errorMsg(funContext, 'One or more of the SurvivalExperiment ',
            'objects in the validationData slot of the PCOSP object are missing ',
            'the PCOSP_prob_good column. Please rerun predictClasses on your ',
            'validation data, then rerun validateModel with the updated ',
            'validation data for your model.'))
    }

    # -- calcualte the ROC sensitivities and specificities
    .calcROC <- function(df) {
        rocRes <- with(df,
            roc(ifelse(prognosis =='good', 1, 0),
                PCOSP_prob_good))
        DT <- as.data.table(coords(rocRes, 'all', transpose=FALSE))
        DT[rev(seq_len(.N))]
    }

    rocDtList <- lapply(survivalDfList, .calcROC)
    for (i in seq_along(rocDtList)) {
        rocDtList[[i]][, cohort := names(rocDtList)[i]]
    }

    rocDT <- rbindlist(rocDtList)
    aucStatsDT <- validationStats(object)[statistic=='AUC'][cohort %in% names(valData)]

    rocDT <- merge.data.table(rocDT, aucStatsDT[, .(cohort, estimate, p_value)],
        by='cohort')

    rocDT[, cohort := paste0(cohort, ': ', estimate, ' (', scientific(p_value, 2), ')')]

    plot <- ggplot(rocDT, aes(x=specificity, y=sensitivity,
                    color=cohort,
                    linetype=p_value < alpha)) +
                geom_segment(aes(x=max(specificity), xend=min(specificity),
                    y=min(sensitivity), yend=max(sensitivity)), colour='grey',
                    size=0.1) +
                geom_step(size=0.75) +
                scale_x_reverse() +
                theme_classic() +
                scale_linetype_manual(breaks=c(TRUE, FALSE), values=c(1, 2)) +
                labs(linetype=paste0('Significant (p<', alpha, ')'),
                    colour='Cohort: AUC (P-value)') +
                guides(linetype=FALSE) +
                theme(legend.justification=c(1,0), legend.position=c(1,0.005)) +
                xlab('Specificity') +
                ylab('Sensitivity')

    if (!missing(xlabel))
        plot + xlab(xlabel)
    if (!missing(ylabel))
        plot + ylab(ylabel)
    if (!missing(title))
        plot + ggtitle(title)

    return(plot)
})
