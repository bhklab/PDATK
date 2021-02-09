#' Render A Density Plot of Model Performance for an `S4` Object
#'
#' @param object An `S4` object representing a statistical or ML model.
#' @param refModel An `S4` object representing a statistical or ML model to
#'   compare `object` against.
#' @param ... Allow additional parameters to be defined for this generic.
#'
#' @return A `ggplot` object containing the plot.
#'
#' @md
#' @export
setGeneric('densityPlotModelComparison',
    function(object, refModel, ...)
        standardGeneric('densityPlotModelComparison'))
#'
#' Render a Density Plot Comparing Model Performance Between Two `PCOSP`,
#'   `RLSModel` or `RGAModel` object.
#'
#' @param object A `PCOSP`, `RLSModel` or `RGAModel` object.
#' @param refModel A `PCOSP`, `RLSModel` or `RGAModel` object to compare
#'   performance against.
#' @param ... Catch unnamed parameters. Not used.
#' @param title Optional `character` vector with plot title.
#' @param xlab Optional `character` vector with x-axis label.
#' @param ylab Optional `character` vector with y-axis label.
#' @param mDataTypeLabels Optional `character` vector whos names are one
#'   or more existing mDataTypes in `object` and `refModel` and whos values
#'   are the desired mDataType labels in the plot facets.
#'
#' @return A `ggplot` object with a density plot of model AUCs for `object`
#'   and a vertical line for the average AUC of `refModel`, faceted by
#'   mDataType.
#'
#' @md
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse melt.data.table transpose setcolorder
#' @importFrom pROC roc
#' @importFrom ggplot2 ggplot geom_density geom_vline labs facet_wrap ggtitle
#'   xlab ylab labeller element_line element_text element_blank
#' @export
setMethod('densityPlotModelComparison',
    signature(object='PCOSP_or_RLS_or_RGA', refModel='PCOSP_or_RLS_or_RGA'),
    function(object, refModel, ..., title, xlab, ylab, mDataTypeLabels)
{
    # Extract the per model prognosis predictions
    valData <- validationData(object)
    metaData <- lapply(valData, metadata)
    predictions <- lapply(metaData, `[[`,
        paste0(class(object)[1], 'predictions'))
    predictions <- lapply(predictions, as.data.table, keep.rownames='model')
    predictions <- lapply(predictions, transpose, keep.names='sample',
        make.names='model')

    # Extract the true model prognosis values
    columnData <- lapply(valData, colData)
    groups <- lapply(columnData, `[[`, 'prognosis')

    # Turn into data.tables with appropriate labels
    for (i in seq_along(predictions)) {
        predictions[[i]][, cohort := names(predictions)[i]]
        groups[[i]] <- data.table(prognosis=groups[[i]],
            sample=rownames(columnData[[i]]), cohort=names(groups)[i],
            mDataType=mcols(valData)$mDataType[i])
    }

    # join the tables
    predictionDT <- rbindlist(predictions)
    groupDT <- rbindlist(groups)
    modelDT <- merge.data.table(predictionDT, groupDT, by=c('cohort', 'sample'))

    # assess predictions for each model in each cohort
    .predictModelAUC <- function(col, groups) {
        suppressMessages({
            roc(response=fifelse(groups == 'good', 1, 0),
                predictor=fifelse(col == 'good', 1, 0))$auc
        })
    }

    # put `object` statistics in a data.table
    overall <- modelDT[, lapply(.SD, .predictModelAUC, groups=prognosis),
        .SDcols=patterns('rank.*')]
    overall[, mDataType := 'combined']

    mDataType <- modelDT[, lapply(.SD, .predictModelAUC, groups=prognosis),
        .SDcols=patterns('rank.*'), by=mDataType]

    modelPerformanceDT <- rbind(overall, mDataType)
    modelPerfDT <- melt.data.table(modelPerformanceDT, value.name='AUC',
        value.vars=patterns('rank*'), variable.name='model', id.vars='mDataType')

    # add `refModel` AUCs by mDataType to the model performance data.table
    refDT <- validationStats(refModel)[statistic=='AUC',
        list(refAUC=mean(estimate)), by=mDataType]
    modelPerfDT <- merge.data.table(modelPerfDT, refDT, by='mDataType')

    # make basic plot
    plot <- ggplot(modelPerfDT, aes(x=AUC, color=mDataType)) +
        geom_density(aes(fill=mDataType), alpha=0.5) +
        geom_vline(aes(xintercept=refAUC, color=mDataType),
                   linetype="dashed", size=1.5) +
        labs(x = "Model AUCs", y="Density") +
        facet_wrap(~ mDataType, ncol=1, strip.position="right") +
        theme(plot.title=element_text(hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "grey"),
              legend.position="none")

    # user modifications to the plot
    if (!missing(title))
        plot <- plot + ggtitle(title)
    if (!missing(xlab))
        plot <- plot + xlab(xlab)
    if (!missing(ylab))
        plot <- plot + ylab(ylab)
    if (!missing(mDataTypeLabels))
        plot <- plot + facet_wrap(~ mDataType, ncol=1, strip.position='right',
            labeller=labeller(mDataType=mDataTypeLabels))

    return(plot)
})