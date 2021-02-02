#' Make A Bar Plot Comparing Perforamnce Between Two `S4` Obects Representing
#'   Mathematical Models.
#'
#' @param model1 An`S4` object containing results of a mathematical model
#' @param model2 An `S4` object containing results of a different mathematical
#'   model, but with the same or overalapping samples.
#'
#' @return A bar plot comparing some aspect of model1 and model2
#'
#' @md
#' @export
setGeneric("barPlotModelComparison", function(model1, model2, ...)
    standardGeneric('barPlotModelComparison'))
#'
#' Make a Bar Plot Comparison Model Performance Between a ClinicalMode
#'   and a PCOSP, RLSModel or RGAModel object.
#'
#' @param model1 A `ClinicalModel` object.
#' @param model2 A `PCOSP` or `RLSModel` or `RGAModel` object.
#' @param stat A `character` vector specifying which statistic to compare the
#'   models using. Options are 'AUC', 'D_index' or 'concordance_index'.
#'
#' @return A `ggplot2` object showing a barplot coloured by the model and
#'   comparing the stat between all cohorts that both models were validated
#'   with.
#'
#' @md
#' @include classUnions.R
#' @importFrom ggplot2 ggplot geom_col aes
#' @export
setMethod('barPlotModelComparison', signature(model1='ClinicalModel',
    model2='PCOSP_or_RLS_or_RGA'), function(model1, model2, stat, ...)
{

    valStats1 <- validationStats(model1)
    valStats1$model <- class(model1)
    valStats2 <- validationStats(model2)
    valStats2$model <- class(model2)

    valStats <- rbind(valStats1, valStats2[cohort %in% valStats1$cohort])
    valStats <- valStats[statistic == stat & !isSummary, ]

    plot <- ggplot(valStats, aes(x=cohort, y=estimate)) +
        geom_col(aes(fill=model), position='dodge')

    plot

})