#' Generic for Plotting Survival Curves from an `S4` Object
#'
#' @param object An `S4` object to plot survival curves for.
#' @param ... Allow new parameters to be defined on this generic.
#'
#' @return A plot, either via side-effects or as the return value.
#'
#' @md
#' @export
setGeneric('plotSurvivalCurves', function(object, ...) standardGeneric('plotSurvivalCurves'))
#'
#' Plot Survival Curves from a Fit `CoxModel` Object
#'
#' @param object A `CoxModel` object with survival curves fit via the
#'   `trainModel` method.
#' @param byCohort `TRUE` to return a single plot object faceted by `facet.by`,
#'  `FALSE` to get a list of individual survival curves per cohort.
#' @param ... Fall through parameters to `survminer::ggsurvplot` function.
#' @param facet.by What column of the object `modelDT` to use for faceting
#'   the survival plot. Defaults to 'cohorts'. Only used if
#'   `byCohort` is `TRUE`.
#'
#' @return A `ggplot` or `list` of ggplot objects containing the survival
#'   curves for each cohort in the `trainData` slot of the `CoxModel`.
#'
#' @md
#' @importFrom survminer ggsurvplot ggsurvplot_facet
#' @importFrom survival survfit Surv
#' @export
setMethod('plotSurvivalCurves', signature(object='CoxModel'),
    function(object, byCohort=TRUE, ..., facet.by='cohort')
{
    funContext <- .context(1)

    if (!('modelDT' %in% names(models(object))))
        stop(.errorMsg(funContext, 'It looks like your ', class(object), ' ',
            ' object has not been trained. Please use trainModel to ',
            'fit the Cox model before trying to plot survival curves!'))

    modelDT <- models(object)$modelDT
    fitCall <- paste0('survfit(Surv(event=event_occurred, time=survival_time) ~',
        paste0(modelParams(object)[[1]], collapse=' + '), ')')
    overallFit <- eval(str2lang(fitCall), envir=modelDT)

    if (byCohort) {
        plot <- ggsurvplot_facet(overallFit, data=modelDT, facet.by='cohort', ...)
    } else {
        plot <- ggsurvplot(overallFit, data=modelDT, ...)
    }

    return(plot)
})