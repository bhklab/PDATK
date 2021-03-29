#' Generic for Plotting Survival Curves from an `S4` Object
#' 
#' @param object An `S4` object to plot surival curves for.
#' @param ... Allow new parameters to be defined on this generic.
#' 
#' @return A plot, either via side-effects or as the return value.
#' 
#' @md
#' @export
setGeneric('plotSurvivalCurves', 
    function(object, ...) standardGeneric('plotSurvivalCurves'))

#'
#' @param object A `CoxModel` object with survival curves fit via the 
#'   `trainModel` method.
#' @param arrangePlots `TRUE` to return a single plot object, `FALSE` to
#'   get a list of individual survival curves.
#' @param ... Fall through parameters to `survminer::ggsurvplot` function.
#' 
#' @return A `ggplot` or `list` of ggplot objects containing the survival
#'   curves for each cohort in the `trainData` slot of the `CoxModel`.
#' 
#' @md
#' @importFrom survminer ggsurvplot arrange_ggsurvplots
#' @export
setMethod('plotSurvivalCurves', signature(object='CoxModel'),
    function(object, arrangePlots=TRUE, ...)
{
    funContext <- .context(1)
    
    if (!all(c('modelData', 'survivalFits') %in% names(models(object))))
        stop(.errorMsg(funContext, 'It looks like your ', class(object), ' ',
            ' object has not been trained. Please use trainModel to ',
            'fit the Cox model before trying to plot survival curves!'))

    modelData <- models(object)$modelData
    survivalFits <- models(object)$survivalFits

    survivalPlots <- ggsurvplot(survivalFits)

    if (!arrangePlots) return(survivalPlots)

    ncol <- ceiling(sqrt(length(survivalPlots)))
    nrow <- ceiling(length(survivalPlots) / ncol)

    return(arrange_ggsurvplots(survivalPlots, ncol=ncol, nrow=nrow, 
        print=FALSE))

})