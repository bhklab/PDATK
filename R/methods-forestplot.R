#' Generate a forest plot from an `S4` object
#'
#' @param object An `S4` object to create a forest plot of.
#' @param ... Allow new parameters to this generic.
#'
#' @return None, draws a forest plot.
#'
#' @export
setGeneric('forestplot', function(object, ...)
    standardGeneric('forestplot'))
#'
#' Render a forest plot from the `validationStats` slot of a `PCOSP` model object.
#'
#' @param object A `PCOSP` model which has been validated with `validateModel`
#'   and therefore has data in the `validationStats` slot.
#' @param stat A `character` vector specifying a statistic to plot.
#' @param groupBy A `character` vector with one or more columns in `validationStats`
#'   to group by. These will be the facets in your
#' @param colourBy A `character` vector specifying the columns in
#'   `validationStats` to colour by.
#' @param ... Fall through arguments to
#'
#' @return A `ggplot2` object
#'
#' @md
#' @import data.table
#' @import ggplot2
#' @export
setMethod('forestplot', signature('PCOSP'),
    function(object, stat, groupBy, ...)
{
    stats <- validationStats(object)[statistic == stat][order(mDataType)]
    plot <- ggplot(stats, aes(y=cohort, x=estimate, xmin=lower, xmax=upper)) +
        geom_pointrange() +
        facet_grid(rows=vars(mDataType), scales='free_y')



})