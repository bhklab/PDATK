# Variable, global to package's namespace.
# This function is not exported to user space and does not need to be documented.
pdatk_options <- settings::options_manager(
    nthread=if (!is.null(getOption('mc.cores'))) {
            getOption('mc.cores')
        } else if ((parallel::detectCores() - 2) > 1) {
            parallel::detectCores() - 2
        } else {
            1
        },
    .allowed=list(nthread=settings::inrange(min=1))
)

#' Reset PDATK Global Options to Default
#'
#' A convenient way to revert the option manager for the PDATK package back
#'   to it's default state.
#'
#' @export
pdatk_reset <- settings::reset(pdatk_options)

#'
#' @param where Optional object to get options from.  If NULL gets the global
#'   options for the package. Default is NULL.
#'
#' @export
setGeneric('manageOptions',
    function(where=NULL, ...) standardGeneric('manageOptions'))
#' @export
setMethod('manageOptions', 'ANY', function(where=NULL, ...) {
    do.call(pdatk_options, c(where, list(...)))
})