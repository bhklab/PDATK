# Variable, global to package's namespace.
# This function is not exported to user space and does not need to be documented.
PDATKOPTIONS <- settings::options_manager(
    nthread=if (!is.null(getOption('mc.cores'))) {
        getOption('mc.cores')
        } else if ((parallel::detectCores() - 2) > 1) {
            parallel::detectCores() - 2
        } else {
            1
        },
    .allowed=list(nthread=settings::inrange(min=1))
)

# User function that gets exported:

#' Set or get options for my package
#'
#' @param ... Option names to retrieve option values or \code{[key]=[value]}
#'   pairs to set options.
#'
#' @section Supported options:
#' The following options are supported
#' \itemize{
#'  \item{\code{nthread}}{(\code{numeric};1) The number of cores to use for
#'      parallel computations. The default tries to read getOption('mc.cores'),
#'      then if that is null tried to use detectCores() - 2 and if that is
#'      less than 2 falls back to 1.}
#' }
#'
#' @importFrom settings stop_if_reserved
#' @export
pdatk_options <- function(...){
  # protect against the use of reserved words.
  settings::stop_if_reserved(...)
  PDATKOPTIONS(...)
}

#' Reset PDATK Global Options to Default
#'
#' A convenient way to revert the option manager for the PDATK package back
#'   to it's default state.
#'
#' @export
pdatk_reset <- settings::reset(PDATKOPTIONS)