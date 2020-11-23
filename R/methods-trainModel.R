#' Train a Model Based on the Data in an S4 Object
#'
#' @param An `S4` object representing an untrained statistical or machine
#'   learning model.
#'
#' @return The same object with the @model slot populated with the fit model
#'
#' @md
#' @export
setGeneric('trainModel', function(object, ...)
    standardGeneric('trainModel'))
#'
#' Train a PCOSP Model Based on The Data the assay `trainMatrix`.
#'
#' @details This function is parallelized using the nthread option from
#'   the `pdatk_options()` function. Please use `pdatk_options(nthread=n)` to
#'   change from the default value, which is inferred when the package
#'   is loaded. Alternatively, you can pass in your own parallel back-end
#'   using the BPARAM parameter, which will fall through to the bplapply
#'   call inside this function.
#'
#' @param object A `PCOSP` object to train.
#' @param ... Fall through arguments to `BiocParallel::bplapply`
#'
#' @return A `PCOSP` object with the trained model in the `model` slot.
#'
#' @importFrom BiocParallel bplapply
#' @md
#' @export
setMethod('trainModel', signature('PCOSP'), function(object, ...) {
    nthread <- pdatk_option(nthread)


})