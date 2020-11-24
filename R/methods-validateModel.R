#' Perform Validation on an `S4` Object Respresenting a Trained Model
#'
#' @param object An `S4` object
#'
#' @md
#' @export
setGeneric('validateModel', function(object, validationData, ...)
    standardGeneric('validateModel'))
#'
#' Evaluate the Performance of a List of Trained KTSP Models from a PCOSP
#'   Object
#'
#' @param object A `PCOSP` model which has been trained using `trainModel`.
#' @param validationData A `CohortList` containing one or more
#'   `SurvivalExperiment`s. The first assay in each `SurvivalExperiment` will
#'   be classified using all top scoring KTSP models in `models(object)`.
#' @param ... Fallthrough arguments to `BiocParallel::bplapply`, use this to
#'   configure the parallelization settings for this function. For example
#'   to specify BPARAM.
#'
#' @seealso BiocParallel::bplapply switchBox::SWAP.KTSP.Classify
#'
#' @return
#'
#' @importFrom BiocParallel bplapply
#' @importFrom switchBox SWAP.KTSP.Classify
#' @md
#' @export
setMethod('validateModel', signature(object='PCOSP',
    validationData='CohortList'), function(object, validationData, ...)
{


})