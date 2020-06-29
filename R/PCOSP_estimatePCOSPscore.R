#' Calculate the probability
#'
##TODO:: HEEWON add function documentation
#'
#' @param valMat A \code{matrix} of expression values for a validation
#'     validation cohort
#' @param selectedModels a \code{list} of models selected using the
#'     `buildPCOSPmodel` function.
#' @param nthread A \code{numeric} vector containing the integer number of
#'     threads to parallelize the calculation across.
#'
#' @return A \code{vector} containing the probabilities per patient
#'
#' @importFrom switchBox SWAP.KTSP.Classify
#' @importFrom BiocParallel bplapply
#' @export
#TODO:: Make this faster!
estimatePCOSPprob <- function(valMat, selectedModels, nthread) {

  # Temporily change number of cores to parallelize over
  opts <- options()
  options("mc.cores"=nthread)
  on.exit(options(opts))

  predictions <- bplapply(selectedModels,
                        function(model, valMat)
                          as.numeric.factor(SWAP.KTSP.Classify(t(valMat), model)),
                        valMat=valMat)

  allPredictions <- do.call(rbind, predictions)
  colnames(allPredictions) <- rownames(valMat)

  predProbabilities <- (length(selectedModels) - colSums(allPredictions)) / length(selectedModels)
  return(predProbabilities)
}