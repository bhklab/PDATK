#' Calculate the D and concordance index for a a list of model predictions and
#'     clinical cohorts
#'
#' @param modelProbabilities A \code{list} of per cohort survival probability
#'     predictions.
#' @param clinicalFeatures A \code{list} of per cohort expression data.
#' @param seqCohorts A \code{character} vector with the names of the
#'     cohorts containing sequencing data. Others are assumed to be
#'     from microarray data.
#' @param model A \code{glm} for predicting survival probabilities.
#'
#' @import survival
#' @export
calculateModelDandCindex <- function(modelProbabilities, clinicalFeatures,
                                     seqCohorts, model=1) {

  namesClinical <- lapply(modelProbabilities[[model]],
                          function(cohort)
                            names(cohort$clinical))

  # Subset clinical cohorts to only sample with clinical predictions
  cFeatures <- structure(lapply(seq_along(namesClinical),
                      function(i, namesClinical, clinicalFeatures)
                          clinicalFeatures[[i]][clinicalFeatures[[i]]$ID %in% namesClinical[[i]],],
                      namesClinical=namesClinical,
                      clinicalFeatures=clinicalFeatures),
                      .Names=names(clinicalFeatures))

  ##TODO:: Determine why we invert the probabiltiies here and if the names
  ##    are still correct?
  clinicalProbs <- lapply(modelProbabilities[[model]], function(cohort) 1 - cohort$clinical)
  PCOSPprobs <-lapply(modelProbabilities[[model]], function(cohort) 1 - cohort$PCOSP)

  clinicalStats <- constructMetaEstimatesDF(clinicalProbs, cFeatures, seqCohorts)
  PCOSPstats <- constructMetaEstimatesDF(PCOSPprobs, cFeatures, seqCohorts)

  return(list(
    "clinical"=clinicalStats,
    "PCOSP"=PCOSPstats
  ))
}