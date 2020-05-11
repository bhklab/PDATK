#' Build PCOSP models with the group labels randomly shuffled
#'
## TODO:: HEEWON Document this function
#'
#' @param trainingCohorts A named \code{list} of training cohorts for which
#'     to fit and select PCOSP models.
#' @param numModels A \code{numeric} vector containing an integer number of
#'     models to fit
#' @param nthread \code{integer} The number of threads to parallelize across
#' @param filePath A \code{character} vector wtih teh path to a directory
#'   to save the model. If you exclude this the function will return the
#'   model object instead.
#' @param fileNames A \code{character} vector with the file name to save under.
#'   If excluded this function returns the object instead.
#'
#' @return \code{list} Either returns the model object or, is \code{saveDir} is
#'   specified it saves to disk instead and return the path
#'
#' @section Warning: This function uses random numbers; remember to
#'   \code{set.seed()} before running to ensure reproducible results
#'
#' @export
buildRandomLabelShufflingModel <- function(trainingCohorts, numModels, nthread,
                                           saveDir) {

  # Set number of threads to parallelize over
  if (!missing(nthread)) {
    ops <- options()
    options("mc.cores"=nthread)
    on.exit(options(ops))
  }

  # Extract cohorts from trainingCohorts
  seqCohort <- trainingCohorts$icgc_seq_cohort
  arrayCohort <- trainingCohorts$icgc_array_cohort

  # Merged common ICGC seq and array trainingCohorts
  commonData <- mergeCommonData(seqCohort, arrayCohort)

  # Training the model on ICGC seq/array common samples cohort
  cohortMatrix <- convertCohortToMatrix(commonData)
  cohortMatrixGroups <- ifelse(as.numeric.factor(commonData$OS) >= 365, 1, 0)

  selectedModels <- .generateRLSmodels(cohortMatrix, cohortMatrixGroups,
                                       numModels)

  # Save to disk or return
  if (!missing(saveDir)) {
    saveRDS(selectedModels,
            file=paste0(file.path(saveDir, 'RLSmodels'), '.rds'))
    return(paste0('Saved model to ', saveDir))
  } else {
    return(selectedModels)
  }
}

##TODO:: See if we can refactor part of this to be reused in reshuffleRandomModels
#' @importFrom caret confusionMatrix
#' @importFrom switchBox SWAP.KTSP.Train
#' @importFrom BiocParallel bplapply
.generateRLSmodels <- function(cohortMatrix, cohortMatrixGroups, numModels) {

  trainingDataRowIdxs <- lapply(rep(40, numModels),
                                .randomSampleIndexShuffle,
                                labels=cohortMatrixGroups,
                                groups=sort(unique(cohortMatrixGroups)))
  system.time({
    trainedModels <- bplapply(trainingDataRowIdxs,
                              function(idx, data)
                                SWAP.KTSP.Train(t(data[idx, ]), levels(idx)),
                              data=cohortMatrix)
  })

  selectedModels <- trainedModels

  return(selectedModels)
}


##TODO:: Generalize this to n dimensions
#' Generate a random sample from each group and randomly shuffle the labels
#'
#' Returns a list of
#'
#' @param n The sample size
#' @param labels A \code{vector} of the group labels for all rows of the
#'
#' @param groups A vector of group labels for the data to sample from
#' @param numSamples The number of samples to take
#'
#' @return A subset of your object with n random samples from each group in
#'   groups. The number of rows returned will equal the number of groups times
#'   the sample size.
#' @keywords internal
.randomSampleIndexShuffle <- function(n, labels, groups) {
  rowIndices <- unlist(mapply(function(x, n, labels) sample(which(labels==x), n, replace=FALSE),
                              x=groups,
                              MoreArgs=list(n=n, labels=labels),
                              SIMPLIFY=FALSE))
  return(structure(rowIndices,
                   .Label=as.factor(sample(labels)[rowIndices])))
}