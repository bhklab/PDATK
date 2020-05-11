#' Build PCOSP Model from input data
#'
##TODO:: HEEWON: What does this do? Write a description (two or three senetences max)
#'
#' The script is used for building Pancreatic cancer overall survival predictor
#'   using unique 89 samples profiled using microarray and sequencing platform.
#'
#' @examples
#' # To ensure reproducible results
#' set.seed(1987)
#'
#' # Load the data
#' data(trainingCohorts)
#'
#' # Return the object
#' selectedModels <- buildPCOSPmodels(trainingCohorts, numModels=10, nthread=1)
#'
#' # OR Save the object to disk
#' buildPCOSPmodel(traingCohort, numModels=10, saveDir=tempdir())
#'
##TODO:: Determine where this dataset came from? Is it the ouput of another
#   package? From a publication?
#' @param trainingCohorts A named \code{list} of training cohorts for which
#'   to fit and select PCOSP models.
#' @param saveDir \code{character} A path to a directory to save the model. If you
#'   exclude this the function will return the model object instead.
#' @param nthread \code{integer} The number of threads to parallelize across
#'
#' @return \code{list} Either returns the model object or, is \code{saveDir} is
#'   specified it saves to disk instead and return the path
#'
#' @section Warning: This function uses random numbers; remember to
#'   \code{set.seed()} before running to ensure reproducible results
#'
#' @export
buildPCOSPmodels <- function(trainingCohorts, numModels, nthread, saveDir) {

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

    selectedModels <- .generateTSPmodels(cohortMatrix, cohortMatrixGroups,
                                         numModels)

    # Save to disk or return
    if (!missing(saveDir)) {
        saveRDS(selectedModels,
                file=paste0(file.path(saveDir, 'PCOSPmodels'), '.rds'))
        return(paste0('Saved model to ', saveDir))
    } else {
        return(selectedModels)
    }
}


##TODO:: See if we can refactor part of this to be reused in reshuffleRandomModels
#' @importFrom caret confusionMatrix
#' @importFrom switchBox SWAP.KTSP.Train
#' @importFrom BiocParallel bplapply
.generateTSPmodels <- function(cohortMatrix, cohortMatrixGroups, numModels) {

    trainingDataRowIdxs <- lapply(rep(40, numModels),
                                .randomSampleIndex,
                                labels=cohortMatrixGroups,
                                groups=sort(unique(cohortMatrixGroups)))
    system.time({
    trainedModels <- bplapply(trainingDataRowIdxs,
                              function(idx, data)
                                  SWAP.KTSP.Train(t(data[idx, ]), levels(idx)),
                              data=cohortMatrix)
    })

    testingDataRowIdxs <- lapply(trainingDataRowIdxs,
                            function(idx, rowIdx, labels)
                                structure(setdiff(rowIdx, idx),
                                          .Label=as.factor(
                                              labels[setdiff(rowIdx, idx)])),
                            rowIdx=seq_len(nrow(cohortMatrix)),
                            labels=cohortMatrixGroups)


    predictions <- bplapply(seq_along(testingDataRowIdxs),
                            function(i, testIdxs, data, models)
                                SWAP.KTSP.Classify(t(data[testIdxs[[i]], ]),
                                                   models[[i]]),
                            testIdxs=testingDataRowIdxs,
                            data=cohortMatrix,
                            models=trainedModels
                            )


    confusionMatrices <- bplapply(seq_along(predictions),
                                  function(i, predictions, labels)
                                           confusionMatrix(predictions[[i]],
                                                           levels(labels[[i]]),
                                                           mode="prec_recall"),
                                       predictions=predictions,
                                       labels=testingDataRowIdxs
                            )

    modelStats <- bplapply(confusionMatrices,
                           function(confMat) confMat$byClass)

    balancedAcc <- unlist(bplapply(modelStats,
                              function(model) model[c('Balanced Accuracy')]))

    selectedModels <- trainedModels[which(balancedAcc > 0.60)]
    return(selectedModels)
}

##TODO:: Generalize this to n dimensions
#' Generate a random sample from each group
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
.randomSampleIndex <- function(n, labels, groups) {
    rowIndices <- unlist(mapply(function(x, n, labels) sample(which(labels==x), n, replace=FALSE),
                         x=groups,
                         MoreArgs=list(n=n, labels=labels),
                         SIMPLIFY=FALSE))
    return(structure(rowIndices,
                     .Label=as.factor(labels[rowIndices])))
}
