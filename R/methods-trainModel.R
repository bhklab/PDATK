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
#' @details This function is parallelized with BiocParallel, thus if you wish
#'   to change the back-end for parallelization, number of threads, or any
#'   other parallelization configuration please pass BPPARAM to bplapply.
#'
#' @param object A `PCOSP` object to train.
#' @param numModels An `integer` specifying the number of models to train.
#'   Defaults to 10. We recommend using 1000+ for good results.
#' @param balancedAccuracy A `float` specifying the balanced accurary required
#'   to consider a model 'top scoring'. Defaults to 0.6. Must be in the
#'   range [0, 1]
#' @param ... Fall through arguments to `BiocParallel::bplapply`
#'
#' @return A `PCOSP` object with the trained model in the `model` slot.
#'
#' @seealso switchBox::SWAP.KTSP.Train BiocParallel::bplapply
#'
#' @importFrom BiocParallel bplapply
#' @md
#' @export
setMethod('trainModel', signature('PCOSP'),
    function(object, numModels=10, balancedAccuracy=0.6, ...)
{
    # Configure local parameters
    opts <- options()
    on.exit(options(opt))

    # Set local seed for random sampling
    set.seed(metadata(object)$randomSeed)

    trainMatrix <- assay(object, 'trainMatrix')
    survivalGroups <- as.factor(colData(object)$survival_group)

    topModels <- .generateTSPmodels(trainMatrix, survivalGroups, numModels,
        balancedAccuracy)


})

##TODO:: See if we can refactor part of this to be reused in reshuffleRandomModels
#' @importFrom caret confusionMatrix
#' @importFrom switchBox SWAP.KTSP.Train
#' @importFrom BiocParallel bplapply
.generateTSPmodels <- function(trainMatrix, survivalGroups, numModels,
    balancedAccurary, ...)
{

    # determine the largest sample size we can take
    sampleSize <- min(sum(survivalGroups == levels(survivalGroups)[1]),
        sum(survivalGroups == levels(survivalGroups[2]))) / 2

    trainingDataColIdxs <- lapply(rep(sampleSize, numModels),
                                .randomSampleIndex,
                                labels=survivalGroups,
                                groups=sort(unique(survivalGroups)))
    system.time({
    trainedModels <- bplapply(trainingDataColIdxs,
                              function(idx, data)
                                  SWAP.KTSP.Train(data[, idx], levels(idx)),
                              data=trainMatrix,
                              ...)
    })

    testingDataColIdxs <- lapply(trainingDataColIdxs,
                                 function(idx, rowIdx, labels)
                                structure(setdiff(rowIdx, idx),
                                          .Label=as.factor(
                                              labels[setdiff(rowIdx, idx)])),
                                 rowIdx=seq_len(ncol(trainMatrix)),
                                 labels=survivalGroups)


    predictions <- bplapply(seq_along(testingDataColIdxs),
                            function(i, testIdxs, data, models)
                                SWAP.KTSP.Classify(data[, testIdxs[[i]]],
                                                   models[[i]]),
                            testIdxs=testingDataColIdxs,
                            data=trainMatrix,
                            models=trainedModels,
                            ...
                            )


    confusionMatrices <- bplapply(seq_along(predictions),
                                  function(i, predictions, labels)
                                           confusionMatrix(predictions[[i]],
                                                           levels(labels[[i]]),
                                                           mode="prec_recall"),
                                       predictions=predictions,
                                       labels=testingDataColIdxs,
                                       ...
                            )

    modelStats <- bplapply(confusionMatrices,
                           function(confMat) confMat$byClass,
                            ...)

    balancedAcc <- unlist(bplapply(modelStats,
                              function(model) model[c('Balanced Accuracy')]),
                              ...)

    selectedModels <- trainedModels[which(balancedAcc > balancedAccuracy)]
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
#'
#' @keywords internal
.randomSampleIndex <- function(n, labels, groups) {
    rowIndices <- unlist(mapply(function(x, n, labels) sample(which(labels==x), n, replace=FALSE),
                         x=groups,
                         MoreArgs=list(n=n, labels=labels),
                         SIMPLIFY=FALSE))
    return(structure(rowIndices,
                     .Label=as.factor(labels[rowIndices])))
}