#' Train a Model Based on the Data in an S4 Object
#'
#' @param object An `S4` object representing an untrained statistical or machine.
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
#' Uses the switchBox SWAP.Train.KTSP function to fit a number of k top scoring
#'   pair models to the data, filtering the results to the best models based
#'   on the specified paramters.
#'
#' @details This function is parallelized with BiocParallel, thus if you wish
#'   to change the back-end for parallelization, number of threads, or any
#'   other parallelization configuration please pass BPPARAM to bplapply.
#'
#' @param object A `PCOSP` object to train.
#' @param numModels An `integer` specifying the number of models to train.
#'   Defaults to 10. We recommend using 1000+ for good results.
#' @param minAccuracy A `float` specifying the balanced accurary required
#'   to consider a model 'top scoring'. Defaults to 0.6. Must be in the
#'   range [0, 1].
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
    function(object, numModels=10, minAccuracy=0.6, ...)
{
    # Configure local parameters
    opts <- options()
    on.exit(options(opts))

    # Set local seed for random sampling
    set.seed(metadata(object)$randomSeed)

    trainMatrix <- assay(object, 'trainMatrix')
    survivalGroups <- as.factor(colData(object)$prognosis)

    topModels <- .generateTSPmodels(trainMatrix, survivalGroups, numModels,
        minAccuracy)

    models(object) <- topModels
    return(object)
})

##TODO:: See if we can refactor part of this to be reused in reshuffleRandomModels
#' @importFrom caret confusionMatrix
#' @importFrom switchBox SWAP.Train.KTSP
#' @importFrom BiocParallel bplapply
#' @importFrom S4Vectors SimpleList
.generateTSPmodels <- function(trainMatrix, survivalGroups, numModels,
    minAccuracy, ...)
{

    # determine the largest sample size we can take
    sampleSize <- min(sum(survivalGroups == levels(survivalGroups)[1]),
        sum(survivalGroups == levels(survivalGroups[2]))) / 2

    # random sample for each model
    trainingDataColIdxs <- lapply(rep(sampleSize, numModels),
                                .randomSampleIndex,
                                labels=survivalGroups,
                                groups=sort(unique(survivalGroups)))

    # train the model
    system.time({
    trainedModels <- bplapply(trainingDataColIdxs,
                              function(idx, data)
                                  SWAP.Train.KTSP(data[, idx], levels(idx)),
                              data=trainMatrix, ...)
    })

    # get the testing data
    testingDataColIdxs <- lapply(trainingDataColIdxs,
                                 function(idx, rowIdx, labels)
                                structure(setdiff(rowIdx, idx),
                                          .Label=as.factor(
                                              labels[setdiff(rowIdx, idx)])),
                                 rowIdx=seq_len(ncol(trainMatrix)),
                                 labels=survivalGroups)

    # make predictions
    predictions <- bplapply(seq_along(testingDataColIdxs),
                            function(i, testIdxs, data, models)
                                SWAP.KTSP.Classify(data[, testIdxs[[i]]],
                                                   models[[i]]),
                            testIdxs=testingDataColIdxs,
                            data=trainMatrix,
                            models=trainedModels,
                            ...)

    # assess the models
    .calcualteConfMatrix <- function(i, predictions, labels) {
        confusionMatrix(predictions[[i]], levels(labels[[i]]),
            mode="prec_recall")
    }

    confusionMatrices <- bplapply(seq_along(predictions), .calcualteConfMatrix,
        predictions=predictions, labels=testingDataColIdxs, ...)

    modelStats <- bplapply(confusionMatrices, `[[`, i='byClass', ...)
    balancedAcc <- unlist(bplapply(modelStats, `[`, i='Balanced Accuracy', ...))

    # sort the models by their accuracy
    keepModels <- balancedAcc > minAccuracy
    selectedModels <- SimpleList(trainedModels[keepModels])
    modelBalancedAcc <- balancedAcc[keepModels]
    selectedModels <- selectedModels[order(modelBalancedAcc, decreasing=TRUE)]
    names(selectedModels) <- paste0('rank', seq_along(selectedModels))

    # capture the accuracies
    mcols(selectedModels)$balancedAcc <-
        modelBalancedAcc[order(modelBalancedAcc, decreasing=TRUE)]
    # capture the function parameters
    metadata(selectedModels) <- list(numModels=numModels,
        minAccuracy=minAccuracy)

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
    .sampleGrp <- function(x, n, labels)
        sample(which(labels==x), n, replace=FALSE)
    rowIndices <- unlist(mapply(.sampleGrp, x=groups,
        MoreArgs=list(n=n, labels=labels), SIMPLIFY=FALSE))
    return(structure(rowIndices,
        .Label=as.factor(labels[rowIndices])))
}