#' Train a Model Based on the Data in an S4 Object
#'
#' @param object An `S4` object representing an untrained statistical or machine.
#'   learning model.
#' @param ... Allow new method to be defined for this generic.
#'
#' @return The same object with the @model slot populated with the fit model
#'
#' @examples
#' data(samplePCOSPmodel)
#' set.seed(getModelSeed(samplePCOSPmodel))
#'
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
#'
#' trainModel(samplePCOSPmodel, numModels=5, minAccuracy=0.6)
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
#'   range \[0, 1\].
#' @param ... Fall through arguments to `BiocParallel::bplapply`. Use this to
#'   configure parallelization options. By default the settings inferred in
#'   `BiocParallel::bpparam()` will be used.
#'
#' @return A `PCOSP` object with the trained model in the `model` slot.
#'
#' @seealso switchBox::SWAP.KTSP.Train BiocParallel::bplapply
#'
#' @examples
#' data(samplePCOSPmodel)
#'
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
#'
#' set.seed(getModelSeed(samplePCOSPmodel))
#' trainModel(samplePCOSPmodel, numModels=2, minAccuracy=0.6)
#'
#' @md
#' @importFrom BiocParallel bplapply
#' @import BiocGenerics
#' @export
setMethod('trainModel', signature('PCOSP'),
    function(object, numModels=10, minAccuracy=0.6, ...)
{

    assays <- assays(object)

    # Lose uniqueness of datatypes here
    trainMatrix <- do.call(rbind, assays)
    survivalGroups <- factor(colData(object)$prognosis, levels=c('good', 'bad'))

    topModels <- .generateTSPmodels(trainMatrix, survivalGroups, numModels,
        minAccuracy, ...)

    models(object) <- topModels
    # Add additiona model paramters to metadta
    metadata(object)$modelParams <- c(metadata(object)$modelParams,
        list(numModels=numModels, minAccuracy=minAccuracy))
    return(object)
})

#' @importFrom caret confusionMatrix
#' @importFrom switchBox SWAP.KTSP.Train SWAP.KTSP.Classify
#' @importFrom BiocParallel bplapply
#' @importFrom S4Vectors SimpleList
#' @keywords internal
.generateTSPmodels <- function(trainMatrix, survivalGroups, numModels,
    minAccuracy, sampleFUN=.randomSampleIndex, ...)
{

    # determine the largest sample size we can take
    sampleSize <- min(sum(survivalGroups == levels(survivalGroups)[1]),
        sum(survivalGroups == levels(survivalGroups)[2])) / 2

    # random sample for each model
    trainingDataColIdxs <- lapply(rep(sampleSize, numModels),
                                sampleFUN,
                                labels=survivalGroups,
                                groups=sort(unique(survivalGroups)))

    # train the model
    system.time({
    trainedModels <- bplapply(trainingDataColIdxs,
        function(idx, data)
            switchBox::SWAP.KTSP.Train(data[, idx], levels(idx)),
        data=trainMatrix, ...)
    })

    # get the testing data
    testingDataColIdxs <- lapply(trainingDataColIdxs,
        function(idx, rowIdx, labels)
            structure(setdiff(rowIdx, idx), .Label=as.factor(
                labels[setdiff(rowIdx, idx)])),
        rowIdx=seq_len(ncol(trainMatrix)),
        labels=survivalGroups)

    # make predictions
    predictions <- bplapply(seq_along(testingDataColIdxs),
        function(i, testIdxs, data, models) {
                    switchBox::SWAP.KTSP.Classify(data[, testIdxs[[i]]],
                               models[[i]])
        },
        testIdxs=testingDataColIdxs,
        data=trainMatrix,
        models=trainedModels,
        ...)

    # assess the models
    .calculateConfMatrix <- function(i, predictions, labels) {
        caret::confusionMatrix(predictions[[i]], levels(labels[[i]]),
            mode="prec_recall")
    }

    confusionMatrices <- bplapply(seq_along(predictions), .calculateConfMatrix,
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
    metadata(selectedModels) <- list(numModels=numModels)

    return(selectedModels)
}

##TODO:: Generalize this to n dimensions
#' Generate a random sample from each group
#'
#' Returns a list of equally size random samples from two or more sample
#'   groupings.
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
        sample(which(labels == x), n, replace=FALSE)
    rowIndices <- unlist(mapply(.sampleGrp, x=groups,
        MoreArgs=list(n=n, labels=labels), SIMPLIFY=FALSE))
    return(structure(rowIndices,
        .Label=as.factor(labels[rowIndices])))
}


# ---- RLSModel methods

#' @inherit trainModel,PCOSP-method
#'
#' @param minAccuracy This parameter should be set to zero, since we do
#'   not expect the permuted models to perform well. Setting this higher will
#'   result in an ensemble with very few models included.
#'
#' @examples
#' data(sampleRLSmodel)
#' set.seed(getModelSeed(sampleRLSmodel))
#'
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
#'
#' trainedRLSmodel <- trainModel(sampleRLSmodel, numModels=2)
#'
#' @md
#' @export
setMethod('trainModel', signature('RLSModel'),
    function(object, numModels=10, minAccuracy=0, ...)
{

    assays <- assays(object)

    # Lose uniqueness of datatypes here
    trainMatrix <- do.call(rbind, assays)
    survivalGroups <- factor(colData(object)$prognosis, levels=c('good', 'bad'))

    RLSmodels <- .generateTSPmodels(trainMatrix, survivalGroups, numModels,
        sampleFUN=.randomSampleIndexShuffle, minAccuracy=minAccuracy, ...)

    models(object) <- RLSmodels
    # Add additional model parameters to metadata
    metadata(object)$modelParams <- c(metadata(object)$modelParams,
        list(numModels=numModels))
    return(object)
})

##TODO:: Generalize this to n dimensions
#' Generate a random sample from each group and randomly shuffle the labels
#'
#' @param n The sample size
#' @param labels A \code{vector} of the group labels for all rows of the
#' @param groups A vector of group labels for the data to sample from
#'
#' @return A subset of your object with n random samples from each group in
#'   groups. The number of rows returned will equal the number of groups times
#'   the sample size.
#'
#' @md
#' @keywords internal
#' @noRd
.randomSampleIndexShuffle <- function(n, labels, groups) {
  rowIndices <- unlist(mapply(
    function(x, n, labels) sample(which(labels==x), n, replace=FALSE),
                              x=groups,
                              MoreArgs=list(n=n, labels=labels),
                              SIMPLIFY=FALSE))
  return(structure(rowIndices,
                   .Label=as.factor(sample(labels)[rowIndices])))
}

#' Train a RGAModel Based on the Data in the assays slot.
#'
#' Uses the switchBox SWAP.Train.KTSP function to fit a number of k top scoring
#'   pair models to the data, filtering the results to the best models based
#'   on the specified paramters.
#'
#' @details This function is parallelized with BiocParallel, thus if you wish
#'   to change the back-end for parallelization, number of threads, or any
#'   other parallelization configuration please pass BPPARAM to bplapply.
#'
#' @param object A `RGAmodel` object to train.
#' @param numModels An `integer` specifying the number of models to train.
#'   Defaults to 10. We recommend using 1000+ for good results.
#' @param minAccuracy A `float` specifying the balanced accuracy required
#'   to consider a model 'top scoring'. Defaults to 0. Must be in the
#'   range 0 to 1.
#' @param ... Fall through arguments to `BiocParallel::bplapply`.
#'
#' @return A `RGAModel` object with the trained model in the `model` slot.
#'
#' @seealso `switchBox::SWAP.KTSP.Train` `BiocParallel::bplapply`
#'
#' @examples
#' data(sampleRGAmodel)
#' set.seed(getModelSeed(sampleRGAmodel))
#'
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
#'
#' trainedRGAmodel <- trainModel(sampleRGAmodel, numModels=2, minAccuracy=0)
#'
#' @md
#' @importFrom BiocParallel bplapply
#' @importFrom SummarizedExperiment assays assays<-
#' @importFrom S4Vectors metadata
#' @export
setMethod('trainModel', signature('RGAModel'),
    function(object, numModels=10, minAccuracy=0, ...)
{

    assays <- assays(object)

    # Lose uniqueness of datatypes here
    trainMatrix <- do.call(rbind, assays)
    survivalGroups <- factor(colData(object)$prognosis, levels=c('good', 'bad'))

    RGAmodels <- .generateTSPmodels(trainMatrix, survivalGroups, numModels,
        minAccuracy, ...)

    geneNames <- rownames(trainMatrix)
    for (i in seq_along(RGAmodels)) {
        RGAmodels[[i]]$TSPs[, 1] <- sample(geneNames, nrow(RGAmodels[[i]]$TSPs))
        RGAmodels[[i]]$TSPs[, 2] <- sample(geneNames, nrow(RGAmodels[[i]]$TSPs))
    }

    models(object) <- RGAmodels
    # Add additional model paramters to metadata
    metadata(object)$modelParams <- c(metadata(object)$modelParams,
        list(numModels=numModels, minAccurary=minAccuracy))
    return(object)
})


# ---- ClinicalModel Methods

#' Fit a GLM Using Clinical Predictors Specified in a `ClinicalModel` Object.
#'
#'
#' @param object A `ClinicalModel` object, with survival data for the model
#'   in the colData slot.
#' @param ... Fall through parameters to the [`stats::glm`] function.
#' @param family Argument to the family parameter of [`stats::glm`]. Defaults
#'   to `binomial(link='logit')`. This parameter must be named.
#'@param na.action Argument to the na.action paramater of [`stats::glm`].
#'   Deafults to 'na.omit', dropping rows with NA values in one or more of the
#'   formula variables.
#'
#' @return A `ClinicalModel` object with a `glm` object in the models slot.
#'
#' @examples
#' data(sampleClinicalModel)
#' set.seed(getModelSeed(sampleClinicalModel))
#'
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
#'
#' trainedClinicalModel <- trainModel(sampleClinicalModel)
#'
#' @md
#' @importFrom stats glm as.formula binomial
#' @importFrom S4Vectors SimpleList
#' @export
setMethod('trainModel', signature(object='ClinicalModel'),
    function(object, ..., family=binomial(link='logit'), na.action=na.exclude)
{

    survivalData <- colData(object)
    survivalData$prognosis <- factor(survivalData$prognosis,
        levels=c('good', 'bad'))
    formula <- as.formula(metadata(object)$modelParams$formula)

    # use eval substitute to ensure the formula is captured in the glm call
    model <- eval(substitute(glm(formula, data=survivalData,
        family=family, na.action=na.action, ...)))

    models(object) <- SimpleList(glm=model)
    return(object)

})

# ---- GeneFuModel Methods

## TODO: Update constructor to allow GeneFuModel data to be input into the object
#' Train a GeneFuModel Object
#'
#' @param object A `GeneFuModel` object to train.
#'
#' @return An error message, since we have not finished implementing this
#'   functionality yet.
#'
#' @md
#' @export
setMethod('trainModel', signature(object='GeneFuModel'), function(object) {

  funContext <- .context(1)
  stop(.errorMsg(funContext, 'Unfortunately we have not implemented model ',
                 'training for a GeneFuModel yet!'))

})


# ---- ConsensusMetaclusteringModel

#' Train A ConsensusMetaclusteringModel
#' 
#' @description
#' Since consensus clustering is an unsuperived learning method, there isn't
#'   really a 'training step' per se. Instead this method computes the consensus
#'   clusters and stores the results in the `models` slot.
#' 
#' @param object A `ConsensusMetaclusteringModel` to train.
#' @param maxK The maximum number of clusters to test. Defaults to 5.
#' @param reps How many random samples should clustering be repeated on? Default
#'   is 10, but 1000+ is recommended for real world use.
#' @param ... Fall through parameters to `BiocParallel::bplapply`. This can
#'   be used to customize your parallelization using BPPARAM or to pass
#'   additional arguments to `ConsensusClusterPlus`.
#' @param distance The distance method to use. Defaults to 'pearson'. See 
#'   `?ConsensusClusterPlus::ConsensusClusterPlus` for more options.
#' @param clusterAlg The clustering algorithm to use. Defaults to 'hc'. See
#'   `?ConesnsusClusterPLus::ConsensusClusterPlus` for more options.
#' @param plot An optional path to output the plots generated by each
#'   call to `ConsensusClusterPlus::ConsensusClusterPlus`. Default is NULL,
#'   which suppresses all plots, otherwise passed to the clustering function.
#' 
#' @return The `ConsensusMetaclusteringModel` with the clustering results in the
#'   `models` slot.
#' 
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @importFrom BiocParallel bplapply
#' 
#' @md
#' @export
setMethod('trainModel', signature(object='ConsensusMetaclusteringModel'), 
    function(object, maxK=5, reps=10, distance='pearson', clusterAlg='hc', 
        plot=NULL, ...)
{
    if (is.null(plot)) pdf(NULL)
    
    clusteringResults <- bplapply(assays(trainData(object)), 
        FUN=ConsensusClusterPlus, maxK=maxK, reps=reps, distance=distance, 
        clusterAlg=clusterAlg, plot=plot, seed=modelParams(object)$randomSeed,
        ...)
    models(object) <- clusteringResults

    if (is.null(plot)) dev.off()

    modelParams(object) <- c(modelParams(object), list(maxK=maxK, reps=reps, 
        distance=distance, clusteringAlg=clusterAlg), list(...))

    return(object)
})

# ---- NetworkCommunitySearchModel

#' Train a NetworkCommunitySearchModel
#' 
#' @param object An `NCSModel` object, created from a 
#'   `ConsensusMetaclusteringModel`.
#' @param alpha A `float` specifying the significance level for cluster
#'   reproducibility. Default is 0.05.
#' @param minRepro A `float` specifying the minimum in-group proportion (IGP) 
#'   for a cluster to be included in the metacluster labels. Default is 0.5.
#' @param minCor A `float` specifying the minimum correlation between a 
#'   centroid and assay cluster to be included in the metacluster labels. 
#'   Default is 0.0.
#' 
#' @md
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder setnames
#' @export
setMethod('trainModel', signature(object='NCSModel'), function(object, 
    alpha=0.05, minRepro=0.5, minCor=0.0) 
{
    modelParams(object) <- SimpleList(list(
        alpha=alpha,
        minRepro=minRepro,
        minCor=minCor
    ))

    models(object)$networkEdges <- models(object)$networkEdges[
        p_value <= alpha &
        ingroup_proportion >= minRepro &
        cor_threshold >= minCor,
    ]

    return(object)
})

# ---- CoxModel-methods

#' 
#' 
#' @param object A `CoxModel` object to fit models for.
#' 
#' @md
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder setnames
#' @importFrom survival coxph survfit survdiff
#' @importFrom stats pchisq
#' @importFrom S4Vectors SimpleList
#' @export
setMethod("trainModel", signature(object='CoxModel'), 
    function(object) 
{
    colDataL <- lapply(experiments(trainData(object)), colData)
    colDataL <- lapply(colDataL, as.data.table, keep.rownames=TRUE)
    colDataL <- lapply(colDataL, na.omit)
    survivalPredictor <- modelParams(object)$survivalPredictor
    if (length(survivalPredictor) > 1) {
        warning('Multiple predictors have not been tested yet...')
        survivalPredictor <- paste0(survivalPredictor, collapse=' + ')
    }
    modelFormula <- paste0('Surv(event=event_occurred, time=survival_time) ~ ', 
        survivalPredictor)
    coxModels <- survivalFits <- survivalDiffs <- vector('list', length(colDataL))
    FUNs <- c('coxph', 'survfit', 'survdiff')
    modelData <- colDataL
    for (i in seq_along(modelData)) {
        callStrings <- paste0(FUNs, '(', modelFormula, ', modelData$', 
            names(colDataL)[i], ')')
        coxModels[[i]] <- eval(str2lang(callStrings[1]))
        survivalFits[[i]] <- eval(str2lang(callStrings[2]))
        survivalDiffs[[i]] <- eval(str2lang(callStrings[3]))
    }
    chisq <- vapply(survivalDiffs, `[[`, 'chisq', FUN.VALUE=numeric(1))
    degFreedom <- vapply(survivalDiffs, function(x) length(x$n) - 1, numeric(1))
    chisqPval <- 1 - pchisq(chisq, degFreedom)

    for (i in seq_along(colDataL)) colDataL[[i]][, cohort := names(colDataL)[i]]
    modelDT <- rbindlist(colDataL, fill=TRUE)
    commonCols <- Reduce(intersect, lapply(colDataL, colnames))
    modelDT <- modelDT[, .SD, .SDcols=commonCols]

    models(object) <- SimpleList(list(
        coxModels=coxModels,
        survivalFits=survivalFits,
        survivalDiffs=survivalDiffs,
        chisqPvalues=chisqPval,
        modelData=modelData,
        modelDT=modelDT
    ))

    return(object)
})