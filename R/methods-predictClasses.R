#' Predict Classes for New Data Based on a Train Classifier Model
#'
#' @param object An `S4` object containing data to predict classes from.
#' @param model An `S4` object containing one or more trained classification
#'   models.
#' @param ... Allow further parameters to be defined on this generic.
#'
#' @return The `S4` object with class predictions added to the metadata.
#'
#' @examples
#' data(sampleTrainedPCOSPmodel)
#' data(samplePCSIsurvExp)
#'
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
#'
#' # Make predictions
#' PCOSPpredSurvExp <- predictClasses(samplePCSIsurvExp,
#'   model=sampleTrainedPCOSPmodel)
#' head(colData(PCOSPpredSurvExp))
#'
#' @md
#' @export
setGeneric('predictClasses',
    function(object, model, ...) standardGeneric('predictClasses'))
#'
#' Predict Survival Prognosis Classes and Risk Scores for A `CohortList` Using
#'   a `PCOSP`, `RLSModel` or `RGAModel` object.
#'
#' @param object A `SurvivalExperiment` object to predict classes for.
#' @param model A trained `PCOSP` model to use for predicting classes.
#' @param ... Fall through arguments to `BiocParallel::bplapply` for configuring
#'   parallelization settings.
#'
#' @seealso [`BiocParallel::bplapply`], [`switchBox::SWAP.KTSP.Classify`]
#'
#' @return A `SurvivalExperiment` with the predictions in its metadata and
#'   a column in colData, `prob_good_survival`, which contains the proportion
#'   of models which predicted good prognosis for each sample.
#'
#' @examples
#' data(sampleTrainedPCOSPmodel)
#' data(samplePCSIsurvExp)
#'
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
#'
#' # Make predictions
#' PCOSPpredSurvExp <- predictClasses(samplePCSIsurvExp,
#'   model=sampleTrainedPCOSPmodel)
#' head(colData(PCOSPpredSurvExp))
#'
#' @md
#' @importFrom SummarizedExperiment assays assay
#' @importFrom CoreGx .warnMsg .errorMsg
#' @importFrom BiocParallel bplapply
#' @importFrom S4Vectors metadata
#' @export
setMethod('predictClasses', signature(object='SurvivalExperiment',
    model='PCOSP_or_RLS_or_RGA'), function(object, model, ...)
{
    funContext <- .context(1)

    # drop NA samples, they mess with calculating statistics
    keepSamples <- rownames(na.omit(colData(object)))
    if (!all(colnames(object) %in% keepSamples)) {
        warning(.warnMsg(funContext, 'Dropped samples with NA survival data!'))
    }
    object <- object[, keepSamples]

    modelList <- models(model)
    if (length(modelList) < 1)
        stop(.errorMsg(funContext, 'There are no models in the PCOSP model ',
            'passed as the model argument. Please ensure you train your model',
            ' with `trainModel` before attempting to predict classes with it.'))

    assayData <- assays(object)
    if (length(assayData) > 1)
        warning(.warnMsg(funContext, 'Please ensure your prediction ',
                         'data only has one assay! Ignoring ',
                         'all but the first assay!'))

    assayMatrix <- assayData[[1]]

    predictionList <- bplapply(modelList, FUN=SWAP.KTSP.Classify,
                               inputMat=assayMatrix)
    # convert factor to character
    predictionListChar <- lapply(predictionList, as.character)
    predictions <- BiocGenerics::do.call(rbind, predictionListChar)
    colnames(predictions) <- colnames(assayMatrix)

    metadata(object)[[paste0(class(model)[1], 'predictions')]] <- predictions
    metadata(object)[[paste0(class(model)[1], 'params')]] <-
        metadata(model)$modelParams

    colData(object)[paste0(class(model)[1], '_prob_good')] <-
        colSums(predictions == 'good') / nrow(predictions)
    colData(object)$prognosis <-
        ifelse(colData(object)$survival_time > 365, 'good', 'bad')

    return(object)
})
#'
#' Predict Survival Prognosis Classes and Risk Scores for A `CohortList` Using
#'   a `PCOSP`, `RLSModel` or `RGAModel` object.
#'
#' @param object A `CohortList` with `SurvivalExperiment`s to predict classes
#'   for.
#' @param model A trained `PCOSP` model to use for predicting classes.
#' @param ... Fall through arguments to `BiocParallel::bplapply` for configuring
#'   parallelization settings.
#'
#' @return A `CohortList` with the model predictions attached to each
#'   `SurvivalExperiment` in the metadata slot and the `prob_good_survival`
#'   column added to the colData slot.
#'
#' @examples
#' data(sampleTrainedPCOSPmodel)
#' data(sampleCohortList)
#'
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
#'
#' # Make predictions
#' PCOSPpredCohortList <- predictClasses(sampleCohortList[seq_len(2)],
#'   model=sampleTrainedPCOSPmodel)
#' head(colData(PCOSPpredCohortList[[1]]))
#'
#' @md
#' @importFrom S4Vectors endoapply mcols metadata
#' @export
setMethod('predictClasses', signature(object='CohortList',
    model='PCOSP_or_RLS_or_RGA'), function(object, model, ...)
{
    predictionResults <- endoapply(object, predictClasses, model=model, ...)
    mcols(predictionResults)$hasPredictions <- TRUE
    metadata(predictionResults)$predictionModel <- model
    return(predictionResults)
})

# ---- ClinicalModel methods

#' Predict Survival Prognosis Classes and Risk Scores for A `SurvivalModel` Using
#'   a `ClinicalModel` Object.
#'
#' @param object A `SurvivalExperiment` object with the correct columns in
#'   `colData` to match the formula for the `ClinicalModel` object.
#' @param model A trained `ClinicalModel` object, as return by `trainModel`.
#' @param ... Fall through parameters to [`stats::predict`].
#' @param na.action The `na.action` paramter passed to [`stats::predict.glm`].
#' @param type The `type` parameter passed to [`stats::predict.glm`]
#'
#' @return A `SurvivalExperiment` with the model predictions in the colData
#'   slot as clinical_prob_good.
#'
#' @examples
#' data(sampleClinicalModel)
#' data(samplePCSIsurvExp)
#'
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
#'
#' # Train Model
#' trainedClinicalModel <- trainModel(sampleClinicalModel)
#'
#' # Make predictions
#' ClinicalPredSurvExp <- predictClasses(samplePCSIsurvExp,
#'   model=trainedClinicalModel)
#' head(colData(ClinicalPredSurvExp))
#'
#' @md
#' @importFrom stats predict glm as.formula
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom CoreGx .errorMsg .warnMsg
#' @importFrom S4Vectors metadata
#' @export
setMethod('predictClasses', signature(object='SurvivalExperiment',
    model='ClinicalModel'), function(object, model, ..., na.action='na.exclude',
        type='response')
{
    funContext <- .context(1)
    # check that the formula is valid and the variables are in the training data
    formula <- as.formula(metadata(model)$modelParams$formula)
    formulaCols <- as.character(formula[seq(2, 3)])
    # split the formula into a vector where each variable is an item
    formulaCols <- unlist(strsplit(formulaCols,
        split='[\\s]*[\\+\\-\\~\\=\\*][\\s]*', perl=TRUE))
    hasFormulaCols <- formulaCols %in% colnames(colData(object))
    if (!all(hasFormulaCols))
        stop(.errorMsg(funContext, 'The columns ', formulaCols[!hasFormulaCols],
            ' are missing from the colData slot of the training data. ',
            'Please only specify valid column names in colData to the formula!'))
    if (length(models(model)) > 1) warning(.warnMsg(funContext, 'There is more ',
        'than one model in your ClinicalModel. Only using the first one...'))

    # Skip rows with levels that aren't in the model; prevents predict.glm for
    #   breaking if there are new levels prediction data
    modelFactorLevels <- models(model)$glm$xlevels
    keepRows <- rep(TRUE, nrow(colData(object)))
    for (name in names(modelFactorLevels)) {
        keep <- colData(object)[[name]] %in% modelFactorLevels[[name]]
        keepRows <- keepRows & keep
    }
    if (!all(keepRows))
        warning(.warnMsg(funContext, 'Rows ', paste0(which(!keepRows),
            collapse=', '), ' have levels that are not in the model, skipping ',
            'these rows...'))
    # Calculate survival probabiltiies
    predictions <- predict(models(model)[[1]], colData(object)[keepRows, ],
        ..., na.action=na.action, type=type)
    # Update the `SurvivalExperiment` object with the predicted probabilities
    colData(object)$clinical_prob_good <- NA
    colData(object)$clinical_prob_good[
        rownames(colData(object)) %in% names(predictions)] <- predictions
    metadata(object)$GLMpredictions <- matrix(
        ifelse(colData(object)$clinical_prob_good > 0.5, 'good', 'bad'),
        byrow=TRUE, nrow=1, dimnames=list('glm', rownames(colData(object))))
    metadata(object)$GLMparams <- metadata(model)$modelParams
    return(object)
})
#'
#' Use a Clinical GLM to Predict Classes for a `CohortList` of
#'   `SurvivalExperment` Objects.
#'
#' @param object A `CohortList` with `SurvivalExperiment`s to predict classes
#'   for. The colData slot in ALL `SurvivalExperiment`s must have column names
#'   which match the formula in the model object.
#' @param model A trained `ClinicalModel` object, as return by `trainModel`.
#' @param ... Fall through parameters to [`stats::predict`].
#' @param na.action The `na.action` paramter passed to [`stats::predict.glm`].
#' @param type The `type` parameter passed to [`stats::predict.glm`]
#'
#' @return A `CohortList` with the model predictions in the colData
#'   slot as clinical_prob_good for each `SurvivalExperiment`, and the
#'   model in the metadata as predictionModel.
#'
#' @examples
#' data(sampleClinicalModel)
#' data(sampleCohortList)
#'
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
#'
#' # Train Model
#' trainedClinicalModel <- trainModel(sampleClinicalModel)
#'
#' # Make predictions
#' ClinicalPredCohortList <- predictClasses(sampleCohortList[c('PCSI', 'TCGA')],
#'   model=trainedClinicalModel)
#' head(colData(ClinicalPredCohortList[[1]]))
#'
#' @md
#' @importFrom S4Vectors endoapply mcols metadata
#' @export
setMethod('predictClasses', signature(object='CohortList',
    model='ClinicalModel'), function(object, model, ..., na.action='na.exclude',
        type='response')
{
    predictionResults <- endoapply(object, predictClasses, model=model, ...,
                                   na.action=na.action, type=type)
    mcols(predictionResults)$hasPredictions <- TRUE
    metadata(predictionResults)$predictionModel <- model
    return(predictionResults)
})


# ---- GeneFuModel methods

#' Use a Gene Signature Based Prediciton Model from the `genefu` Package
#'   to Predict Signature Scores for Each Sample.
#'
#' @param object A `SurvivalExperiment` to predict classes for.
#' @param model A `GeneFuModel` object to predict classes with.
#' @param ... Fall through parameters to `genefu::sig.score`.
#' @param annot A named parameter with annotations mapping from the gene
#'   identifiers in the genefu model.
#'
#' @details
#' A signature score should be interpreted as unit-less continuous risk predictor.
#'
#' @return The `SurvivalExperiment` passed to the `object` argument with
#'   the `genefu_score` column added to the objects `colData` slot.
#'
#' @md
#' @importFrom genefu sig.score
#' @importFrom CoreGx .warnMsg
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors metadata
#' @export
setMethod('predictClasses', signature(object='SurvivalExperiment',
    model='GeneFuModel'), function(object, model, ..., annot=NA)
{
    funContext <- .context(1)

    if (length(models(model)) > 1) {
        warning(.warnMsg(funContext, 'There is more than one model in the ',
            'models slot of the GeneFuModel. Currently only single model ',
            'predictions are supported, ignoring other models.'))
    }

    # Calculate survival releative risks (not probabiltiies)
    predictions <- genefu::sig.score(x=models(model)[[1]],
                                     t(assay(object, 1)), annot=annot, ...)
    # Add a column to colData for each model in the GeneFuModel
    colData(object)$genefu_score <- predictions$score
    metadata(object)$genefuPredictions <- 'GeneFu models do not make explicit
        category predictions, instead they return a signature score which should
        be treated as a unitless continuous risk predictor. To make a class
        prediction, the user must select a cut-off, which is highly subjective.'
    metadata(object)$genefuParams <- c(metadata(model)$modelParams,
        predictions[-1])

    return(object)
})
#'
#' Use a Gene Signature Based Prediciton Model from the `genefu` Package
#'   to Predict Signature Scores for Each Sample
#'
#' @param object A `CohortList` with `SurvivalExperiment`s to predict classes
#'   for.
#' @param model A trained `GeneFuModel` object.
#' @param ... Fall through parameters to [`genefu::sig.score`].
#' @param annot The `annot` parameter passed to [`genefu::sig.score`].
#'   Defaults to NA, which assumes your assay rowname match the gene labels
#'   in the model.
#'
#' @return A `CohortList` with the model predictions in the colData
#'   slot as genefu_score for each `SurvivalExperiment`, and the
#'   model in the metadata as predictionModel.
#'
#' @md
#' @importFrom S4Vectors endoapply mcols metadata
#' @export
setMethod('predictClasses', signature(object='CohortList',
    model='GeneFuModel'), function(object, model, ..., annot=NA)
{
    predictionResults <- endoapply(object, predictClasses, model=model, ...,
                                   annot=annot)
    mcols(predictionResults)$hasPredictions <- TRUE
    metadata(predictionResults)$predictionModel <- model
    return(predictionResults)
})


# ---- ConsensusMetaclusteringModel

#' Compute the Optimal Clustering Solution for a Trained 
#'   ConsensusMetaclusteringModel
#' 
#' Compute the optimal clustering solution out of possibilities generated
#'   with trainModel. Assigns the cluster labels to the `MultiAssayExperiment`
#'   object.
#' 
#' @param object A `MutliAssayExperiment` object
#' @param subinterval A `numeric` vector of two float values, the first
#'   being the lower and second being the upper limit of the subinteral
#'   to compare cluster ambiguity over. Default is c(0.1, 0.9), i.e. comparing
#'   the 10th and 90th percentile of cluster consensus to calculate the 
#'   ambiguity of a given clustering solution. This is the value used to
#'   selected the optimal K value from the potential solutions for each
#'   assay in the training data.
#' 
#' @return A `object` `ConsensusMetaclusteringModel`, with class predictions
#'   assigned to the colData of `trianData`
#' 
#' @md
#' @export
setMethod('predictClasses', signature(object='ConsensusMetaclusteringModel'), 
    function(object, subinterval=c(0.1, 0.9))
{
    ## TODO:: Refactor into some helpers so function is shorter
    
    funContext <- .context(1)
    if (length(models(object)) < 1) stop(.errorMsg(funContext, 'The ',
        class(object)[1], ' object does not have any clustering results.
        Please run trainModel first, before trying to predictClasses'))
    
    # -- Find optimal K value
    assayClusters <- models(object)
    .getFromClusteringResults <- function(x, i='consensusMatrix') lapply(x[-1], 
        FUN=`[[`, i=i)
    consensusMatrices <- lapply(assayClusters, FUN=.getFromClusteringResults)
    .getLowerTris <- function(x) lapply(x, lower.tri)
    lowerTris <- lapply(consensusMatrices, .getLowerTris)
    .subsetLowerTris <- function(x, lowerTris) mapply(`[`, x, lowerTris, 
        SIMPLIFY=FALSE)
    lowerTriConMatrices <- mapply(.subsetLowerTris, consensusMatrices, 
        lowerTris, SIMPLIFY=FALSE)
    # Make and empirical cummulative densitiy function for each consensus matrix
    .getECDFS <- function(x) lapply(x, FUN=ecdf)
    assayECDFS <- lapply(lowerTriConMatrices, FUN=.getECDFS)

    # Compare the 1st and 9th deciles as a metric of cluster ambiguity
    .calcPropAmbiguousClusters <- function(ecdfs, subinterval) {
        vapply(ecdfs, function(ecdf, subinterval) 
                ecdf(subinterval[2]) - ecdf(subinterval[1]),
            subinterval=subinterval,
            FUN.VALUE=numeric(1))
    }
    clusterAmbiguities <- lapply(assayECDFS, FUN=.calcPropAmbiguousClusters,
        subinterval=subinterval)
    assayOptimalK <- vapply(clusterAmbiguities, which.min, numeric(1))

    # -- Extract the class predctions for the optimal K
    assayClusterLabels <- lapply(assayClusters, FUN=.getFromClusteringResults,
        i='consensusClass')
    optimalClusterLabels <- mapply(`[[`, x=assayClusterLabels, i=assayOptimalK,
        SIMPLIFY=FALSE)
    
    # -- Assign class predictions to the experiments in tranData
    trainingData <- trainData(object)
    experimentColData <- lapply(experiments(trainingData), FUN=colData)
    updatedExpColData <- mapply(`[[<-`, x=experimentColData, i='cluster_label', 
        value=optimalClusterLabels, SIMPLIFY=FALSE)
    updatedExperiments <- mendoapply(`colData<-`, x=experiments(trainingData),
        value=updatedExpColData)
    mcols(updatedExperiments)$optimalK <- assayOptimalK + 1
    experiments(trainingData) <- updatedExperiments
    trainData(object) <- trainingData

    # -- Calculate the cluster centroids and assign to the models slot
    assayList <- assays(trainingData)
    uniqueClusterLabels <- lapply(optimalClusterLabels, unique)
    ## TODO:: Clean up this function implementation
    .calcClusterCentroid <- function(assay, uniqueLabels, labels) {
        sapply(uniqueLabels, function(i, x, labels) {
            rowMeans(x[, which(labels == i), drop=FALSE], na.rm=TRUE)
        }, x=assay, labels=labels)
    }
    clusterCentroids <- mapply(FUN=.calcClusterCentroid, assay=assayList, 
        uniqueLabels=uniqueClusterLabels, labels=optimalClusterLabels, 
        SIMPLIFY=FALSE)
    models(object) <- SimpleList(clusterCentroids)
    metadata(object)$consensuClusteringRaw <- assayClusters

    return(object)
})


# ---- NetworkCommunitySearchModel

#' Predict Metacluster Labels for a NetworkCommunitySearchModel
#' 
#' @param object A `NCSModel` which has been trained.
#' 
#' @return The `object` model with 
#' 
#' @md
#' @importFrom igraph graph_from_edgelist layout_with_fr as.undirected
#'     fastgreedy.community
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder setnames
#'   tstrsplit
#' @importFrom MultiAssayExperiment experiments experiments<-
#' @importFrom S4Vectors endoapply mendoapply merge DataFrame
#' @export
setMethod('predictClasses', signature(object='NCSModel'), function(object) {

    # -- Extract network edge data
    signifEdgeDT <- models(object)$networkEdges

    edgeMatrix <- as.matrix(signifEdgeDT[,
        .(centroid_cluster=paste0(centroid_cohort, '-', centroid_K), 
            assay_cluster=paste0(assay_cohort, '-', assay_K))
    ])

    # -- Construct graph from network edges and use the graph to predict
    #  cross-cohort meta-clusters
    graph <- graph_from_edgelist(edgeMatrix)
    coords <- layout_with_fr(graph)
    ugraph <- as.undirected(graph)
    metaclusters <- fastgreedy.community(ugraph, 
        weights=signifEdgeDT$cor_threshold)

    # -- Assign raw graph results to the NCSModel object
    models(object) <- c(SimpleList(list(graphData=list(
        graph=graph,
        coords=coords,
        ugraph=ugraph,
        metaclusters=metaclusters
    ))), models(object))

    # -- Format the metacluster results into a data.table
    metaclusterDT <- data.table(
        tmp=metaclusters$names,
        metacluster_label=metaclusters$membership
    )
    metaclusterDT[, c('cohort', 'cluster_label') := tstrsplit(tmp, '-')]
    metaclusterDT[, `:=`(tmp=NULL, cluster_label=as.integer(cluster_label))]
    setcolorder(metaclusterDT, c('cohort', 'cluster_label', 'metacluster_label'))

    models(object)$metaclusterDT <- metaclusterDT

    # -- Assign metacluster lables to each SummarizedExperiment
    cohortMAE <- trainData(object)
    colDataL <- lapply(experiments(cohortMAE), colData)
    metaclusterL <- split(metaclusterDT, by='cohort')[names(colDataL)]
    metaclusterL <- lapply(metaclusterL, function(DT) DataFrame(DT[, cohort := NULL]))
    metaclustColDataL <- mendoapply(merge, colDataL, metaclusterL,
        by='cluster_label', all.x=TRUE)
    colDataRownames <- lapply(colDataL, rownames)
    metaclustColDataL <- mendoapply(`rownames<-`, x=metaclustColDataL, 
        value=colDataRownames)
    experiments(cohortMAE) <- mendoapply(`colData<-`, x=experiments(cohortMAE),
        value=metaclustColDataL)
    trainData(object) <- cohortMAE

    # -- Find cohort clusters with no metacluster assignment
    notMetaclustered <- lapply(metaclustColDataL, 
        function(x) unique(x[is.na(x$metacluster_label), ]$cluster_label))
    notMetaclustered <- notMetaclustered[!isEmpty(notMetaclustered)]
    notMetaclusteredDT <- data.table(
        cohort=names(notMetaclustered),
        cluster_label=notMetaclustered
    )

    # -- Add not clustered information into models
    models(object)$notMetaclusteredDT <- notMetaclusteredDT

    return(object)
})