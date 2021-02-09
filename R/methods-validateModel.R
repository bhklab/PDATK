#' Perform Validation on an `S4` Object Respresenting a Trained Model
#'
#' @param model An `S4` object.
#' @param valData `Any` Data to verify the model with.
#' @param ... Allow new parameters to be defined for this generic.
#'
#' @return The `S4` object with added model performance metadata.
#'
#' @examples
#' data(samplePCOSPmodel)
#' data(sampleCohortList)
#'
#' # Train Model
#' trainedPCOSPmodel <- trainModel(samplePCOSPmodel, numModels=10,
#'   minAccuracy=0.6)
#'
#' # Make predictions
#' PCOSPpredCohortList <- predictClasses(sampleCohortList,
#'   model=trainedPCOSPmodel)
#'
#' # Validate model
#' validatedPCOSPmodel <- validateModel(trainedPCOSPmodel,
#'   valData=PCOSPpredCohortList)
#'
#' @md
#' @export
setGeneric('validateModel', function(model, valData, ...)
    standardGeneric('validateModel'))
#'
#' Evaluate the Performance of a List of Trained KTSP Models from a PCOSP
#'   Model
#'
#' @param model A `PCOSP` model which has been trained using `trainModel`.
#' @param valData A `CohortList` containing one or more
#'   `SurvivalExperiment`s. The first assay in each `SurvivalExperiment` will
#'   be classified using all top scoring KTSP models in `models(model)`.
#' @param ... Fallthrough arguments to `BiocParallel::bplapply`, use this to
#'   configure the parallelization settings for this function. For example
#'   to specify BPARAM.
#'
#' @seealso [`BiocParallel::bplapply`], [`switchBox::SWAP.KTSP.Classify`]
#'
#' @return The `model` object with the validationStats and validationData
#'   slots occupied.
#'
#' @examples
#' data(samplePCOSPmodel)
#' data(sampleCohortList)
#'
#' # Train Model
#' trainedPCOSPmodel <- trainModel(samplePCOSPmodel, numModels=10,
#'   minAccuracy=0.6)
#'
#' # Make predictions
#' PCOSPpredCohortList <- predictClasses(sampleCohortList,
#'   model=trainedPCOSPmodel)
#'
#' # Validate model
#' validatedPCOSPmodel <- validateModel(trainedPCOSPmodel,
#'   valData=PCOSPpredCohortList)
#'
#' @md
#' @import survcomp
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder
#' @importFrom BiocParallel bplapply
#' @importFrom switchBox SWAP.KTSP.Classify
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom stats pnorm qnorm
#' @importFrom survival strata
#' @importFrom S4Vectors metadata mcols
#' @export
setMethod('validateModel', signature(model='PCOSP_or_RLS_or_RGA',
    valData='CohortList'), function(model, valData, ...)
{
    # determine if the validation data already has predictions and if
    #   if the predictions were made using the same model
    if ('hasPredictions' %in% colnames(mcols(valData))) {
        if (all(mcols(valData)$hasPredictions)) {
            if (all.equal(model, metadata(valData)$predictionModel) == TRUE) {
                predCohortList <- valData
            } else {
                warning(.warnMsg(.context(6), 'The validationData argument ',
                    'has predictions, but the prediction model does not match',
                    'the model argument. Recalculating classes...'))
                predCohortList <- predictClasses(valData, model=model)
            }
        } else {
          warning(.warnMsg(.context(6), 'One or more of the
                SurvivalExperiments in valData does not have model
                model predictions, recalculating...'))
            predCohortList <- predictClasses(valData, model=model)
        }
    } else {
        predCohortList <- predictClasses(valData, model=model)
    }

    # validate the model against the validation data
    valPCOSPmodelList <-
        bplapply(predCohortList, validateModel, model=model, ...)
    validatedPCOSPmodel <- valPCOSPmodelList[[1]]
    validationDT <- rbindlist(lapply(valPCOSPmodelList, validationStats))
    validationDT[, `:=`(cohort=rep(names(predCohortList), each=3),
        mDataType=rep(mcols(predCohortList)$mDataType, each=3))]
    validationStats(validatedPCOSPmodel) <- copy(validationDT)

    # calculate the per molecular data type statistics
    byMolecDT <- validationDT[,
        j=c(combine.est(estimate, se, hetero=FALSE, na.rm=TRUE), n=sum(n),
            isSummary=TRUE),
        by=.(mDataType, statistic)]
    byMolecDT[, cohort := toupper(mDataType)]

    overallDT <- validationDT[,
        j=c(combine.est(estimate, se, hetero=TRUE, na.rm=TRUE), n=sum(n),
            isSummary=TRUE),
        by=statistic]

    # calculate the overall statistics
    overallDT[, `:=`(mDataType='combined', cohort='OVERALL')]
    setcolorder(overallDT, colnames(byMolecDT))

    # get lower, upper, p-value and n for the aggregate statistics
    combinedDT <- rbindlist(list(byMolecDT, overallDT))
    combinedDT[, lower := estimate + qnorm(0.025, lower.tail=TRUE) * se,
        by=.(statistic, mDataType)]
    combinedDT[, upper := estimate + qnorm(0.025, lower.tail=FALSE) * se,
        by=.(statistic, mDataType)]

    .dIndexMetaPValue <- function(estimate, se)
        2 * pnorm(-abs(log(estimate) / se))
    .conIndexMetaPValue <- function(estimate, se)
        2 * pnorm((estimate - 0.5) / se, lower.tail=estimate < 0.5)

    combinedDT[,
        p_value := fifelse(statistic == 'D_index',
            .dIndexMetaPValue(estimate, se),
            .conIndexMetaPValue(estimate, se)),
        by=.(statistic, mDataType)]

    allValStatsDT <- rbindlist(list(validationDT, combinedDT), fill=TRUE)

    validationStats(validatedPCOSPmodel) <- allValStatsDT
    metadata(validatedPCOSPmodel)$isValidated <- TRUE
    validationData(validatedPCOSPmodel) <- predCohortList
    return(validatedPCOSPmodel)
})
#'
#' Validate a PCOSP model with a single SurvivalExperiment object.
#'
#' @param model A `PCOSP` model which has been trained using `trainModel`.
#' @param valData A `SurvivalExperiment` to validate the model with.
#'
#' @return The `PCOSPmodel` with the validation statistics in the `validationStats`
#'   slot and the validation data in the `validationData` slot.
#'
#' @examples
#' data(samplePCOSPmodel)
#' data(samplePCSIsurvExp)
#'
#' # Train Model
#' trainedPCOSPmodel <- trainModel(samplePCOSPmodel, numModels=10,
#'   minAccuracy=0.6)
#'
#' # Make predictions
#' PCOSPpredSurvExp <- predictClasses(samplePCSIsurvExp,
#'   model=trainedPCOSPmodel)
#'
#' # Validate model
#' validatedPCOSPmodel <- validateModel(trainedPCOSPmodel,
#'   valData=PCOSPpredSurvExp)
#'
#' @md
#' @include classUnions.R
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder
#' @import survcomp
#' @import reportROC
#' @import survival
#' @import verification
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @export
setMethod('validateModel', signature(model='PCOSP_or_RLS_or_RGA',
    valData='SurvivalExperiment'), function(model, valData)
{

    # determine if we need to rerun the classification model
    if (identical(metadata(model)$modelParams,
        metadata(valData)[[paste0(class(model), 'params')]]))
    {
        survivalDF <- colData(valData)
        predSurvExp <- valData
    } else {
        predSurvExp <- predictClasses(valData, model)
        survivalDF <- colData(predSurvExp)
    }

    # convert prognosis to numeric for the ROC stats
    survivalDF <- within(survivalDF,
        prognosis <- ifelse(prognosis == 'good', 1L, 0L)
    )

    # calculate AUROC statistics
    riskProbCol <- paste0(class(model)[1], '_prob_good')
    aucStats <- with(survivalDF,
        c(as.numeric(reportROC(prognosis, get(riskProbCol),
                plot=FALSE)[c('AUC', 'AUC.SE', 'AUC.low', 'AUC.up')]),
            roc.area(prognosis, get(riskProbCol))$p.value, nrow(survivalDF)))
    names(aucStats) <- c('estimate', 'se', 'lower', 'upper', 'p.value', 'n')

    # calculate the validation statistcs
    validationStats <- with(survivalDF,
        list(
            dIndex=D.index(x=1 - get(riskProbCol), surv.time=days_survived,
                surv.event=is_deceased, na.rm=TRUE, alpha=0.5,
                method.test='logrank'),
            cIndex=concordance.index(x=1 - get(riskProbCol),
                surv.time=days_survived, surv.event=is_deceased,
                method='noether', na.rm=TRUE),
            AUC=as.list(aucStats)
        )
    )

    # assemble into a data.frame
    valStatsDF <- data.frame(
        statistic=c('D_index', 'concordance_index', 'AUC'),
        estimate=c(validationStats$dIndex$d.index,
            validationStats$cIndex$c.index, validationStats$AUC$estimate),
        se=vapply(validationStats, `[[`, i='se', FUN.VALUE=numeric(1)),
        lower=vapply(validationStats, `[[`, i='lower', FUN.VALUE=numeric(1)),
        upper=vapply(validationStats, `[[`, i='upper', FUN.VALUE=numeric(1)),
        p_value=vapply(validationStats, `[[`, i='p.value', FUN.VALUE=numeric(1)),
        n=vapply(validationStats, `[[`, i='n', FUN.VALUE=numeric(1)),
        isSummary=FALSE
    )

    valStatsDF$cohort <- class(model)

    validationStats(model) <- valStatsDF
    validationData(model) <- CohortList(list(predSurvExp),
        mDataTypes=metadata(predSurvExp)$mDataType)
    metadata(model)$isValidated <- TRUE
    return(model)
})

# ---- ClinicalModel methods

#' Validate a ClinicalModel object with a single SurvivalExperiment object.
#'
#' @param model A `ClinicalModel` object which has been trained using
#'  `trainModel`.
#' @param valData A `SurvivalExperiment` to validate the model with.
#'
#' @return The `ClinicalModel` with the validation statistics in the
#'   `validationStats` slot and the validation data in the
#'   `validationData` slot.
#'
#' @examples
#' data(sampleClinicalModel)
#' data(sampleCohortList)
#'
#' # Train Model
#' trainedClinicalModel <- trainModel(sampleClinicalModel)
#'
#' # Make predictions
#' clinicalPredCohortList <- predictClasses(sampleCohortList[c('PCSI', 'TCGA')],
#'   model=trainedClinicalModel)
#'
#' # Validate model
#' validatedClinicalModel <- validateModel(trainedClinicalModel,
#'   valData=clinicalPredCohortList)
#'
#' @md
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder
#' @import survcomp
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @importFrom reportROC reportROC
#' @importFrom survival strata
#' @importFrom verification roc.area
#' @export
setMethod('validateModel', signature(model='ClinicalModel',
    valData='SurvivalExperiment'), function(model, valData)
{
    # determine if we need to rerun the classification model
    if (identical(metadata(model)$modelParams, metadata(valData)$GLMparams))
    {
        survivalDF <- colData(valData)[, c('sample_name', 'days_survived',
            'is_deceased', 'clinical_prob_good', 'prognosis')]
        predSurvExp <- valData
    } else {
        predSurvExp <- predictClasses(model, valData)
        survivalDF <- colData(predSurvExp)[, c('sample_name', 'days_survived',
            'is_deceased', 'clinical_prob_good', 'prognosis')]
    }

    # convert prognosis to numeric for the ROC stats
    survivalDF <- within(survivalDF,
        prognosis <- ifelse(prognosis == 'good', 1L, 0L)
    )

    # calculate AUROC statistics
    aucStats <- with(survivalDF,
        c(as.numeric(reportROC(prognosis, clinical_prob_good,
                plot=FALSE)[c('AUC', 'AUC.SE', 'AUC.low', 'AUC.up')]),
          roc.area(prognosis, clinical_prob_good)$p.value, nrow(survivalDF)))
    names(aucStats) <- c('estimate', 'se', 'lower', 'upper', 'p.value', 'n')

    # calculate the validation statistcs
    ## FIXME:: Should this be 1 - clinical_prob_good?
    validationStats <- with(survivalDF,
        list(
            dIndex=D.index(x=clinical_prob_good, surv.time=days_survived,
                surv.event=is_deceased, na.rm=TRUE, alpha=0.5,
                method.test='logrank'),
            cIndex=concordance.index(x=clinical_prob_good,
                surv.time=days_survived, surv.event=is_deceased,
                method='noether', na.rm=TRUE),
            AUC=as.list(aucStats)
        )
    )

    # assemble into a data.frame
    valStatsDF <- data.frame(
        statistic=c('D_index', 'concordance_index', 'AUC'),
        estimate=c(validationStats$dIndex$d.index,
            validationStats$cIndex$c.index, validationStats$AUC$estimate),
        se=vapply(validationStats, `[[`, i='se', FUN.VALUE=numeric(1)),
        lower=vapply(validationStats, `[[`, i='lower', FUN.VALUE=numeric(1)),
        upper=vapply(validationStats, `[[`, i='upper', FUN.VALUE=numeric(1)),
        p_value=vapply(validationStats, `[[`, i='p.value', FUN.VALUE=numeric(1)),
        n=vapply(validationStats, `[[`, i='n', FUN.VALUE=numeric(1)),
        isSummary=FALSE
    )

    valStatsDF$cohort <- class(model)

    validationStats(model) <- valStatsDF
    validationData(model) <- CohortList(list(predSurvExp),
        mDataTypes=metadata(predSurvExp)$mDataType)
    return(model)
})
## TODO:: Refactor this into a helper method or extend the class union to include ClinicalModel
#' @inherit validateModel,PCOSP_or_RLS_or_RGA,CohortList-method
#'
#' @param model A trained `ClinicalModel` object, as returned by the `trainModel`
#'   method.
#'
#' @examples
#' data(sampleClinicalModel)
#' data(samplePCSIsurvExp)
#'
#' # Train Model
#' trainedClinicalModel <- trainModel(sampleClinicalModel)
#'
#' # Make predictions
#' clinicalPredSurvExp <- predictClasses(samplePCSIsurvExp,
#'   model=trainedClinicalModel)
#'
#' # Validate model
#' validatedClincalModel <- validateModel(trainedClinicalModel,
#'   valData=clinicalPredSurvExp)
#'
#' @md
#' @importFrom survival strata
#' @export
setMethod('validateModel', signature(model='ClinicalModel',
    valData='CohortList'), function(model, valData, ...)
{
    # determine if the validation data already has predictions and if
    #   if the predictions were made using the same model
    if ('hasPredictions' %in% colnames(mcols(valData))) {
        if (all(mcols(valData)$hasPredictions)) {
            if (all.equal(model, metadata(valData)$predictionModel)) {
                predCohortList <- valData
            } else {
                warning(.warnMsg(.context(), 'The validationData argument ',
                    'has predictions, but the prediction model does not match',
                    'the model argument. Recalculating classes...'))
                predCohortList <- predictClasses(valData, model=model)
            }
        } else {
          warning(.warnMsg(.context(), 'One or more of the
                SurvivalExperiments in valData does not have model
                model predictions, recalculating...'))
            predCohortList <- predictClasses(valData, model=model)
        }
    } else {
        predCohortList <- predictClasses(valData, model=model)
    }

    # validate the model against the validation data
    validatedModelList <-
        bplapply(predCohortList, validateModel, model=model, ...)
    validatedModel <- validatedModelList[[1]]
    validationDT <- rbindlist(lapply(validatedModelList, validationStats))
    validationDT[, `:=`(cohort=rep(names(predCohortList), each=3),
        mDataType=rep(mcols(predCohortList)$mDataType, each=3))]
    validationStats(validatedModel) <- copy(validationDT)

    # calculate the per molecular data type statistics
    byMolecDT <- validationDT[,
        j=c(combine.est(estimate, se, hetero=FALSE, na.rm=TRUE), n=sum(n),
            isSummary=TRUE),
        by=.(mDataType, statistic)]
    byMolecDT[, cohort := toupper(mDataType)]

    overallDT <- validationDT[,
        j=c(combine.est(estimate, se, hetero=TRUE, na.rm=TRUE), n=sum(n),
            isSummary=TRUE),
        by=statistic]

    # calculate the overall statistics
    overallDT[, `:=`(mDataType='combined', cohort='OVERALL')]
    setcolorder(overallDT, colnames(byMolecDT))

    # get lower, upper, p-value and n for the aggregate statistics
    combinedDT <- rbindlist(list(byMolecDT, overallDT))
    combinedDT[, lower := estimate + qnorm(0.025, lower.tail=TRUE) * se,
        by=.(statistic, mDataType)]
    combinedDT[, upper := estimate + qnorm(0.025, lower.tail=FALSE) * se,
        by=.(statistic, mDataType)]

    .dIndexMetaPValue <- function(estimate, se)
        2 * pnorm(-abs(log(estimate) / se))
    .conIndexMetaPValue <- function(estimate, se)
        2 * pnorm((estimate - 0.5) / se, lower.tail=estimate < 0.5)

    combinedDT[,
        p_value := fifelse(statistic == 'D_index',
            .dIndexMetaPValue(estimate, se),
            .conIndexMetaPValue(estimate, se)),
        by=.(statistic, mDataType)]

    allValStatsDT <- rbindlist(list(validationDT, combinedDT), fill=TRUE)

    validationStats(validatedModel) <- allValStatsDT
    metadata(validatedModel)$isValidated <- TRUE
    validationData(validatedModel) <- predCohortList
    return(validatedModel)
})

# ---- GeneFuModel Methods

#' Validate a `GenefuModel` object with a single `SurvivalExperiment` object.
#'
#' @param model A `GenefuModel` object which has been trained using
#'  `trainModel`.
#' @param valData A `SurvivalExperiment` to validate the model with.
#'
#' @return The `GeneModel` with the validation statistics in the
#'   `validationStats` slot and the validation data in the
#'   `validationData` slot.
#'
#' @md
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder
#' @import survcomp
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @importFrom reportROC reportROC
#' @importFrom verification roc.area
#' @importFrom survival strata
#' @export
setMethod('validateModel', signature(model='GeneFuModel',
    valData='SurvivalExperiment'), function(model, valData)
{
    survivalDF <- colData(valData)
    predSurvExp <- valData

    # deal with missing prognosis column
    if (!('prognosis' %in% colnames(survivalDF))) {
        warning(.warnMsg(.context(), 'The prognosis column is missing from',
            'the validation SurvivalExperiment, calculating based on ',
            'minDaysSurvived value in modelParams...'))
        survivalDF <- within(survivalDF,
            prognosis <- ifelse(days_survived >
                metadata(model)$modelParams$minDaysSurvived, 'good', 'bad'))
    }

    # convert prognosis to numeric for the ROC stats
    survivalDF <- within(survivalDF, {
        prognosis <- ifelse(prognosis == 'good', 1L, 0L)
    })

    # calculate the validation statistics
    validationStats <- with(survivalDF,
        list(
            dIndex=D.index(x=genefu_score, surv.time=days_survived,
                surv.event=is_deceased, na.rm=TRUE, alpha=0.5,
                method.test='logrank'),
            cIndex=concordance.index(x=genefu_score, surv.time=days_survived,
                surv.event=is_deceased, method='noether', na.rm=TRUE)
        )
    )

    # assemble into a data.frame
    valStatsDF <- data.frame(
        statistic=c('D_index', 'concordance_index'),
        estimate=c(validationStats$dIndex$d.index,
            validationStats$cIndex$c.index),
        se=vapply(validationStats, `[[`, i='se', FUN.VALUE=numeric(1)),
        lower=vapply(validationStats, `[[`, i='lower', FUN.VALUE=numeric(1)),
        upper=vapply(validationStats, `[[`, i='upper', FUN.VALUE=numeric(1)),
        p_value=vapply(validationStats, `[[`, i='p.value', FUN.VALUE=numeric(1)),
        n=vapply(validationStats, `[[`, i='n', FUN.VALUE=numeric(1)),
        isSummary=FALSE
    )

    valStatsDF$cohort <- class(model)

    validationStats(model) <- valStatsDF
    validationData(model) <- CohortList(list(predSurvExp),
        mDataTypes=metadata(predSurvExp)$mDataType)
    return(model)
})


## TODO:: Refactor this into a helper method or extend the class union to include ClinicalModel
#'
#' @inherit validateModel,PCOSP_or_RLS_or_RGA,CohortList-method
#' @param model A `GeneFuModel` with a `DataFrame` of gene coefficients in
#'   the models slot.
#'
#' @md
#' @export
setMethod('validateModel', signature(model='GeneFuModel',
    valData='CohortList'), function(model, valData, ...)
{
    # determine if the validation data already has predictions and if
    #   if the predictions were made using the same model
    if ('hasPredictions' %in% colnames(mcols(valData))) {
        if (all(mcols(valData)$hasPredictions)) {
            if (all.equal(model, metadata(valData)$predictionModel)) {
                predCohortList <- valData
            } else {
                warning(.warnMsg(.context(), 'The validationData argument ',
                    'has predictions, but the prediction model does not match',
                    'the model argument. Recalculating classes...'))
                predCohortList <- predictClasses(valData, model=model)
            }
        } else {
          warning(.warnMsg(.context(), 'One or more of the
                SurvivalExperiments in valData does not have model
                model predictions, recalculating...'))
            predCohortList <- predictClasses(valData, model=model)
        }
    } else {
        predCohortList <- predictClasses(valData, model=model)
    }

    # validate the model against the validation data
    ## FIXME:: remove comment from dots
    validatedModelList <-
        bplapply(predCohortList, validateModel, model=model)#, ...)
    validatedModel <- validatedModelList[[1]]
    validationDT <- rbindlist(lapply(validatedModelList, validationStats))
    validationDT[, `:=`(
        cohort=rep(names(predCohortList), each=length(unique(statistic))),
        mDataType=rep(mcols(predCohortList)$mDataType,
            each=length(unique(statistic)))
    )]
    validationStats(validatedModel) <- copy(validationDT)

    # calculate the per molecular data type statistics
    byMolecDT <- validationDT[,
        j=c(combine.est(estimate, se, hetero=FALSE, na.rm=TRUE), n=sum(n),
            isSummary=TRUE),
        by=.(mDataType, statistic)]
    byMolecDT[, cohort := toupper(mDataType)]

    overallDT <- validationDT[,
        j=c(combine.est(estimate, se, hetero=TRUE, na.rm=TRUE), n=sum(n),
            isSummary=TRUE),
        by=statistic]

    # calculate the overall statistics
    overallDT[, `:=`(mDataType='combined', cohort='OVERALL')]
    setcolorder(overallDT, colnames(byMolecDT))

    # get lower, upper, p-value and n for the aggregate statistics
    combinedDT <- rbindlist(list(byMolecDT, overallDT))
    combinedDT[, lower := estimate + qnorm(0.025, lower.tail=TRUE) * se,
        by=.(statistic, mDataType)]
    combinedDT[, upper := estimate + qnorm(0.025, lower.tail=FALSE) * se,
        by=.(statistic, mDataType)]

    .dIndexMetaPValue <- function(estimate, se)
        2 * pnorm(-abs(log(estimate) / se))
    .conIndexMetaPValue <- function(estimate, se)
        2 * pnorm((estimate - 0.5) / se, lower.tail=estimate < 0.5)

    combinedDT[,
        p_value := fifelse(statistic == 'D_index',
            .dIndexMetaPValue(estimate, se),
            .conIndexMetaPValue(estimate, se)),
        by=.(statistic, mDataType)]

    allValStatsDT <- rbindlist(list(validationDT, combinedDT), fill=TRUE)

    validationStats(validatedModel) <- allValStatsDT
    mcols(predCohortList)$isValidated <- TRUE
    validationData(validatedModel) <- predCohortList
    return(validatedModel)
})