#' Perform Validation on an `S4` Object Representing a Trained Model
#'
#' @param model An `S4` object.
#' @param valData `Any` Data to verify the model with.
#' @param ... Allow new parameters to be defined for this generic.
#'
#' @return The `S4` object with added model performance metadata.
#'
#' @examples
#' data(sampleTrainedPCOSPmodel)
#' data(samplePCOSPpredList)
#'
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
#'
#' # Validate model
#' validatedPCOSPmodel <- validateModel(sampleTrainedPCOSPmodel,
#'   valData=samplePCOSPpredList[[1]])
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
#' data(sampleTrainedPCOSPmodel)
#' data(samplePCOSPpredList)
#'
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
#'
#' # Validate model
#' validatedPCOSPmodel <- validateModel(sampleTrainedPCOSPmodel,
#'   valData=samplePCOSPpredList)
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
    funContext <- .context(1)

    # determine if the validation data already has predictions and if
    #   if the predictions were made using the same model
    if ('hasPredictions' %in% colnames(mcols(valData))) {
        if (all(mcols(valData)$hasPredictions)) {
            if (all.equal(model, metadata(valData)$predictionModel) == TRUE) {
                predCohortList <- valData
            } else {
                warning(.warnMsg(funContext, 'The validationData argument ',
                    'has predictions, but the prediction model does not match',
                    'the model argument. Recalculating classes...'))
                predCohortList <- predictClasses(valData, model=model)
            }
        } else {
          warning(.warnMsg(funContext, 'One or more of the
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
    hasSubtypes <- 'subtype' %in% colnames(validationDT)
    rep <- if (hasSubtypes) length(unique(validationDT$subtype)) * 3 else 3
    validationDT[, `:=`(cohort=rep(names(predCohortList), each=rep),
        mDataType=rep(mcols(predCohortList)$mDataType, each=rep))]
    validationStats(validatedPCOSPmodel) <- copy(validationDT)

    by <- c('mDataType', 'statistic')
    if (hasSubtypes) by <- c(by, 'subtype')
    # calculate the per molecular data type statistics
    byMolecDT <- validationDT[,
        j=c(combine.est(estimate, se, hetero=FALSE, na.rm=TRUE), n=sum(n),
            isSummary=TRUE),
        by=by]
    byMolecDT[, cohort := toupper(mDataType)]

    by <- c('statistic')
    if (hasSubtypes) by <- c(by, 'subtype')
    overallDT <- validationDT[,
        j=c(combine.est(estimate, se, hetero=TRUE, na.rm=TRUE), n=sum(n),
                isSummary=TRUE),
            by=by]

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
        2 * pnorm(-abs(estimate / se))
    .conIndexMetaPValue <- function(estimate, se)
        2 * pnorm((estimate - 0.5) / se, lower.tail=estimate < 0.5)

    by <- c('statistic', 'mDataType')
    if (hasSubtypes) by <- c(by, 'subtype')

    combinedDT[,
        p_value := fifelse(statistic == 'log_D_index',
            .dIndexMetaPValue(estimate, se),
            .conIndexMetaPValue(estimate, se)),
        by=by]

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
#' @return
#' The `PCOSPmodel` with the validation statistics in the `validationStats`
#' slot and the validation data in the `validationData` slot.
#'
#' @examples
#' data(sampleTrainedPCOSPmodel)
#' data(samplePCOSPpredList)
#'
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
#'
#' # Validate model
#' validatedPCOSPmodel <- validateModel(sampleTrainedPCOSPmodel,
#'   valData=samplePCOSPpredList)
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
        survivalDT <- as.data.table(colData(valData))
        predSurvExp <- valData
    } else {
        predSurvExp <- predictClasses(valData, model)
        survivalDT <- as.data.table(colData(predSurvExp))
    }

    # convert prognosis to numeric for the ROC stats
    survivalDT[, prognosis := ifelse(prognosis == 'good', 1L, 0L)]

    # make the name of the score column from the class of the model
    riskProbCol <- paste0(class(model)[1], '_prob_good')

    # calcualte normal validation statistics
    ## TODO:: Can we clean this up a bit?
    valStatsDT <- .calculateUntransposedValStatsDT(survivalDT, riskProbCol)
    valStatsDT <- transpose(valStatsDT, keep.names='statistic')
    colnames(valStatsDT) <-  c('statistic', 'estimate', 'se', 'lower', 'upper',
        'p_value', 'n')
    valStatsDT$cohort <- class(model)

    # calculate subtype specific validation statistics
    if ('hasSubtypes' %in% names(metadata(valData)) && metadata(valData)$hasSubtype) {
        subtypeStatsDT <- .calculateUntransposedValStatsDT(survivalDT,
            riskProbCol, by='subtype')
        subtypeList <- list()
        for (subtypeName in unique(survivalDT$subtype)) {
            DT <- transpose(subtypeStatsDT[subtype == subtypeName, -'subtype'],
                keep.names='statistic')
            DT[, subtype := subtypeName]
            subtypeList <- c(subtypeList, list(DT))
        }
        subtypeStatsDT <- rbindlist(subtypeList)
        colnames(subtypeStatsDT) <- c('statistic', 'estimate', 'se', 'lower',
            'upper', 'p_value', 'n', 'subtype')
        subtypeStatsDT$cohort <- class(model)
        valStatsDT$subtype <- 'all'
        valStatsDT <- rbind(valStatsDT, subtypeStatsDT)
        metadata(model)$hasSubtpyes <- TRUE
    }

    valStatsDT[, isSummary := FALSE]
    validationStats(model) <- valStatsDT
    validationData(model) <- CohortList(list(predSurvExp),
        mDataTypes=metadata(predSurvExp)$mDataType)
    metadata(model)$isValidated <- TRUE
    return(model)
})

#' Calculate the Survival Validation Statistics
#'
#' @param DT A `data.table` produced from coercing the colData of a
#'   `SurvivalExperiment` object with [`as.data.table`]
#' @param riskProbCol A character vector witht he name of the column with risk
#'   probabitlies in it.
#'
#' @return A `data.table` where the columns are AUC, log_D_index,
#'   concordance_index, and any arguments specified to by in ...
#'
#' @md
#' @noRd
#' @keywords internal
.calculateUntransposedValStatsDT <- function(DT, riskProbCol, ...) {
    DT[, .(AUC=.safe_AUC(prognosis, get(riskProbCol), .N),
        log_D_index=.safe_D.index(get(riskProbCol), survival_time, event_occurred, .N),
        concordance_index=.safe_concordance.index(get(riskProbCol),
            survival_time, event_occurred, .N)
        ),
        ...
    ]
}

#' Calculate AUC or return NAs if it fails
#'
#' @noRd
#' @keywords internal
.safe_AUC <- function(prognosis, risk, n) {
    as.numeric(
        c(tryCatch({
            x <- reportROC(prognosis, risk, plot=FALSE)[c('AUC', 'AUC.SE', 'AUC.low',
                'AUC.up')]
            if (length(x) != 4) stop('report ROC failed') else x
            }, error=function(e) { print(e); return(rep(NA_real_, 4)) }),
        tryCatch({ roc.area(prognosis, risk)$p.value },
            error=function(e) { print(e); return(NA_real_) }),
        c(n=n))
        )
}

#' Calculate log D.index or return NAs if it fails
#'
#' @noRd
#' @keywords internal
.safe_D.index <- function(risk, survival_time, event_occurred, n) {
    as.numeric(
        tryCatch({
            dIndex <- D.index(x=1 - risk, surv.time=survival_time,
                surv.event=event_occurred, na.rm=TRUE, alpha=0.5,
                method.test='logrank')[c('coef', 'se', 'lower', 'upper',
                    'p.value', 'n')]
            dIndex[['lower']] <- log(dIndex[['lower']])
            dIndex[['upper']] <- log(dIndex[['upper']])
            dIndex
                },
            error=function(e) { print(e); return(c(rep(NA, 5), n))}
        )
    )
}

#' Calculate concordance.index or return NAs if it fails
#'
#' @noRd
#' @keywords internal
.safe_concordance.index <- function(risk, survival_time, event_occurred, n) {
    as.numeric(
        tryCatch({
            concordance.index(x=1 - risk,
                surv.time=survival_time, surv.event=event_occurred,
                method='noether', na.rm=TRUE)[c('c.index', 'se', 'lower',
                    'upper', 'p.value', 'n')] },
            error=function(e) {print(e); return(c(rep(NA, 5), n)) }))
}

# ---- ClinicalModel methods

#' Validate a ClinicalModel object with a single SurvivalExperiment object.
#'
#' @param model A `ClinicalModel` object which has been trained using
#'  `trainModel`.
#' @param valData A `SurvivalExperiment` to validate the model with.
#'
#' @return The `ClinicalModel` with the validation statistics in the
#'   `validationStats` slot and the validation data in the `validationData` slot.
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
        survivalDF <- colData(valData)[, c('sample_name', 'survival_time',
            'event_occurred', 'clinical_prob_good', 'prognosis')]
        predSurvExp <- valData
    } else {
        predSurvExp <- predictClasses(model, valData)
        survivalDF <- colData(predSurvExp)[, c('sample_name', 'survival_time',
            'event_occurred', 'clinical_prob_good', 'prognosis')]
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
            dIndex=D.index(x=clinical_prob_good, surv.time=survival_time,
                surv.event=event_occurred, na.rm=TRUE, alpha=0.5,
                method.test='logrank'),
            cIndex=concordance.index(x=clinical_prob_good,
                surv.time=survival_time, surv.event=event_occurred,
                method='noether', na.rm=TRUE),
            AUC=as.list(aucStats)
        )
    )

    # assemble into a data.frame
    valStatsDF <- data.frame(
        statistic=c('log_D_index', 'concordance_index', 'AUC'),
        estimate=c(validationStats$dIndex$coef,
            validationStats$cIndex$c.index, validationStats$AUC$estimate),
        se=vapply(validationStats, `[[`, i='se', FUN.VALUE=numeric(1)),
        lower=log(vapply(validationStats, `[[`, i='lower', FUN.VALUE=numeric(1))),
        upper=log(vapply(validationStats, `[[`, i='upper', FUN.VALUE=numeric(1))),
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
#' # Set parallelization settings
#' BiocParallel::register(BiocParallel::SerialParam())
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
    funContext <- .context(1)
    # determine if the validation data already has predictions and if
    #   if the predictions were made using the same model
    if ('hasPredictions' %in% colnames(mcols(valData))) {
        if (all(mcols(valData)$hasPredictions)) {
            if (all.equal(model, metadata(valData)$predictionModel)) {
                predCohortList <- valData
            } else {
                warning(.warnMsg(funContext, 'The validationData argument ',
                    'has predictions, but the prediction model does not match',
                    'the model argument. Recalculating classes...'))
                predCohortList <- predictClasses(valData, model=model)
            }
        } else {
            warning(.warnMsg(funContext, 'One or more of the
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
        2 * pnorm(-abs(estimate / se))
    .conIndexMetaPValue <- function(estimate, se)
        2 * pnorm((estimate - 0.5) / se, lower.tail=estimate < 0.5)

    combinedDT[,
        p_value := fifelse(statistic == 'log_D_index',
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
#'   `validationStats` slot and the validation data in the `validationData` slot.
#'
#' @md
#' @importFrom data.table data.table as.data.table merge.data.table rbindlist
#'   `:=` copy .N .SD fifelse merge.data.table transpose setcolorder
#' @import survcomp
#' @importFrom CoreGx .warnMsg
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @importFrom reportROC reportROC
#' @importFrom verification roc.area
#' @importFrom survival strata
#' @export
setMethod('validateModel', signature(model='GeneFuModel',
    valData='SurvivalExperiment'), function(model, valData)
{
    funContext <- .context(1)

    survivalDF <- colData(valData)
    predSurvExp <- valData

    # deal with missing prognosis column
    if (!('prognosis' %in% colnames(survivalDF))) {
        warning(.warnMsg(funContext, 'The prognosis column is missing from',
            ' the validation SurvivalExperiment, calculating based on ',
            'minDaysSurvived value in modelParams...'))
        survivalDF <- within(survivalDF,
            prognosis <- ifelse(survival_time >
                metadata(model)$modelParams$minDaysSurvived, 'good', 'bad'))
    }

    # convert prognosis to numeric for the ROC stats
    survivalDF <- within(survivalDF, {
        prognosis <- ifelse(prognosis == 'good', 1L, 0L)
    })

    # calculate the validation statistics
    validationStats <- with(survivalDF,
        list(
            dIndex=D.index(x=genefu_score, surv.time=survival_time,
                surv.event=event_occurred, na.rm=TRUE, alpha=0.5,
                method.test='logrank'),
            cIndex=concordance.index(x=genefu_score, surv.time=survival_time,
                surv.event=event_occurred, method='noether', na.rm=TRUE)
        )
    )

    # assemble into a data.frame
    valStatsDF <- data.frame(
        statistic=c('log_D_index', 'concordance_index'),
        estimate=c(validationStats$dIndex$coef,
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
    funContext <- .context(1)
    # determine if the validation data already has predictions and if
    #   if the predictions were made using the same model
    if ('hasPredictions' %in% colnames(mcols(valData))) {
        if (all(mcols(valData)$hasPredictions)) {
            if (all.equal(model, metadata(valData)$predictionModel)) {
                predCohortList <- valData
            } else {
                warning(.warnMsg(funContext, 'The validationData argument ',
                    'has predictions, but the prediction model does not match',
                    'the model argument. Recalculating classes...'))
                predCohortList <- predictClasses(valData, model=model)
            }
        } else {
          warning(.warnMsg(funContext, 'One or more of the
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
        2 * pnorm(-abs(estimate / se))
    .conIndexMetaPValue <- function(estimate, se)
        2 * pnorm((estimate - 0.5) / se, lower.tail=estimate < 0.5)

    combinedDT[,
        p_value := fifelse(statistic == 'log_D_index',
            .dIndexMetaPValue(estimate, se),
            .conIndexMetaPValue(estimate, se)),
        by=.(statistic, mDataType)]

    allValStatsDT <- rbindlist(list(validationDT, combinedDT), fill=TRUE)

    validationStats(validatedModel) <- allValStatsDT
    metadata(validatedModel)$isValidated <- TRUE
    validationData(validatedModel) <- predCohortList
    return(validatedModel)
})


# ---- ConsensusMetaclusteringModel

#' 
#' @param object A `ConsensusMetaclusteringModel` object with cluster_labels in 
#'   assigned to each experiment, as returned by `predictClasses`.
#' @param valData A `ConsensusMetaclusteringModel` object with cluster_labels
#'   assigned to each experiment, as returned by `predictClasses`. This
#'   consensus cluster should contain outgroup cohorts, such as normal 
#'   patients to be compared against the disease cohorts being used for 
#'   class discovery.
#' @param ... Fallthrough parameters to `BiocParallel::bpmapply`. This can
#'   also be used to customize the call to `stats::cor.test` used for 
#'   calculating the cluster thresholds.
#' 
#' @return The `ConsensusMetaclusteringModel` from object, with the training
#'   data from `valData` in the `validationData` slot, the models from the
#'   `valData` object appended to the `models` of object, and the 
#'   `validationStats` slot populated with pair-wise comparisons between
#'   all experiments in both `object` and `valData`.
#' 
#' @importFrom BiocParallel bpmapply
#' @importFrom data.table rbindlist
#' 
#' @md
#' @export
setMethod('validateModel', signature(model='ConsensusMetaclusteringModel', 
    valData='ConsensusMetaclusteringModel'), function(model, valData, ...) 
{
    funContext <- .context(1)
    
    # Prepare comparisons for clustering results between every cohort
    validationData(model) <- SimpleList(experiments=trainData(valData),
        models=models(valData))
    valCohorts <- names(trainData(valData))

    repeats <- modelParams(model)$reps
    centroids <- c(models(model), models(valData))
    cohorts <- c(experiments(trainData(model)), 
        experiments(trainData(valData)))
    assays <- lapply(cohorts, assay, 1)
    clusterLabels <- lapply(cohorts, function(x) x$cluster_label)

    # Find all non-self comparisons between the cohorts
    cohortPairs <- .findAllCohortPairs(names(cohorts))
    metadata(model)$cohortPairs <- cohortPairs

    # Extract and expand the data for each cohort
    comparisonCentroids <- centroids[cohortPairs$x]
    comparisonAssays <- assays[cohortPairs$y]
    comparisonClasses <-  clusterLabels[cohortPairs$y]
    comparisons <- rownames(cohortPairs)
    centroidOptimalK <- mcols(cohorts)$optimalK[cohortPairs$x]
    assayOptimalK <- mcols(cohorts)$optimalK[cohortPairs$y]

    # Use cor.test to estimate association between the feature values
    #  in a centroid vs a cohort. Compare each subcluster individually for
    #  all assays.
    thresholdDtList <- bpmapply(FUN=.calculateMSMthresholds, 
        comparison=comparisons, centroid=comparisonCentroids, 
        assay=comparisonAssays, classes=comparisonClasses, 
        centroidK=centroidOptimalK, assayK=assayOptimalK, SIMPLIFY=FALSE)

    thresholdDT <- rbindlist(thresholdDtList)
    modelParams(model) <- c(modelParams(model), 
        list(thresholdDT=thresholdDT))

    clusterReproDtList <- bpmapply(FUN=.calcClusterRepro,
        comparison=comparisons, centroid=comparisonCentroids, 
        assay=comparisonAssays, MoreArgs=list(rep=repeats), 
        SIMPLIFY=FALSE)

    reproDT <- rbindlist(clusterReproDtList)

    validationStats(model) <- rbind(reproDT, thresholdDT, fill=TRUE)

    return(model)
})

#' @importFrom clusterRepro clusterRepro
#' @importFrom data.table data.table
#' 
#' @md
#' @keywords internal
.calcClusterRepro <- function(comparison, centroid, assay, rep)
{
    message(comparison)
    reproStats <- as.data.table(clusterRepro(centroid, assay, rep))
    cohorts <- unlist(strsplit(comparison, '-'))
    reproDT <- data.table(
        'metric'='clusterRepro',
        'comparison'=comparison,
        'centroid_cohort'=cohorts[1],
        'assay_cohort'=cohorts[2],
        'centroid_K'=seq_len(ncol(centroid)),
        'estimate'=reproStats$Actual.IGP,
        'p_value'=reproStats$p.value,
        'assay_N'=reproStats$Actual.Size
    )
    return(reproDT)
}

#' @importFrom data.table data.table rbindlist
#' 
#' @param ... Fallthrough arguments to `stats::cor.test`
#' 
#' @md
#' @keywords internal
.calculateMSMthresholds <- function(comparison, centroid, assay, classes, 
    centroidK, assayK, ...)
{
    cohorts <- unlist(strsplit(comparison, '-'))
    sharedGenes <- intersect(rownames(centroid), rownames(assay))
    centroid <- na.omit(centroid[sharedGenes, , drop=FALSE])
    assay <- na.omit(assay[sharedGenes, , drop=FALSE])
    corList <- vector("list", centroidK * assayK)
    k <- 1
    for (i in seq_len(centroidK)) {
      for (j in seq_len(assayK)) {
        classIdxs <- which(classes == j)
        if (length(classIdxs) < 1) stop(.errorMsg(funContext, 'No samples have',
            'the class ', j, ' in ', cohorts[2], '. Something has gone wrong!'))
        localCentroid <- centroid[, i]
        corTestResults <- lapply(classIdxs,
                              function(idx, localCentroid, data) {
                                  cor.test(
                                    localCentroid,
                                    data[, idx],
                                    ...)
                                  },
                              localCentroid=localCentroid,
                              data=assay)
        estimate <- mean(vapply(corTestResults, function(x) x$estimate, numeric(1)), na.rm=TRUE)
        p_value <- mean(vapply(corTestResults, function(x) x$p.value, numeric(1)), na.rm=TRUE)
        corList[[k]] <- data.table(
            'metric'='cor.test',
            'comparison'=comparison,
            'centroid_cohort'=cohorts[1],
            'assay_cohort'=cohorts[2],
            'centroid_K'=i,
            'assay_K'=j,
            'estimate'=estimate,
            'p_value'=p_value)
        k <- k + 1
      }
    }
    rbindlist(corList, fill=TRUE)
}

#' Find all non-self pair-wise combinations of cohorts
#'
#' @param clusterNames A `character` vector of cohort names
#'
#' @return A `data.frame` with the index of all non-self pair-wise cohort combinations. Rownames
#'     are the names of the two clusters being compared.
#'
#' @md
#' @keywords internal
.findAllCohortPairs <- function(clusterNames) {
    pairs <- expand.grid(x=seq_along(clusterNames), y=seq_along(clusterNames))
    namePairs <- expand.grid(x=clusterNames, y=clusterNames)

    # Remove self comparisons
    allPairs <- pairs[-which(pairs[, 1] == pairs[, 2]), ]
    allNames <- namePairs[-which(namePairs[, 1] == namePairs[, 2]), ]

    # Paste together pair names
    pairNames <- mapply(paste, allNames[, 1], allNames[, 2], MoreArgs=list(sep="-"))

    # Assign names as rownames to pair data.frame
    rownames(allPairs) <- pairNames

    return(allPairs)
}