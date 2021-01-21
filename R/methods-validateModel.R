#' Perform Validation on an `S4` Object Respresenting a Trained Model
#'
#' @param model An `S4` object
#' @param valData `Any` Data to verify the model with.
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
#' @return
#'
#' @importFrom BiocParallel bplapply
#' @importFrom switchBox SWAP.KTSP.Classify
#' @import data.table
#' @md
#' @export
setMethod('validateModel', signature(model='PCOSP_or_RLS_or_RGA',
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
                SurvivalExperiments on validationData does not have model
                model predictions, recalculating...'))
            predCohortList <- predictClasses(valData, model=model)
        }
    } else {
        predCohortList <- predictClasses(valData, model=model)
    }

    # validate the model against the validation data
    valPCOSPmodelList <-
        bplapply(predCohortList, validateModel, model=model)#, ...)
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
    mcols(predCohortList)$isValidated <- TRUE
    validationData(validatedPCOSPmodel) <- predCohortList
    return(validatedPCOSPmodel)
})

#' Validate a PCOSP model with a single SurvivalExperiment object.
#'
#' @param model A `PCOSP` model which has been trained using `trainModel`.
#' @param validationData A `SurvivalExperiment` to validate the model with.
#'
#' @return The `PCOSPmodel` with the validation statistics in the `validationStats`
#'   slot and the validation data in the `validationData` slot.
#'
#' @md
#' @import S4Vectors
#' @import data.table
#' @importFrom survcomp D.index concordance.index combine.est
#' @importFrom reportROC reportROC
#' @importFrom verification roc.area
#' @export
setMethod('validateModel', signature(model='PCOSP_or_RLS_or_RGA',
    valData='SurvivalExperiment'), function(model, valData)
{
    # determine if we need to rerun the classification model
    if (identical(metadata(model)$modelParams, metadata(valData)$PCOSPparams))
    {
        survivalDF <- colData(valData)[, c('sample_name', 'days_survived',
            'is_deceased', 'PCOSP_prob_good', 'prognosis')]
        predSurvExp <- valData
    } else {
        predSurvExp <- predictClasses(model, valData)
        survivalDF <- colData(predSurvExp)[, c('sample_name', 'days_survived',
            'is_deceased', 'PCOSP_prob_good', 'prognosis')]
    }

    # convert prognosis to numeric for the ROC stats
    survivalDF <- within(survivalDF,
        prognosis <- ifelse(prognosis == 'good', 1, 0)
    )

    # calculate AUROC statistics
    aucStats <- with(survivalDF,
        c(as.numeric(reportROC(prognosis, PCOSP_prob_good,
                plot=FALSE)[c('AUC', 'AUC.SE', 'AUC.low', 'AUC.up')]),
            roc.area(prognosis, PCOSP_prob_good)$p.value, nrow(survivalDF)))
    names(aucStats) <- c('estimate', 'se', 'lower', 'upper', 'p.value', 'n')

    # calculate the validation statistcs
    validationStats <- with(survivalDF,
        list(
            dIndex=D.index(x=1 - PCOSP_prob_good, surv.time=days_survived,
                surv.event=is_deceased, na.rm=TRUE, alpha=0.5,
                method.test='logrank'),
            cIndex=concordance.index(x=1 - PCOSP_prob_good, surv.time=days_survived,
                surv.event=is_deceased, method='noether', na.rm=TRUE),
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

    validationStats(PCOSPmodel) <- valStatsDF
    validationData(PCOSPmodel) <- CohortList(list(predSurvExp),
        mDataTypes=metadata(predSurvExp)$mDataType)
    return(PCOSPmodel)
})