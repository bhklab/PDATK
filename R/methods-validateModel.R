#' Perform Validation on an `S4` Object Respresenting a Trained Model
#'
#' @param model An `S4` object
#'
#' @md
#' @export
setGeneric('validateModel', function(model, validationData, ...)
    standardGeneric('validateModel'))
#'
#' Evaluate the Performance of a List of Trained KTSP Models from a PCOSP
#'   Model
#'
#' @param model A `PCOSP` model which has been trained using `trainModel`.
#' @param validationData A `CohortList` containing one or more
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
setMethod('validateModel', signature(model='PCOSP',
    validationData='CohortList'), function(model, validationData, ...)
{
    # determine if the validation data already has predictions and if
    #   if the predictions were made using the same model
    if ('hasPredictions' %in% colnames(mcols(validationData))) {
        if (all(mcols(validationData)$hasPredictions)) {
            if (all.equal(model, metadata(validationData)$predictionModel)) {
                predCohortList <- validationData
            } else {
                warning(.warnMsg(.context(), 'The validationData argument ',
                    'has predictions, but the prediction model does not match',
                    'the model argument. Recalculating classes...'))
                predCohortList <- predictClasses(validationData, model=model)
            }
        } else {

          warning(.warnMsg(.context(), 'One or more of the
                SurvivalExperiments on validationData does not have model
                model predictions, recalculating...'))
            predCohortList <- predictClasses(validationData, model=model)
        }
    } else {
        predCohortList <- predictClasses(validationData, model=model)
    }

    # validate the model against the validation data
    valPCOSPmodelList <-
        bplapply(predCohortList, validateModel, model=model,...)
    validatedPCOSPmodel <- valPCOSPmodelList[[1]]
    validationDT <- rbindlist(lapply(valPCOSPmodelList, validationStats))
    validationDT[, `:=`(cohort=rep(names(predCohortList), each=2),
        mDataType=rep(mcols(predCohortList)$mDataType, each=2))]
    validationStats(validatedPCOSPmodel) <- copy(validationDT)

    valDataList <- Reduce(c, lapply(valPCOSPmodelList, validationData))
    names(valDataList) <- names(predCohortList)
    validationData(validatedPCOSPmodel) <- valDataList

    # calculate the per molecular data type statistics
    byMolecDT <- validationDT[,
        j=c(combine.est(estimate, se, hetero=FALSE, na.rm=TRUE), n=sum(n)),
        by=.(mDataType, statistic)]
    byMolecDT[, cohort := toupper(mDataType)]

    overallDT <- validationDT[,
        j=c(combine.est(estimate, se, hetero=TRUE, na.rm=TRUE), n=sum(n)),
        by=statistic]

    # calculate the overall statistics
    overallDT[, `:=`(mDataType='combined', cohort='ALL')]
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
    validationData(validatedPCOSPmodel) <- validationCohortList
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
#' @importFrom survcomp D.index concordance.index combine.est
#' @import S4Vectors
#' @export
setMethod('validateModel', signature(model='PCOSP',
    validationData='SurvivalExperiment'), function(model, validationData)
{
    # determine if we need to rerun the classification model
    if (identical(metadata(model)$modelParams, metadata(validationData)$PCOSPparams))
    {
        survivalDF <- colData(validationData)[, c('sample_name', 'days_survived',
            'is_deceased', 'PCOSP_prob_good')]
        predSurvExp <- validationData
    } else {
        predSurvExp <- predictClasses(model, validationData)
        survivalDF <- colData(predSurvExp)[, c('sample_name', 'days_survived',
            'is_deceased', 'PCOSP_prob_good')]
    }

    # calculate the validation statistcs
    validationStats <- with(survivalDF,
        list(
            dIndex=D.index(x=1 - PCOSP_prob_good, surv.time=days_survived,
                surv.event=is_deceased, na.rm=TRUE, alpha=0.5,
                method.test='logrank'),
            cIndex=concordance.index(x=1 - PCOSP_prob_good, surv.time=days_survived,
                surv.event=is_deceased, method='noether', na.rm=TRUE)
        )
    )

    # assemble into a data.frame
    valStatsDF <- data.frame(
        statistic=c('D_index', 'concordance_index'),
        estimate=c(validationStats$dIndex$d.index, validationStats$cIndex$c.index),
        se=vapply(validationStats, `[[`, i='se', FUN.VALUE=numeric(1)),
        lower=vapply(validationStats, `[[`, i='lower', FUN.VALUE=numeric(1)),
        upper=vapply(validationStats, `[[`, i='upper', FUN.VALUE=numeric(1)),
        p_value=vapply(validationStats, `[[`, i='p.value', FUN.VALUE=numeric(1)),
        n=vapply(validationStats, `[[`, i='n', FUN.VALUE=numeric(1))
    )

    validationStats(PCOSPmodel) <- valStatsDF
    validationData(PCOSPmodel) <- CohortList(list(validationData),
        mDataTypes=metadata(validationData)$mDataType)
    return(PCOSPmodel)
})