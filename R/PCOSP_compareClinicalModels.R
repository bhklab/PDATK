#' Compare clinical models between patient cohorts
#'
## TODO:: HEEWON - description
#'
#' @param clinicalFeatures A
#' @param cohortClasses  A
#' @param cohorts A
#' @param models A \code{character} vector names for cohorts in `clinicalFeatures`
#'     to fit a linear model for. All models will be compared against each
#'     cohort in `cohorts`, or all cohorts if `cohorts` is not specified.
#' @param formula A
#'
#' @return A \code{list} with the first level representing the model the
#'   data was fit to, the second level representing the cohort the model was
#'   compared to, and the third level the associated statistics for that
#'   comparison.
#'
#' @importFrom reportROC reportROC
#' @export
compareClinicalModels <- function(clinicalFeatures, cohortClasses, cohorts,
                                  models, formula="binary_grp ~ Age + Sex +
                                  T_status + N + M + Grade")
{
    fitModels <- summarizeClinicalModels(clinicalFeatures, cohorts=models)

    cFeatures <- clinicalFeatures[cohorts]

    cClasses <- cohortClasses[cohorts]

    cohortProbs <- calculateCohortProbabilties(fitModels, cFeatures)

    modelComparisonStats <- .estimateModelAUCs(cFeatures, cohortProbs)

    return(modelComparisonStats)
}

#' Fit a generalized linear model for each clinical cohort in a list
#'
#' @param clinicalFeatures A list of \code{data.frame}s containing the clinical
#'     features for each cohort
#' @param cohorts A \code{character} vector of cohorts to subset from
#'      the clinicalFeatures argument. If excluded all cohorts are used.
#' @param formula A \code{character} vector containing the model formula
#'      expressed as a string. If excluded defaults to
#'      "binary_grp ~ Age + Sex + T_status +N + M + Grade".
#' @param ... Fallthrough parameters to `glm` function. If used, each cohort
#'      is passed as the data argument and formula as the formula. All other
#'      parameters must be specified in `...`.
#'
#' @return A named \code{list} containing the fitted model, model summary and anova
#'     for each cohort (named by `cohorts` argument or as
#'     `names(clincalFeatures)` if cohorts is absent).
#'
#' @export
summarizeClinicalModels <- function(clinicalFeatures, cohorts,
                                    formula="binary_grp ~ Age + Sex +
                                    T_status + N + M + Grade",
                                    ...)
{
    if (!missing(cohorts)) {
        cFeatures <- clinicalFeatures[cohorts]
        names(cFeatures) <- cohorts
    } else {
        cFeatures <- clinicalFeatures
    }

    ##TODO:: Add error handling check for formula with data columns
    if (!missing(...)) {
        models <- lapply(cFeatures,
                          function(data, formula)
                              glm(as.formula(formula), data=data, ...),
                          formula=formula)
    } else {
        models <- lapply(cFeatures,
                          function(data, formula)
                              glm(as.formula(formula), data=data,
                                  na.action=na.exclude,
                                  family=binomial(link="logit")),
                          formula=formula)
    }

    return(models)
}

#' Calculate the clinical and PCOSP probabiltie for a clinical cohort and
#'     a given clinical model
#'
#' @param clinicalCohorts A \code{list} of clinical cohorts to evaluate
#'     a predictive model against
#' @param model A \code{glm} model trained on a training cohort
#'
#' @keywords internal
.predictClinicalCohortProb <- function(clinicalCohort, model) {

    rownames(clinicalCohort) <- clinicalCohort$ID

    clinical <- predict(model, clinicalCohort, na.action=na.exclude,
                           type="response")

    pcosp <- 1 - clinicalCohort$pred_prob[which(clinicalCohort$ID %in%
                                                          names(clinical))]
    return(list(
        "clinical"=clinical,
        "PCOSP"=pcosp
    ))
}

#' Calculate the AUROC for a list of expression data matrices per clinical cohort
#'
#' @param cFeatures A \code{list} of per cohort clinical data.
#' @param cohortProbs A \code{list} of survival probability predictions
#'    for each cohort in `cFeatures`.
#'
#' @keywords internal
.estimateModelAUCs <- function (cFeatures, cohortProbs) {

    res <- list()
    for (i in seq_along(cohortProbs)) {
        res[[i]] <- structure(lapply(seq_along(cFeatures),
               function(j, cFeatures, cohortProbs)
                   .calcModelAUC(cFeatures[[j]], cohortProbs[[i]][[j]]),
               cFeatures=cFeatures,
               cohortProbs=cohortProbs),
               .Names=names(cFeatures))
    }
    return(structure(res, .Names=names(cohortProbs)))
}

#' Calculates the AUROCs for a clinical model
#'
#' @param coh A \code{matrix} or \code{data.frame} of
#'     expression data for a clinical cohort
#' @param preds A \code{vector} of predicted survival
#'     probabilities for the given cohort.
#'
#' @importFrom reportROC reportROC
#' @importFrom verification roc.area
.calcModelAUC <- function(coh, preds) {

    rownames(coh) <- coh$ID

    list(
        "roc1"=list(
            "roc"=reportROC(coh[names(preds$clinical),]$binary_grp, preds$clinical, plot=FALSE),
            "pval"=roc.area(coh[names(preds$clinical),]$binary_grp, preds$clinical)$p.value
        ),
        "roc2"=list(
            "roc"=reportROC(coh[names(preds$clinical),]$binary_grp, preds$PCOSP, plot=FALSE),
            "pval"=roc.area(coh[names(preds$clinical),]$binary_grp, preds$PCOSP)$p.value
        )
    )
}

#' Calcualte two different prediction probabilties under a given model for each
#'   cohort
#'
#' @param fitModels A \code{list} of clinical models, as returned by
#'     `summarizeClinicalModels`
#' @param clinicalCohorts A \code{list} of clinical cohorts to calculate
#'     probabilties under the given models
#'
#' @return A \code{list} of
#'
#' @import survival
#' @export
calculateCohortProbabilties <- function(fitModels, clinicalCohorts) {
    lapply(fitModels,
           function(model, clinicalCohorts)
               lapply(clinicalCohorts,
                      .predictClinicalCohortProb,
                      model=model),
           clinicalCohorts=clinicalCohorts)
}

#' Plot a bar plot comparing models for sequencing, microarray and overall cohorts
#'
##TODO:: HEEWON description
#'
#' @param modelComparisonStats A \code{list} of model comparison statistics as
#'     returned by `compareClinicalModels`.
#' @param model A \code{character} vector the name of the cohort the model
#'     was trained on, or a \code{numeric} vector with the integer index of
#'     the model. Defaults to 1, only need this if a model was trained on
#'     more than one cohort.
#' @param names A \code{character} vector of names for each paired barplot
#' @param colours A \code{character} vector of colours for the paired boxplots
#' @param filePath A
#' @param fileName A
#' @param ... Fallthrough arguments to `barplot`. Default values are used if
#'     this is exlcuded.
#'
#' @return Nothing, draws a plot
#'
#' @export
barplotModelComparison <- function(modelComparisonStats, model=1, names, colours, filePath, fileName, ...) {
    data <-data.frame(
        lapply(modelComparisonStats[[model]],
               function(cohort)
                   as.numeric(c(cohort$roc1$roc$AUC, cohort$roc2$roc$AUC))
        ))
    colnames(data) <- names
    rownames(data) <- c("Clinicopathological mode", "PCOSP")

    if (!missing(filePath) && !missing(fileName)) {
        pdf(file=file.path(filePath, paste0(fileName, ".pdf")))
        if (!missing(...)) {
            barplot(as.matrix(data), ...)
        } else {
            barplot(as.matrix(data), main="", ylim=c(0, 0.8), ylab="AUCs",
                    beside=TRUE, col=colours, cex.main=1.4, border="NA",
                    space=rep(c(0.6, 0.08), length(modelComparisonStats[[model]])))
        }
        dev.off()
    }

    if (!missing(...)) {
       barplot(as.matrix(data), ...)
    } else {
        barplot(as.matrix(data), main="", ylim=c(0, 0.8), ylab="AUCs",
                beside=TRUE, col=colours, cex.main=1.4, border="NA",
                space=rep(c(0.6, 0.08), length(modelComparisonStats[[model]])))
    }

}

#'
#'
#'
#' @importFrom survcomp combine.est
#' @export
metaEstimateComparisonAUCs <- function(modelComparisonStats, model=1, seqCohorts) {

    comparisonStats <- modelComparisonStats[[model]]

    isSeq <- grepl(paste(seqCohorts, collapse="|"), names(comparisonStats))

    pred1Estimates <- .metaEstimateModelAUCs(comparisonStats,
                                             subsets=list(isSeq, !isSeq,
                                                          c(which(isSeq), which(!isSeq))),
                                             subsetNames=c("Sequencing",
                                                           "Microarray",
                                                           "Overall"),
                                             na.rm=c(TRUE, FALSE, TRUE),
                                             hetero=c(TRUE, TRUE, TRUE),
                                             prediction=1)

    pred2Estimates <- .metaEstimateModelAUCs(comparisonStats,
                                             subsets=list(isSeq, !isSeq, c(which(isSeq), which(!isSeq))),
                                             subsetNames=c("Sequencing",
                                                           "Microarray",
                                                           "Overall"),
                                             na.rm=c(TRUE, FALSE, TRUE),
                                             hetero=c(TRUE, TRUE, TRUE),
                                             prediction=2)

    pred1Pvals <- lapply(pred1Estimates, function(subset)
        2*pnorm((subset$prediction$estimate - 0.5)/subset$prediction$se,
              lower.tail = subset$prediction$estimate < 0.5))

    pred2Pvals <- lapply(pred2Estimates, function(subset)
        2*pnorm((subset$prediction$estimate - 0.5)/subset$prediction$se,
                lower.tail = subset$prediction$estimate < 0.5))

    return(list(
        "clinical"=pred1Pvals ,
        "PCOSP"=pred2Pvals
    ))
}

#'
#'
#'
#'
#'
.metaEstimateModelAUCs <- function(comparisonStats, subsets, subsetNames,
                                   prediction, na.rm, hetero) {
    structure(lapply(seq_along(subsets),
           function(i, subsets, comparisonStats, prediction, na.rm, hetero)
               .calculateModelAUCs(comparisonStats[subsets[[i]]],
                                   modelIdx=prediction, na.rm=na.rm[i],
                                   hetero=hetero[i]),
           comparisonStats=comparisonStats,
           na.rm=na.rm,
           hetero=hetero,
           subsets=subsets,
           prediction=prediction),
           .Names=subsetNames)
}


#'
#'
#'
#'
#'
#'
.calculateModelAUCs <- function(comparisonStats, modelIdx, na.rm, hetero) {
        list(
            "prediction"=combine.est(
            vapply(comparisonStats,
                   function(cohort) as.numeric(cohort[[modelIdx]]$roc$AUC),
                   FUN.VALUE=numeric(1)),
            vapply(comparisonStats,
                   function(cohort) as.numeric(cohort[[modelIdx]]$roc$AUC.SE),
                   FUN.VALUE=numeric(1)),
            na.rm=na.rm,
            hetero=hetero),
            "cohorts"=names(comparisonStats)
    )
}
