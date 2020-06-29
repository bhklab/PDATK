#' Return
#'
#' @param validationCohorts A
#' @param validationStats A
#'
#' @return A named \code{list} of \code{list}s containing AUC, standard error of AUC,
#'     and the assocaited p-value for each cohort as well as summaries for
#'     sequencing cohorts, microarray cohorts as well as an overall summary.
#'     The list is named the same as validationCohorts.
#'
#' @export
aucMetaEstimates <- function(validationCohorts, validationStats,
                                 seqCohorts) {

    formattedValCohorts <- formatValidationCohorts(validationCohorts)

    PCOSPscores <- structure(lapply(validationStats$probabilities,
                                    function(cohort)
                                        cohort),
                             .Names=names(validationStats$probabilities))

    ## FIXME:: Can we ignore the wilcox.test ties warning?
    aucStats <- suppressWarnings(AUCstats(formattedValCohorts, PCOSPscores))

    isSeq <- names(formattedValCohorts) %in% seqCohorts

    aucStats$Overall <- .metaEstimateAUCstats(aucStats, hetero=TRUE)

    aucStats$Sequencing <- .metaEstimateAUCstats(aucStats[isSeq], hetero=FALSE)

    aucStats$Microarray <- .metaEstimateAUCstats(aucStats[!isSeq], hetero=FALSE)

    return(aucStats)
}

#' Calculate AUC and its standard error for a list of validation cohorts
#'
#' @param formattedValCohorts A \code{list} of validation cohorts.
#'    formatted using `formatValidationCohorts`.
#' @param probList A \code{list} of per cohort survival probability prediction.
#'
#' @importFrom pROC roc
#' @importFrom verification roc.area
#' @importFrom reportROC reportROC
#' @export
AUCstats <- function(formattedValCohorts, probList) {

    suppressMessages(
    structure(lapply(seq_along(formattedValCohorts),
                     function(i, cohorts, scores)
                         list("AUC"=roc(cohorts[[i]]$grp,
                                         scores[[i]][cohorts[[i]]$grpIndex])$auc[1],
                              "aucSE"=as.numeric(reportROC(cohorts[[i]]$grp,
                                                 scores[[i]][cohorts[[i]]$grpIndex],
                                                 plot=FALSE)$AUC.SE),
                              "pValue"=roc.area(cohorts[[i]]$grp,
                                                 1 - scores[[i]][cohorts[[i]]$grpIndex])$p.value
                              ),
                   cohorts=formattedValCohorts,
                   scores=probList),
                   .Names=names(formattedValCohorts)))
}

#' Generate a meta estimate of AUC for from a list of cohort statistics
#'
#' @param aucStats A \code{list} of AUC statistics for each cohort.
#' @param hetero A \code{boolean} indicating whether the cohorts are contain
#'     both sequencing and array cohorts.
#'
#' @return A \code{list} containing the meta-estimate for AUC, AUC standard
#'     error and the associated p-value.
#'
#' @importFrom survcomp combine.est
#' @keywords internal
.metaEstimateAUCstats <- function(aucStats, hetero) {
    metaEstimate <- combine.est(
        vapply(aucStats, function(stat) stat$AUC, FUN.VALUE=numeric(1)),
        vapply(aucStats, function(stat) stat$aucSE, FUN.VALUE=numeric(1)),
        hetero=hetero)
    list(
        "AUC"=metaEstimate$estimate,
        "aucSE"=metaEstimate$se,
        "pValue"=2 * pnorm((metaEstimate$estimate - 0.5)/metaEstimate$se,
                         lower.tail = metaEstimate$estimate < 0.5)
    )
}

