#' Validate the Meta Estimate Calculation
#'
##TODO:: HEEWON - what does this function do?
#'
#' The function calculates D-index and Concordance index for independent cohorts
#'   and calculates the meta-estimate of C-index and D-index. The forestplot of
#'   C-index and D-index was plotted for all the cohorts.
#'
#' @example
#' # Load validation data
#' data(validationCohorts)
#' data(selectedModels)
#'
#TODO:: Set nthread to 1
#' # Validate the meta-estimate calculation
#' validationStats <- validateMetaEstimateCalculation(validationCohorts[1], selectedModels,
#'     seqCohorts=c("PCSI", "TCGA", "Kirby"), nthread=2)
#'
#' @param validationCohorts A named \code{list} of validation cohorts.
#' @param selectedModels A \code{list} of selected models from the
#'     `buildPCOSPmodels` function.
#' @param seqCohorts A \code{character} vector specifying which cohorts
#'     are from sequencing data. It is assumed that non-sequencing cohorts
#'     are from micro-array cohorts.
#' @param nthread A \code{numeric} vector with an integer number of thread to
#'     parallelize across.
#' @param saveDir A \code{character} vector specifying the directory in which
#'     to save results. If missing the function returns the data.
#'
#' @return A named \code{list} of statistical meta-estimates for the combined,
#'   sequencing and microarray cohorts supplied in validationCohorts. Top level
#'   list labels indicate which cohorts, with the subsequent list level
#'   indicating the type of statistic (i.e., dindex vs concordance index)
#'
#' @importFrom reportROC reportROC
#' @importFrom pROC roc
#' @importFrom verification roc.area
#' @importFrom survcomp D.index concordance.index combine.est
#' @import survival
#' @export
calculateValidationStats <- function(validationCohorts, selectedModels, seqCohorts,
                                nthread, saveDir) {

  probList <- calculatePCOSPscores(validationCohorts, selectedModels, nthread)

  constructMetaEstimatesDF(probList, validationCohorts, seqCohorts)
}

#' Calculate Dindex, Cindex and cohort meta-estimates from a list of cohorts
#'     and associated probabilities.
#'
#' @param probList A \code{list} of probability vectors, one for each cohort
#' @param cohortData A \code{list} of matrices or data.frames containg per cohort
#'   survival data. Columns OS_Stats and OS are assumed to be present!
#' @param seqCohorts A \code{character} vector with the names of the cohorts
#'   which contain sequencing data. Must matcht he names attribute of  probList
#'   and cohortData.
#' @param hetero A \code{logocal} vector indicating if the `hetero` parameter
#'   should be TRUE or FALSE in the combined summary, sequencing summary and
#'   microarray summary, respectively. Defaults to TRUE, FALSE, FALSE.
#'
#' @return A \code{list} with items dIndex, cIndex, probabilities and isSequencing where
#'   `dIndex`/`cIndex`` contains a data.frame with statistics for each cohort and summary
#'   and `probabilites`/`isSequencing` contains the probabilities used in the statistical
#'   calcualte and wheter or not a respective cohort contains sequening data.
#'
#' @importFrom survcomp D.index concordance.index combine.est
#' @importFrom reportROC reportROC
#' @importFrom pROC roc
#' @importFrom verification roc.area
#' @import survival
#' @export
constructMetaEstimatesDF <- function(probList, cohortData, seqCohorts, hetero=c(TRUE, FALSE, FALSE)) {

  ## Dindex estimate calculation
  DindexList <- .estimateDindex(probList, cohortData)

  ## Concordance index calculation
  concordanceIndexList <- .estimateConcordanceIndex(probList,
                                                    cohortData)

  ## Calculating meta-estimates of D-index and Concordance index

  ##  Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR OVERALL DATA
  combinedStats <- metaEstimateStats(DindexList, concordanceIndexList,
                                     hetero=hetero[1])

  ## Determine which cohorts are from sequencing data
  isSeq = grepl(paste(seqCohorts, collapse="|"), names(cohortData))

  ## Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR sequencing cohort
  sequencingStats <- metaEstimateStats(DindexList[isSeq],
                                       concordanceIndexList[isSeq],
                                       hetero=hetero[2])

  ## Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR microarray cohort
  arrayStats <- metaEstimateStats(DindexList[!isSeq],
                                  concordanceIndexList[!isSeq],
                                  hetero=hetero[3])

  # Extract statistics for Dindex and concordanceIndex into a list of data.frames
  list(
    "dIndex"=data.frame(
      "mean"=c(vapply(DindexList, function(x) log2(x$d.index), FUN.VALUE=numeric(1)),
               vapply(list(sequencingStats, arrayStats, combinedStats),
                      function(x) log2(x$metaEstimate$dIndex$estimate), FUN.VALUE=numeric(1))),
      "lower"=c(vapply(DindexList, function(x) log2(x$lower), FUN.VALUE=numeric(1)),
                vapply(list(sequencingStats, arrayStats, combinedStats),
                       function(x) log2(x$lowerTail$dIndex), FUN.VALUE=numeric(1))),
      "upper"=c(vapply(DindexList, function(x) log2(x$upper), FUN.VALUE=numeric(1)),
                vapply(list(sequencingStats, arrayStats, combinedStats),
                       function(x) log2(x$upperTail$dIndex), FUN.VALUE=numeric(1))),
      "pval"=c(vapply(DindexList, function(x) x$p.value, FUN.VALUE=numeric(1)),
               vapply(list(sequencingStats, arrayStats, combinedStats),
                      function(x) x$pValue$dIndex, FUN.VALUE=numeric(1))),
      row.names=c(names(DindexList), "Sequencing", "Microarray", "Overall")
    ),
    "cIndex"=data.frame(
      "mean"=c(vapply(concordanceIndexList, function(x) x$c.index, FUN.VALUE=numeric(1)),
               vapply(list(sequencingStats, arrayStats, combinedStats),
                      function(x) x$metaEstimate$cIndex$estimate, FUN.VALUE=numeric(1))),
      "lower"=c(vapply(concordanceIndexList, function(x) x$lower, FUN.VALUE=numeric(1)),
                vapply(list(sequencingStats, arrayStats, combinedStats),
                       function(x) x$lowerTail$cIndex, FUN.VALUE=numeric(1))),
      "upper"=c(vapply(concordanceIndexList, function(x) x$upper, FUN.VALUE=numeric(1)),
                vapply(list(sequencingStats, arrayStats, combinedStats),
                       function(x) x$upperTail$cIndex, FUN.VALUE=numeric(1))),
      "pval"=c(vapply(concordanceIndexList, function(x) x$p.value, FUN.VALUE=numeric(1)),
               vapply(list(sequencingStats, arrayStats, combinedStats),
                      function(x) x$pValue$cIndex, FUN.VALUE=numeric(1))),
      row.names=c(names(concordanceIndexList), "Sequencing", "Microarray", "Overall")
    ),
    "probabilities"=probList,
    "isSequencing"=isSeq
  )
}

#' Calculate the cohort-wise PCOSP score from a list of validation cohorts
#'
#' @param validationCohorts A \code{list} of validation cohorts
#' @param selectedModels A \code{list} of selected models, as returned by the
#'   `buildPCOSPmodels` function.
#' @param nthread A \code{numeric} vector containing the integer number of threads
#'     to a parallelize across.
#'
#' @return A \code{list} of PCOSP scores for each validation cohorot
#'
#' @importFrom survcomp D.index concordance.index combine.est
#' @importFrom reportROC reportROC
#' @importFrom pROC roc
#' @importFrom verification roc.area
#' @import survival
#' @export
calculatePCOSPscores <- function(validationCohorts, selectedModels, nthread) {

  formattedValCohorts <- formatValidationCohorts(validationCohorts)

  ## Extract matrices from the cohort list
  cohortMatrixList <- lapply(formattedValCohorts, function(cohort) cohort$mat)
  ##TODO:: This is not used here, do we need it?
  cohortGroupList <- lapply(formattedValCohorts, function(cohort) cohort$grp)

  # Estimate PCOSP scores from cohort matrixes
  PCOSPscoreList <- lapply(cohortMatrixList,
                           function(cohortMat, selectedModels, nthread)
                             structure(estimatePCOSPprob(cohortMat,
                                                         selectedModels,
                                                         nthread),
                                       .Names=rownames(cohortMat)),
                           selectedModels=selectedModels, nthread=nthread)
  return(PCOSPscoreList)
}

#' Predict meta-estimates of Dindex and Cindex from an equal length list of
#'     each.
#'
#' @param DindexList A
#' @param concordanceIndexList A
#' @param hetero A
#'
#' @importFrom survcomp D.index concordance.index combine.est
#' @import survival
#' @export
metaEstimateStats <- function(DindexList, concordanceIndexList, hetero) {
  DindexMetaEstimate <- .metaEstimateDindex(DindexList, hetero)
  concordanceIndexMetaEstimate <-
    .metaEstimateConcordanceIndex(concordanceIndexList, hetero)
  stats <- .zipLists(DindexMetaEstimate,
                     concordanceIndexMetaEstimate,
                     sublistNames=c("dIndex", "cIndex"))
  return(stats)
}

#' Estimate the D-index for a list of validation cohorts
#'
#' @param probList A named \code{list} of survival probabilities as calculated
#'   with `estimatePCOSPscore`
#' @param validationCohorts A named \code{list} of validation cohorts
#'
#' @return A \code{list} of per cohort D index statistics
#'    calculated with the `D.index` function from the
#'    `survcomp` package.
#'
#' @importFrom survcomp D.index
#' @import survival
#' @importFrom survival strata
#'
.estimateDindex <- function(probList, validationCohorts) {
  require(survival)
  structure(lapply(seq_along(probList), function(i, probList, validationCohorts) {
    cohort <- validationCohorts[[i]]
    probScore <- probList[[i]]
    D.index(x=probScore,
            surv.time=as.numeric.factor(cohort[, 'OS']),
            surv.event=as.numeric.factor(cohort[, 'OS_Status']),
            na.rm=TRUE,
            alpha=0.05,
            method.test="logrank")
  },
  validationCohorts=validationCohorts,
  probList=probList),
  .Names=names(probList))
}

#' Estimate the concordance index for a list of validation cohorts
#'
#' @param probList A named \code{list} of survival probabilities as calculated
#'   with `estimatePCOSPscore`
#' @param validationCohorts A named \code{list} of validation cohorts
#'
#' @return A \code{list} of per cohort concordance index statistics
#'    calculated with the `concordance.index` function from the
#'    `survcomp` package.
#'
#' @importFrom survcomp concordance.index
#' @import survival
#' @importFrom survival strata
.estimateConcordanceIndex <- function(probList, validationCohorts) {
  require(survival)
  structure(lapply(seq_along(validationCohorts), function(i) {
    cohort <- validationCohorts[[i]]
    PCOSPscore <- probList[[i]]
    concordance.index(x=PCOSPscore,
            surv.time=as.numeric.factor(cohort[, 'OS']),
            surv.event=as.numeric.factor(cohort[, 'OS_Status']),
            method="noether")
  }), .Names=names(probList))
}

#' Meta-estimate the Dindex for a list of validation cohorts
#'
#' @param DindexList A \code{list} of per cohort D indexes
#' @param hetero A \code{logical} vector passed to the
#'    hetero argument of `combine.est` from the `survcomp`
#'    package. See `?combine.est` for more details.
#'
#' @importFrom survcomp combine.est
#' @import survival
#'
.metaEstimateDindex <- function(DindexList, hetero) {

  Dindexes <- vapply(DindexList, function(cohort) cohort$d.index,
                     FUN.VALUE=numeric(1))
  DindexSEs <- vapply(DindexList, function(cohort) cohort$se,
                      FUN.VALUE=numeric(1))

  DindexMetaEstimate <- combine.est(Dindexes, DindexSEs, na.rm=TRUE, hetero=hetero)

  ##TODO:: Define a metaestimate S4 object? Can then dispatch plots on it
  return(list(
    "metaEstimate"= DindexMetaEstimate,
    "lowerTail"=DindexMetaEstimate$estimate + qnorm(0.025, lower.tail=TRUE) *
      DindexMetaEstimate$se,
    "upperTail"=DindexMetaEstimate$estimate + qnorm(0.025, lower.tail=FALSE) *
      DindexMetaEstimate$se,
    "se"=DindexMetaEstimate$se,
    "pValue"=2*pnorm(-abs(log(DindexMetaEstimate$estimate)/DindexMetaEstimate$se)),
    "cohortNames"=names(DindexList)
  ))
}

#' Meta-estimate the concordance index for a list of validation cohorts
#'
#' @param concordanceIndexList A \code{list} of per cohort concordance indexes
#' @param hetero A \code{logical} vector passed to the
#'    hetero argument of `combine.est` from the `survcomp`
#'    package. See `?combine.est` for more details.
#'
#' @importFrom survcomp concordance.index combine.est
#' @import survival
#'
.metaEstimateConcordanceIndex <- function(concordanceIndexList, hetero) {

  conIndexes <- vapply(concordanceIndexList, function(cohort) cohort$c.index,
                     FUN.VALUE=numeric(1))
  conIndexSEs <- vapply(concordanceIndexList, function(cohort) cohort$se,
                      FUN.VALUE=numeric(1))

  conIndexMetaEstimate <- combine.est(conIndexes, conIndexSEs, na.rm=TRUE, hetero=hetero)

  ##TODO:: Define a metaestimate S4 object? Can then dispatch plots on it
  return(list(
    "metaEstimate"= conIndexMetaEstimate,
    "lowerTail"=conIndexMetaEstimate$estimate + qnorm(0.025, lower.tail=TRUE) *
      conIndexMetaEstimate$se,
    "upperTail"=conIndexMetaEstimate$estimate + qnorm(0.025, lower.tail=FALSE) *
      conIndexMetaEstimate$se,
    "se"=conIndexMetaEstimate$se,
    "pValue"=2*pnorm((conIndexMetaEstimate$estimate - 0.5)/conIndexMetaEstimate$se,
                     lower.tail= conIndexMetaEstimate$estimate < 0.5),
    "cohortNames"=names(concordanceIndexList)
  ))
}