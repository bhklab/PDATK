#' Return a list of lists of classification data.frames, top level is sequencing, microarray, overall;
#'     subsequent level is both, basal, classical
#'
#' @param probList A \code{list} of per cohort survival probability predictions
#' @param cohortData A \code{list} of per cohort survival data
#' @param cohortClasses A \code{list} of per cohort subtype classifications.
#' @param seqCohorts A \code{character} vector of names for the cohorts
#'    containing sequencing data. Other cohorts are assumed to contain
#'    microarray data.
#' @param hetero A \code{logical} vector passed to `combine.est` for the
#'    combined, sequencing and microarray cohorts, respectively.
#'
#' @import survival
#' @export
calculateCohortSubtypeStats <- function(probList, cohortData, cohortClasses, seqCohorts, hetero) {

    isSeq= names(probList) %in% seqCohorts

    if (missing(hetero)) {
        list(
            "Sequencing"=constructSubtypeMetaEstimatesDF(probList[isSeq], cohortData[isSeq], cohortClasses[isSeq]),
            "Microarray"=constructSubtypeMetaEstimatesDF(probList[!isSeq], cohortData[!isSeq], cohortClasses[!isSeq]),
            "Overall"=constructSubtypeMetaEstimatesDF(probList, cohortData, cohortClasses)
        )
    } else {
        list(
            "Sequencing"=constructSubtypeMetaEstimatesDF(probList[isSeq], cohortData[isSeq], cohortClasses[isSeq], hetero=hetero[1]),
            "Microarray"=constructSubtypeMetaEstimatesDF(probList[!isSeq], cohortData[!isSeq], cohortClasses[!isSeq], hetero=hetero[2]),
            "Overall"=constructSubtypeMetaEstimatesDF(probList, cohortData, cohortClasses, hetero=hetero[3])
        )
    }
}


##FIXME:: Refactor into validation stats calculations
##FIXME:: Refactor this monster to construct the DF's in a separate function call
#' Calculate stats for combined, basal and classical subtypes for each cohort and return a list of data.frames
#'
#' @param probList A \code{list} of per cohort survival probability predictions
#' @param cohortData A \code{list} of per cohort survival data
#' @param cohortClasses A \code{list} of per cohort subtype classifications.
#' @param hetero A \code{logical} vector passed to `combine.est` for the
#'    combined, basal and classical subtypes, respectively.
#'
#' @import survival
#' @export
constructSubtypeMetaEstimatesDF <- function(probList, cohortData, cohortClasses, hetero=c(TRUE, TRUE, TRUE)) {

    ### All

    ## Dindex estimate calculation
    DindexList <- .estimateDindex(probList, cohortData)

    ## Concordance index calculation
    concordanceIndexList <- .estimateConcordanceIndex(probList,
                                                      cohortData)
    ### Basal Subtype

    basalData <- .extractSubtypeData(probList, cohortData, cohortClasses, "basal")

    DindexListBasal <- .estimateDindex(basalData$probs, basalData$surv)

    concordanceIndexListBasal <- .estimateConcordanceIndex(basalData$probs, basalData$surv)

    ### Classical subtype

    classicalData <- .extractSubtypeData(probList, cohortData, cohortClasses, "classical")

    DindexListClassical <- .estimateDindex(classicalData$probs, classicalData$surv)

    concordanceIndexListClassical <- .estimateConcordanceIndex(classicalData$probs, classicalData$surv)

    ## Calculating meta-estimates of D-index and Concordance index

    ##  Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR OVERALL DATA
    combinedStats <- metaEstimateStats(DindexList, concordanceIndexList,
                                       hetero=hetero[1])

    ## Basal
    basalStats <- metaEstimateStats(DindexListBasal, concordanceIndexListBasal,
                                    hetero=hetero[2])

    ## Classical
    classicalStats <- metaEstimateStats(DindexListClassical, concordanceIndexListClassical,
                                        hetero=hetero[3])

    # Extract statistics for Dindex and concordanceIndex into a list of data.frames
    list(
        "combined"=list(
            "dIndex"=data.frame(
                "mean"=c(vapply(DindexList, function(x) log2(x$d.index), FUN.VALUE=numeric(1)),
                         vapply(list(combinedStats),
                                function(x) log2(x$metaEstimate$dIndex$estimate), FUN.VALUE=numeric(1))),
                "lower"=c(vapply(DindexList, function(x) log2(x$lower), FUN.VALUE=numeric(1)),
                          vapply(list(combinedStats),
                                 function(x) log2(x$lowerTail$dIndex), FUN.VALUE=numeric(1))),
                "upper"=c(vapply(DindexList, function(x) log2(x$upper), FUN.VALUE=numeric(1)),
                          vapply(list(combinedStats),
                                 function(x) log2(x$upperTail$dIndex), FUN.VALUE=numeric(1))),
                "pval"=c(vapply(DindexList, function(x) x$p.value, FUN.VALUE=numeric(1)),
                         vapply(list(combinedStats),
                                function(x) x$pValue$dIndex, FUN.VALUE=numeric(1))),
                row.names=c(names(DindexList), "combined")
            ),
            "cIndex"=data.frame(
                "mean"=c(vapply(concordanceIndexList, function(x) x$c.index, FUN.VALUE=numeric(1)),
                         vapply(list(combinedStats),
                                function(x) x$metaEstimate$cIndex$estimate, FUN.VALUE=numeric(1))),
                "lower"=c(vapply(concordanceIndexList, function(x) x$lower, FUN.VALUE=numeric(1)),
                          vapply(list(combinedStats),
                                 function(x) x$lowerTail$cIndex, FUN.VALUE=numeric(1))),
                "upper"=c(vapply(concordanceIndexList, function(x) x$upper, FUN.VALUE=numeric(1)),
                          vapply(list(combinedStats),
                                 function(x) x$upperTail$cIndex, FUN.VALUE=numeric(1))),
                "pval"=c(vapply(concordanceIndexList, function(x) x$p.value, FUN.VALUE=numeric(1)),
                         vapply(list(combinedStats),
                                function(x) x$pValue$cIndex, FUN.VALUE=numeric(1))),
                row.names=c(names(concordanceIndexList), "combined")
                ),
            "probabilities"=probList,
            "survival"=cohortData,
            "subtype"="combined"
        ),
        "basal"=list(
            "dIndex"=data.frame(
                "mean"=c(vapply(DindexListBasal, function(x) log2(x$d.index), FUN.VALUE=numeric(1)),
                         vapply(list(basalStats),
                                function(x) log2(x$metaEstimate$dIndex$estimate), FUN.VALUE=numeric(1))),
                "lower"=c(vapply(DindexListBasal, function(x) log2(x$lower), FUN.VALUE=numeric(1)),
                          vapply(list(basalStats),
                                 function(x) log2(x$lowerTail$dIndex), FUN.VALUE=numeric(1))),
                "upper"=c(vapply(DindexListBasal, function(x) log2(x$upper), FUN.VALUE=numeric(1)),
                          vapply(list(basalStats),
                                 function(x) log2(x$upperTail$dIndex), FUN.VALUE=numeric(1))),
                "pval"=c(vapply(DindexListBasal, function(x) x$p.value, FUN.VALUE=numeric(1)),
                         vapply(list(basalStats),
                                function(x) x$pValue$dIndex, FUN.VALUE=numeric(1))),
                row.names=c(names(DindexListBasal), "basal")
            ),
            "cIndex"=data.frame(
                "mean"=c(vapply(concordanceIndexListBasal, function(x) x$c.index, FUN.VALUE=numeric(1)),
                         vapply(list(basalStats),
                                function(x) x$metaEstimate$cIndex$estimate, FUN.VALUE=numeric(1))),
                "lower"=c(vapply(concordanceIndexListBasal, function(x) x$lower, FUN.VALUE=numeric(1)),
                          vapply(list(basalStats),
                                 function(x) x$lowerTail$cIndex, FUN.VALUE=numeric(1))),
                "upper"=c(vapply(concordanceIndexListBasal, function(x) x$upper, FUN.VALUE=numeric(1)),
                          vapply(list(basalStats),
                                 function(x) x$upperTail$cIndex, FUN.VALUE=numeric(1))),
                "pval"=c(vapply(concordanceIndexListBasal, function(x) x$p.value, FUN.VALUE=numeric(1)),
                         vapply(list(basalStats),
                                function(x) x$pValue$cIndex, FUN.VALUE=numeric(1))),
                row.names=c(names(concordanceIndexListBasal), "basal")
            ),
            "probabilities"=basalData$probs,
            "survival"=basalData$surv,
            "subtype"="basal"
        ),
        "classical"=list(
            "dIndex"=data.frame(
                "mean"=c(vapply(DindexListClassical, function(x) log2(x$d.index), FUN.VALUE=numeric(1)),
                         vapply(list(classicalStats),
                                function(x) log2(x$metaEstimate$dIndex$estimate), FUN.VALUE=numeric(1))),
                "lower"=c(vapply(DindexListClassical, function(x) log2(x$lower), FUN.VALUE=numeric(1)),
                          vapply(list(classicalStats),
                                 function(x) log2(x$lowerTail$dIndex), FUN.VALUE=numeric(1))),
                "upper"=c(vapply(DindexListClassical, function(x) log2(x$upper), FUN.VALUE=numeric(1)),
                          vapply(list(classicalStats),
                                 function(x) log2(x$upperTail$dIndex), FUN.VALUE=numeric(1))),
                "pval"=c(vapply(DindexListClassical, function(x) x$p.value, FUN.VALUE=numeric(1)),
                         vapply(list(classicalStats),
                                function(x) x$pValue$dIndex, FUN.VALUE=numeric(1))),
                row.names=c(names(DindexListClassical), "classical")
            ),
            "cIndex"=data.frame(
                "mean"=c(vapply(concordanceIndexListClassical, function(x) x$c.index, FUN.VALUE=numeric(1)),
                         vapply(list(classicalStats),
                                function(x) x$metaEstimate$cIndex$estimate, FUN.VALUE=numeric(1))),
                "lower"=c(vapply(concordanceIndexListClassical, function(x) x$lower, FUN.VALUE=numeric(1)),
                          vapply(list(classicalStats),
                                 function(x) x$lowerTail$cIndex, FUN.VALUE=numeric(1))),
                "upper"=c(vapply(concordanceIndexListClassical, function(x) x$upper, FUN.VALUE=numeric(1)),
                          vapply(list(classicalStats),
                                 function(x) x$upperTail$cIndex, FUN.VALUE=numeric(1))),
                "pval"=c(vapply(concordanceIndexListClassical, function(x) x$p.value, FUN.VALUE=numeric(1)),
                         vapply(list(classicalStats),
                                function(x) x$pValue$cIndex, FUN.VALUE=numeric(1))),
                row.names=c(names(concordanceIndexListClassical), "classical")
            ),
            "probabilities"=classicalData$prob,
            "surival"=classicalData$surv,
            "subtype"="classical"
        )
    )
}

#' Subset a list of probabilities and survival data to samples in a subtype
#'
#' @param probList A \code{list} of per cohort survival probability predictions
#' @param cohortData A \code{list} per cohort survival data
#' @param cohortSubtypes A \code{list} of subtypes per clinical cohort
#' @param subtype A \code{character} vector indicating which subtype
#'    to subset on.
#'
#' @keywords internal
.extractSubtypeData <- function(probList, cohortData, cohortSubtypes, subtype) {
    # Get sample names for the selected subtype
    sampleSubtypes <- lapply(cohortSubtypes,
                             function(cohort, subtype) cohort[cohort[, 2] == subtype, 1],
                             subtype=subtype)

    # Subset cohort probalities to subtype samples
    subProbList <- structure(lapply(seq_along(probList),
                          function(i, cohorts, sampSubtypes)
                              cohorts[[i]][names(cohorts[[i]]) %in% sampSubtypes[[i]]],
                          sampSubtypes=sampleSubtypes,
                          cohorts=probList),
                          .Names=names(probList))

    # Subset survival data to subtype samples
    subSurvival <- structure(lapply(seq_along(cohortData),
                          function(i, cohorts, sampSubtypes)
                              cohorts[[i]][rownames(cohorts[[i]]) %in% sampSubtypes[[i]], ],
                          sampSubtypes=sampleSubtypes,
                          cohorts=cohortData),
                          .Names=names(cohortData))
    return(list(
        "probs"=subProbList,
        "surv"=subSurvival
    ))
}

#'
#'
#'
#'
#' @export
.formatCohortSubtypeStatsForPlot <- function(cohortSubtypeStats, stat) {
    # Extract data.frame for each platform
    seq <- cohortSubtypeStats$Sequencing
    array <- cohortSubtypeStats$Microarray
    overall <- cohortSubtypeStats$Overall

    # Merge statistics with rownames indicating subtype
    overallStat <- lapply(overall,
                          function(platform, stat) {
                              mat <- t(as.matrix(platform[[stat]], rownames.force=TRUE))
                              rownames(mat) <- paste0(rownames(mat), "_", platform$subtype)
                              mat
                          },
                          stat=stat)

    overallStatMat <- do.call(rbind, overallStat)

    seqStat <- lapply(seq,
                          function(platform, stat) {
                              mat <- t(as.matrix(platform[[stat]], rownames.force=TRUE))
                              rownames(mat) <- paste0(rownames(mat), "_", platform$subtype)
                              mat
                          },
                          stat=stat)

    seqStatMat <- do.call(rbind, seqStat)

    arrayStat <- lapply(array,
                          function(platform, stat) {
                              mat <- t(as.matrix(platform[[stat]], rownames.force=TRUE))
                              rownames(mat) <- paste0(rownames(mat), "_", platform$subtype)
                              mat
                          },
                          stat=stat)

    arrayStatMat <- do.call(rbind, arrayStat)

    return(list(
        "Overall"=overallStatMat,
        "Sequencing"=seqStatMat,
        "Mircoarray"=arrayStatMat
    ))
}






























