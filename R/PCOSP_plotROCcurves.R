#' Draw an ROC curve for each cohort in validationCohorts
#'
#' @param validationCohorts A
#' @param validationStats A
#' @param colours A
#' @param legendLoc A
#' @param filePath A
#' @param fileName A
#'
#' @return None, draws a plot
#'
#' @export
plotROCcurves <- function(formattedValidationCohorts, PCOSPscores,
                          colours, legendLoc="bottomright", filePath, fileName) {

    if (!missing(filePath) && !missing(fileName)) {
        pdf(file.path(filePath, paste0(fileName, ".pdf")))
        # Set default colours
        if (missing(colours)) {
            if (length(PCOSPscores < 13)) {
                colours <- brewer.pal(n=length(PCOSPscores), name="Set3")
            } else {
                colours <- c(brewer.pal(n=12, name="Set3"),
                             brewer.pal(n=length(PCOSPscores) - 12), name="Set2")
            }
        }

        ROCs <- .calculateROCs(formattedValCohorts, PCOSPscores)

        isSeq <- validationStats$isSequencing

        lineTypes <- ifelse(isSeq, 1, 3)

        for (i in seq_along(ROCs)) {
            if (i == 1) plot(ROCs[[i]], lwd=4, col=colours[i], lty=lineTypes[i])
            else plot(ROCs[[i]], lwd=4, col=colours[i], lty=lineTypes[i], add=TRUE)
        }

        cohortAUCstats <- AUCstats(formattedValCohorts, PCOSPscores)
        AUCs <- round(vapply(cohortAUCstats, function(cohortStats) cohortStats$AUC,
                             FUN.VALUE=numeric(1)), 2)
        pValues <- scales::scientific(vapply(cohortAUCstats, function(cohortStats) cohortStats$pValue,
                                             FUN.VALUE=numeric(1)), 2)
        cohortNames <- names(formattedValidationCohorts)

        legend(
            legendLoc,
            legend=paste0(cohortNames, ": ", AUCs, " (P = ", pValues, ")"),
            fill=colours,
            y.intersp=1,
            cex=0.9,
            bty="n"
        )
        dev.off()
    }
    if (interactive()) {

    }
    # Set default colours
    if (missing(colours)) {
        if (length(PCOSPscores < 13)) {
            colours <- brewer.pal(n=length(PCOSPscores), name="Set3")
        } else {
            colours <- c(brewer.pal(n=12, name="Set3"),
                         brewer.pal(n=length(PCOSPscores) - 12), name="Set2")
        }
    }

    ROCs <- .calculateROCs(formattedValCohorts, PCOSPscores)

    isSeq <- validationStats$isSequencing

    lineTypes <- ifelse(isSeq, 1, 3)

    for (i in seq_along(ROCs)) {
        if (i == 1) plot(ROCs[[i]], lwd=4, col=colours[i], lty=lineTypes[i])
        else plot(ROCs[[i]], lwd=4, col=colours[i], lty=lineTypes[i], add=TRUE)
    }

    cohortAUCstats <- AUCstats(formattedValCohorts, PCOSPscores)
    AUCs <- round(vapply(cohortAUCstats, function(cohortStats) cohortStats$AUC,
                   FUN.VALUE=numeric(1)), 2)
    pValues <- scales::scientific(vapply(cohortAUCstats, function(cohortStats) cohortStats$pValue,
                      FUN.VALUE=numeric(1)), 2)
    cohortNames <- names(formattedValidationCohorts)

    legend(
        legendLoc,
        legend=paste0(cohortNames, ": ", AUCs, " (P = ", pValues, ")"),
        fill=colours,
        y.intersp=1,
        cex=0.9,
        bty="n"
    )
}

#' @importFrom pROC roc
.calculateROCs <- function(formattedValCohorts, PCOSPscores) {
    structure(lapply(seq_along(formattedValCohorts),
                     function(i, cohorts, scores)
                         roc(cohorts[[i]]$grp,
                             scores[[i]][cohorts[[i]]$grpIndex]),
                     cohorts=formattedValCohorts,
                     scores=PCOSPscores),
              .Names=names(formattedValCohorts))
}
