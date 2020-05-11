#' Draw a forest plot using the statistics calcualted in validation metaestimate
#'
#' @examples
#'
#' @param validationStatsDF A \code{data.frame} containing statistics for Dindex
#'     or concordance index, as available in the \code{list} returned by
#'     `validateMetaEstimateCalculation`.
#' @param stat A \code{character} vector containing the statistic to
#'     use in the forest plot. Options are "dIndex" or "cIndex".
#' @param isSummary A \code{logical} vector indicating which cohorts contain
#'     summary data. Must be one boolean value per row in the data.frame.
#' @param filePath A \code{character} vector containing the file path to write
#'     the plotting results to.
#' @param fileName A \code{character} vector containing the desired file name
#'     with the desired file extension. Supported extensions include .pdf, .png,
#'     .svg and .jpeg, for more information see documentation for the `ggsave()`
#'     function from `ggplot2`.
#' @param ... Additional arguments passed to forestplot.
#'
#' @import scales
#' @import forestplot
#' @importFrom grid unit grid.grabExpr grid.draw gpar
#' @importFrom ggplot2 ggsave
#' @export
forestPlotMetaEstimate <- function(validationStats, stat, isSummary, filePath,
                                   fileName, ...) {

    # Extract necessary statistics for plotting
    validationStatsDF <- validationStats[[stat]]
    validationStatsDF <- rbind(rep(NA, ncol(validationStatsDF)),
                               validationStatsDF)
    PCOSPscoreList <- validationStats$probabilities
    isSeq <- validationStats$isSequencing
    isSummary <- c(TRUE, isSummary)

    # Construct the forest plot table
    labelText <- data.frame(
        "cohort"=c("Cohorts", rownames(validationStatsDF)[-1]),
        "pvalue"=c("P value", scales::scientific(validationStatsDF$pval[-1], 2))
    )

    # Extract box sizes
    lengthSeq <- length(unlist(PCOSPscoreList[isSeq]))
    lengthArray <- length(unlist(PCOSPscoreList[!isSeq]))
    boxSizes <- c(NA, vapply(PCOSPscoreList, function(score) length(unlist(score)),
                         FUN.VALUE=numeric(1)),
                  lengthSeq, lengthArray,
                  lengthSeq + lengthArray) / 1000

    # Match correct plot function to call
    if(missing(...)) {
        if (stat == "dIndex") {
            plot <- .forestPlotDindex(labelText, validationStatsDF, isSeq,
                                      isSummary, boxSizes)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex(labelText, validationStatsDF, isSeq,
                                      isSummary, boxSizes)
        } else {
            stop(paste0("There is no statistic called: ", stat))
        }
    # Allow user to specify custom parameters
    } else {
        if (stat == "dIndex") {
            plot <- .forestPlotDindex(labelText, validationStatsDF, isSeq,
                                      isSummary, boxSizes, ...)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex(labelText, validationStatsDF, isSeq,
                                      isSummary, boxSizes, ...)
        } else {
            stop(paste0("There is no statistic called: ", stat))
        }
    }
    # Decide whether to plot to device or save to disk
    if (missing(filePath) || missing(fileName)) {
        grid.draw(plot)
    } else {
        grid.draw(plot)
        ggsave(filename=fileName, path=filePath, plot=plot)
    }
}


.forestPlotDindex <- function(labelText, validationStatsDF, isSeq, isSummary, boxSizes,
                              ...) {

    if (!missing(...)) {
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                               validationStatsDF[, c("mean", "lower", "upper")],
                               ...))
    } else {
        # Set plot colouring functions
        # ##TODO:: Determine if there is a more readable way to write this?
        normFun <- local({
            i = 0
            isSeq=isSeq
            b_clrs=ifelse(isSeq, "#FF7F00", "#1F78B4")
            l_clrs=ifelse(isSeq, "#FF7F00", "#1F78B4")
            function(..., clr.line, clr.marker){
                i <<- i + 1
                fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
            }
        })
        sumFun <- local({
            i = 0
            s_clrs =c("#FF7F00","#1F78B4","grey57")
            function(..., col){
                i <<- i + 1
                fpDrawSummaryCI(...,col=s_clrs[i])
            }
        })
        # Make the plot
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                               validationStatsDF[, c("mean", "lower", "upper")],
                               xlab="Log2 HR",
                               is.summary=isSummary,
                               clip=c(-1, 2.5),
                               txt_gp=fpTxtGp(label=gpar(fontfamily="Helvetica"),
                                              ticks=gpar(cex=0.8),
                                              xlab=gpar(fontfamily="Helvetica",
                                                        cex=1)),
                               col=fpColors(box="black"),
                               title=" ",
                               new_page=FALSE,
                               fn.ci_norm=normFun,
                               fn.ci_sum=sumFun,
                               zero=0,
                               graphwidth=unit(2, "inches"),
                               align=c("l"),
                               pch=16,
                               boxsize=boxSizes + 0.2))
    }
    return(plot)
}

.forestPlotCindex <- function(labelText, validationStatsDF, isSeq, isSummary, boxSizes,
                              ...) {
    if (!missing(...)) {
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     validationStatsDF[, c("mean", "lower", "upper")],
                                                     ...))
    } else {
        # Set plot colouring functions
        # ##TODO:: Determine if there is a more readable way to write this?
        normFun <- local({
            i = 0
            isSeq=isSeq
            b_clrs=ifelse(isSeq, "#FF7F00", "#1F78B4")
            l_clrs=ifelse(isSeq, "#FF7F00", "#1F78B4")
            function(..., clr.line, clr.marker){
                i <<- i + 1
                fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
            }
        })
        sumFun <- local({
            i = 0
            s_clrs =c("#FF7F00","#1F78B4","grey57")
            function(..., col){
                i <<- i + 1
                fpDrawSummaryCI(...,col=s_clrs[i])
            }
        })
        # Make the plot
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     validationStatsDF[, c("mean", "lower", "upper")],
                                                     xlab="C-index",
                                                     is.summary=isSummary,
                                                     clip=c(0.3, 0.8),
                                                     txt_gp=fpTxtGp(label=gpar(fontfamily="Helvetica"),
                                                                    ticks=gpar(cex=0.9),
                                                                    xlab=gpar(fontfamily="Helvetica",
                                                                              cex=1)),
                                                     col=fpColors(box="black"),
                                                     title=" ",
                                                     new_page=FALSE,
                                                     fn.ci_norm=normFun,
                                                     fn.ci_sum=sumFun,
                                                     zero=0.5,
                                                     graphwidth=unit(2, "inches"),
                                                     align=c("l"),
                                                     pch=16,
                                                     boxsize=boxSizes + 0.2))
    }
    return(plot)
}


#' Draw a forest plot using the statistics calcualted in validation metaestimate
#'
#' @examples
#'
#' @param validationStatsDF A \code{data.frame} containing statistics for Dindex
#'     or concordance index, as available in the \code{list} returned by
#'     `validateMetaEstimateCalculation`.
#' @param stat A \code{character} vector containing the statistic to
#'     use in the forest plot. Options are "dIndex" or "cIndex".
#' @param isSummary A \code{logical} vector indicating which cohorts contain
#'     summary data. Must be one boolean value per row in the data.frame.
#' @param filePath A \code{character} vector containing the file path to write
#'     the plotting results to.
#' @param fileName A \code{character} vector containing the desired file name
#'     with the desired file extension. Supported extensions include .pdf, .png,
#'     .svg and .jpeg, for more information see documentation for the `ggsave()`
#'     function from `ggplot2`.
#' @param ... Additional arguments passed to forestplot.
#'
#' @importFrom scales scientific
#' @importFrom forestplot forestplot fpTxtGp fpColors fpDrawNormalCI fpDrawSummaryCI
#' @importFrom grid unit grid.grabExpr grid.draw gpar
#' @importFrom ggplot2 ggsave
#' @export
forestPlotModelComparison <- function(clinicalModelStats, stat, isSummary, filePath,
                                       fileName, ...) {

    indexClinical <- as.matrix(clinicalModelStats$clinical[[stat]])
    indexPCOSP <- as.matrix(clinicalModelStats$PCOSP[[stat]])

    spacing <- rep(NA, 4)
    plotData <- matrix(nrow=0, ncol=4)
    for (i in seq_len(nrow(indexClinical))) {
        plotData <-
            rbind(
                plotData,
                spacing,
                indexClinical[i, 1:4],
                indexPCOSP[i, 1:4],
            )
    }

    names <- c("Clinical model", "PCOSP", "")
    labels <- as.vector(vapply(rownames(indexPCOSP), function(name, names) c(name, names),
                               names=names,
                               FUN.VALUE=character(4)))

    summaries <- c(TRUE, unlist(lapply(isSummary,
                               function(is) {
                                   if (is) {
                                       c(rep(TRUE, 4))
                                       } else { c(rep(TRUE,1), rep(FALSE, 3))
                                           }
                                   })))


    isSeq <- clinicalModelStats$clinical$isSequencing

    # Construct the forest plot table
    labelText <- data.frame(
        "cohort"=c("Cohorts", labels),
        "pvalue"=c("P value", NA, c(scales::scientific(plotData[-1,][, "pval"], 2)))
    )

    plotData <- rbind(rep(NA, 4), plotData)

    # Extract box sizes
    lengthSeq <- length(unlist(clinicalModelStats$clinical$probabilities[isSeq]))
    lengthArray <- length(unlist(clinicalModelStats$clinical$probabilities[!isSeq]))

    # Match correct plot function to call
    if(missing(...)) {
        if (stat == "dIndex") {
            plot <- .forestPlotDindex2(labelText, plotData, summaries)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex2(labelText, plotData, summaries)
        } else {
            stop(paste0("There is no statistic called: ", stat))
        }
        # Allow user to specify custom parameters
    } else {
        if (stat == "dIndex") {
            plot <- .forestPlotDindex2(labelText, plotData,
                                       summaries, ...)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex2(labelText, plotData,
                                       summaries, ...)
        } else {
            stop(paste0("There is no statistic called: ", stat))
        }
    }
    # Decide whether to plot to device or save to disk
    if (missing(filePath) || missing(fileName)) {
        grid.draw(plot)
    } else {
        grid.draw(plot)
        ggsave(filename=fileName, path=filePath, plot=plot)
    }
}


## FIXME:: Refactor these into one function with more parameters!
.forestPlotCindex2 <- function(labelText, plotData, isSummary,
                              ...) {

    if (!missing(...)) {
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     plotData[, c("mean", "lower", "upper")],
                                                     is.summary=isSummary,
                                                     ...))
    } else {
        # Set plot colouring functions
        # ##TODO:: Determine if there is a more readable way to write this?
        fn <- local({
            i = 0

            b_clrs =  c("palevioletred1","darkgrey","palevioletred1","darkgrey", "palevioletred1","darkgrey","palevioletred1","darkgrey")
            l_clrs =    c("palevioletred1","darkgrey","palevioletred1","darkgrey", "palevioletred1","darkgrey","palevioletred1","darkgrey")
            function(..., clr.line, clr.marker){
                i <<- i + 1
                fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
            }
        })

        fn1 <- local({
            i = 0

            s_clrs =c("palevioletred1","darkgrey","palevioletred1","darkgrey","black","black")
            function(..., col){
                i <<- i + 1
                fpDrawSummaryCI(...,col=s_clrs[i])
            }
        })
        # Make the plot
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     plotData[, c("mean", "lower", "upper")],
                                                     xlab="Concordance index",
                                                     is.summary=isSummary,
                                                     clip=c(0,3.0),
                                                     cex=10,
                                                     fn.ci_norm=fn,
                                                     fn.ci_sum=fn1,
                                                     zero=0.5,
                                                     graphwidth=unit(2, "inches"),
                                                     align=c("l"),
                                                     new_page = FALSE,
                                                     txt_gp=fpTxtGp(label=gpar(fontfamily="Helvetica"),
                                                                    ticks=gpar(cex=0.8),
                                                     xlab =gpar(fontfamily = "Helvetica",
                                                                cex = 1)),
                                                     col = fpColors(text="black")))
    }
    return(plot)
}

#' @importFrom scales scientific
#' @importFrom forestplot forestplot fpTxtGp fpColors fpDrawNormalCI fpDrawSummaryCI
#' @importFrom grid unit grid.grabExpr grid.draw gpar
#' @importFrom ggplot2 ggsave
.forestPlotDindex2 <- function(labelText, plotData, isSummary, ...) {
    if (!missing(...)) {
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     plotData[, c("mean", "lower", "upper")],
                                                     is.summary=isSummary,
                                                     ...))
    } else {
        # Set plot colouring functions
        # ##TODO:: Determine if there is a more readable way to write this?
        fn <- local({
            i = 0
            b_clrs=c("palevioletred1","darkgrey", "palevioletred1","darkgrey",
                     "palevioletred1","darkgrey", "palevioletred1","darkgrey")
            l_clrs=c("palevioletred1", "darkgrey", "palevioletred1", "darkgrey",
                     "palevioletred1","darkgrey", "palevioletred1","darkgrey")
            function(..., clr.line, clr.marker){
                i <<- i + 1
                fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
            }
        })

        fn1 <- local({
            i = 0
            s_clrs =c("palevioletred1", "darkgrey", "palevioletred1", "darkgrey",
                      "black", "black")
            function(..., col){
                i <<- i + 1
                fpDrawSummaryCI(...,col=s_clrs[i])
            }
        })
        # Make the plot
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     plotData[, c("mean", "lower", "upper")],
                                                     xlab="Log2 D-index",
                                                     is.summary=isSummary,
                                                     clip=c(-2,3.0),
                                                     cex=9,
                                                     fn.ci_norm=fn,
                                                     fn.ci_sum=fn1,
                                                     zero=0,
                                                     graphwidth=unit(2, "inches"),
                                                     align=c("l"),
                                                     new_page = FALSE,
                                                     txt_gp = fpTxtGp(label=gpar(fontfamily = "Helvetica"),
                                                                      ticks=gpar(cex=0.8),
                                                     xlab=gpar(fontfamily="Helvetica",
                                                               cex=1)),
                                                     col = fpColors(text="black")))
    }
    return(plot)
}

#' Draw a forest plot using the statistics calcualted in validation metaestimate
#'
#' @examples
#'
#' @param validationStatsDF A \code{data.frame} containing statistics for Dindex
#'     or concordance index, as available in the \code{list} returned by
#'     `validateMetaEstimateCalculation`.
#' @param stat A \code{character} vector containing the statistic to
#'     use in the forest plot. Options are "dIndex" or "cIndex".
#' @param names A \code{character} vector of names for each classifier in the plot.
#'     Assumes the names are in the same order as the items in `classifierStats`.
#' @param filePath A \code{character} vector containing the file path to write
#'     the plotting results to.
#' @param fileName A \code{character} vector containing the desired file name
#'     with the desired file extension. Supported extensions include .pdf, .png,
#'     .svg and .jpeg, for more information see documentation for the `ggsave()`
#'     function from `ggplot2`.
#' @param ... Additional arguments passed to forestplot.
#'
#' @importFrom scales scientific
#' @importFrom forestplot forestplot fpTxtGp fpColors fpDrawNormalCI fpDrawSummaryCI
#' @importFrom grid unit grid.grabExpr grid.draw gpar
#' @importFrom ggplot2 ggsave
#' @export
forestPlotClassifierModelComparision <- function(classifierStats, stat, names, filePath,
                                       fileName, ...) {

    names(classifierStats) <- names

    microarray <- lapply(classifierStats, function(class) class[[stat]]["Microarray", ])
    microarrayMat <- as.matrix(do.call(rbind, microarray))

    sequencing <- lapply(classifierStats, function(class) class[[stat]]["Sequencing", ])
    seqMat <- as.matrix(do.call(rbind, sequencing))

    overall <- lapply(classifierStats, function(class) class[[stat]]["Overall", ])
    allMat <- as.matrix(do.call(rbind, overall))

    ## FIXME:: Generalize to N classifiers
    plotData <-
        rbind(
            rep(NA, 4),
            "Microarray Cohorts"=rep(NA, 4),
            microarrayMat,
            rep(NA, 4),
            "Sequencing Cohorts"=rep(NA, 4),
            seqMat,
            rep(NA, 4),
            "Overall"=rep(NA, 4),
            allMat
        )

    summaries <- rep(TRUE, nrow(plotData))

    isSeq <- classifierStats[[1]]$isSequencing

    # Construct the forest plot table
    labelText <- data.frame(
        "cohort"=c("Classifier", rownames(plotData)),
        "pvalue"=c("P value", NA, c(scales::scientific(plotData[-1,][, "pval"], 2)))
    )

    plotData <- rbind(rep(NA, 4), plotData)

    # Match correct plot function to call
    if(missing(...)) {
        if (stat == "dIndex") {
            plot <- .forestPlotDindex3(labelText, plotData, summaries)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex3(labelText, plotData, summaries)
        } else {
            stop(paste0("There is no statistic called: ", stat))
        }
        # Allow user to specify custom parameters
    } else {
        if (stat == "dIndex") {
            plot <- .forestPlotDindex3(labelText, plotData,
                                       summaries, ...)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex3(labelText, plotData,
                                       summaries, ...)
        } else {
            stop(paste0("There is no statistic called: ", stat))
        }
    }
    # Decide whether to plot to device or save to disk
    if (missing(filePath) || missing(fileName)) {
        grid.draw(plot)
    } else {
        grid.draw(plot)
        ggsave(filename=fileName, path=filePath, plot=plot)
    }
}


## FIXME:: Refactor these into one function with more parameters!
#' @importFrom scales scientific
#' @importFrom forestplot forestplot fpTxtGp fpColors fpDrawNormalCI fpDrawSummaryCI
#' @importFrom grid unit grid.grabExpr grid.draw gpar
#' @importFrom ggplot2 ggsave
.forestPlotCindex3 <- function(labelText, plotData, isSummary,
                               ...) {

    if (!missing(...)) {
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     plotData[, c("mean", "lower", "upper")],
                                                     is.summary=isSummary,
                                                     ...))
    } else {
        # Set plot colouring functions
        # ##TODO:: Determine if there is a more readable way to write this?
        fn1 <- local({
            i = 0

            s_clrs =c(c("#666666","#666666","#666666","#E7298A"),
                      c("#666666","#666666","#666666","#E7298A"),
                      c("#666666","#666666","#666666","#E7298A"))
            function(..., col){
                i <<- i + 1
                fpDrawSummaryCI(...,col=s_clrs[i])
            }
        })
        # Make the plot
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     plotData[, c("mean", "lower", "upper")],
                                                     xlab="Concordance index",
                                                     is.summary=isSummary,
                                                     new_page=FALSE,
                                                     fn.ci_sum = fn1,
                                                     clip=c(0.4,0.9),
                                                     txt_gp=fpTxtGp(label=gpar(fontfamily = "Helvetica"),
                                                                    ticks=gpar(cex=0.8),
                                                                    xlab =gpar(fontfamily="Helvetica",
                                                                               cex=1)),
                                                     col=fpColors(text="black"),
                                                     title="",
                                                     zero=0.5,
                                                     graphwidth=unit(2, "inches"),
                                                     align=c("l"),
                                                     boxsize = 0.25))
    }
    return(plot)
}

#' @importFrom scales scientific
#' @importFrom forestplot forestplot fpTxtGp fpColors fpDrawNormalCI fpDrawSummaryCI
#' @importFrom grid unit grid.grabExpr grid.draw gpar
#' @importFrom ggplot2 ggsave
.forestPlotDindex3 <- function(labelText, plotData, isSummary, ...) {
    if (!missing(...)) {
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     plotData[, c("mean", "lower", "upper")],
                                                     is.summary=isSummary,
                                                     ...))
    } else {
        # Set plot colouring functions
        # ##TODO:: Determine if there is a more readable way to write this?
        fn1 <- local({
            i = 0

            s_clrs =c(c("#666666","#666666","#666666","#E7298A"),
                      c("#666666","#666666","#666666","#E7298A"),
                      c("#666666","#666666","#666666","#E7298A"))
            function(..., col){
                i <<- i + 1
                fpDrawSummaryCI(...,col=s_clrs[i])
            }
        })
        # Make the plot
        plot <- grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     plotData[, c("mean", "lower", "upper")],
                                                     is.summary=isSummary,
                                                     xlab="Log2 D-index",
                                                     new_page=FALSE,
                                                     clip=c(-1,4),
                                                     txt_gp=fpTxtGp(label=gpar(fontfamily = "Helvetica"),
                                                                    ticks=gpar(cex=0.8),
                                                                    xlab=gpar(fontfamily="Helvetica",
                                                                              cex=1)),
                                                     col = fpColors(text="black"),
                                                     title=" ",
                                                     zero=0,
                                                     graphwidth=unit(2, "inches"),
                                                     align=c("l"),
                                                     fn.ci_sum=fn1,
                                                     boxsize = 0.25))
    }
    return(plot)
}

#' Draw a forest plot using the statistics calcualted in validation metaestimate
#'
#'
#' @param validationStatsDF A \code{data.frame} containing statistics for Dindex
#'     or concordance index, as available in the \code{list} returned by
#'     `validateMetaEstimateCalculation`.
#' @param stat A \code{character} vector containing the statistic to
#'     use in the forest plot. Options are "dIndex" or "cIndex".
#' @param names A \code{character} vector of names for each classifier in the plot.
#'     Assumes the names are in the same order as the items in `classifierStats`.
#' @param filePath A \code{character} vector containing the file path to write
#'     the plotting results to.
#' @param fileName A \code{character} vector containing the desired file name
#'     with the desired file extension. Supported extensions include .pdf, .png,
#'     .svg and .jpeg, for more information see documentation for the `ggsave()`
#'     function from `ggplot2`.
#' @param ... Additional arguments passed to forestplot.
#'
#' @importFrom scales scientific
#' @importFrom forestplot forestplot fpTxtGp fpColors fpDrawNormalCI fpDrawSummaryCI
#' @importFrom grid unit grid.grabExpr grid.draw gpar
#' @importFrom ggplot2 ggsave
#' @export
forestPlotCohortSubtypeComparison <- function(cohortSubtypeStats, stat, cohortNames, summaryNames, filePath,
                                                 fileName, ...) {

    formattedSubtypeStats <- .formatCohortSubtypeStatsForPlot(cohortSubtypeStats, stat)

    # Get per cohort Data
    overall <- as.data.frame(formattedSubtypeStats$Overall)
    cohortPlotData <- lapply(overall,
                             function(column, isBasal, isClassic)
                                 rbind(
                                     rep(NA, 4),
                                     column[!(isClassic | isBasal)],
                                     column[isBasal],
                                     column[isClassic]
                                 ),
                             isBasal=grepl("basal", rownames(overall)),
                             isClassic=grepl("classical", rownames(overall)))

    cohortPlotMat <- rbind(
        do.call(rbind, cohortPlotData[names(cohortPlotData) != "combined"]),
        rep(NA, 4))
    rownames(cohortPlotMat) <- cohortNames

    # Overall summary
    overallSummary <- cohortPlotData$combined

    # Sequencing summary
    sequencing <- as.data.frame(formattedSubtypeStats$Sequencing)
    sequencingSummary <- lapply(sequencing,
                              function(column, isBasal, isClassic)
                                  rbind(
                                      rep(NA, 4),
                                      column[!(isClassic | isBasal)],
                                      column[isBasal],
                                      column[isClassic]
                                  ),
                              isBasal=grepl("basal", rownames(overall)),
                              isClassic=grepl("classical", rownames(overall)))$combined

    # Microarray summary
    array <- as.data.frame(formattedSubtypeStats$Mircoarray)
    arraySummary <- lapply(array,
                          function(column, isBasal, isClassic)
                              rbind(
                                  rep(NA, 4),
                                  column[!(isClassic | isBasal)],
                                  column[isBasal],
                                  column[isClassic]
                              ),
                          isBasal=grepl("basal", rownames(overall)),
                          isClassic=grepl("classical", rownames(overall)))$combined

    # Combine summaries
    summaryData <- rbind(
        sequencingSummary,
        arraySummary,
        overallSummary
    )
    rownames(summaryData) <- summaryNames

    ## FIXME:: Generalize to N classifiers
    plotData <-
        rbind(
            cohortPlotMat,
            summaryData
        )
    colnames(plotData) <- c("mean", "lower", "upper", "pval")

    # Construct the forest plot table
    labelText <- data.frame(
        "cohort"=c("Cohorts", rownames(plotData)),
        "pvalue"=c("P value", NA, c(scales::scientific(plotData[-1,][, "pval"], 2)))
    )

    plotData <- rbind(rep(NA, 4), plotData)
    summaries <- c(TRUE, TRUE, rep(FALSE,40),rep(TRUE,13))

    # Match correct plot function to call
    if(missing(...)) {
        if (stat == "dIndex") {
            plot <- .forestPlotDindex4(labelText, plotData, summaries)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex4(labelText, plotData, summaries)
        } else {
            stop(paste0("There is no statistic called: ", stat))
        }
        # Allow user to specify custom parameters
    } else {
        if (stat == "dIndex") {
            plot <- .forestPlotDindex4(labelText, plotData,
                                       summaries, ...)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex4(labelText, plotData,
                                       summaries, ...)
        } else {
            stop(paste0("There is no statistic called: ", stat))
        }
    }
    # Decide whether to plot to device or save to disk
    if (missing(filePath) || missing(fileName)) {
        grid.draw(plot)
    } else {
        grid.draw(plot)
        ggsave(filename=fileName, path=filePath, plot=plot)
    }
}




## FIXME:: Refactor these into one function with more parameters!
#' @importFrom scales scientific
#' @importFrom forestplot forestplot fpTxtGp fpColors fpDrawNormalCI fpDrawSummaryCI
#' @importFrom grid unit grid.grabExpr grid.draw gpar
#' @importFrom ggplot2 ggsave
.forestPlotCindex4 <- function(labelText, plotData, isSummary,
                               ...) {

    if (!missing(...)) {
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     plotData[, c("mean", "lower", "upper")],
                                                     is.summary=isSummary,
                                                     ...))
    } else {
        # Set plot colouring functions
        # ##TODO:: Determine if there is a more readable way to write this?
        fn <- local({
            i = 0

            l_clrs = c(c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#E7298A"))
            b_clrs = c(c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#E7298A"))

            function(..., clr.line, clr.marker){
                i <<- i + 1
                fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
            }
        })
        fn1 <- local({
            i = 0

            s_clrs =c("#666666", "#1B9E77", "#E7298A", "#666666", "#1B9E77", "#E7298A", "#666666", "#1B9E77", "#E7298A")
            function(..., col){
                i <<- i + 1
                fpDrawSummaryCI(..., col=s_clrs[i])
            }
        })
        # Make the plot
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     plotData[, c("mean", "lower", "upper")],
                                                     xlab="Concordance index",
                                                     is.summary=isSummary,
                                                     new_page = FALSE,
                                                     clip=c(0.3, 0.8),
                                                     txt_gp = fpTxtGp(label = gpar(fontfamily = "Helvetica"),
                                                                      ticks = gpar(cex=0.5),
                                                                      xlab  = gpar(fontfamily = "Helvetica", cex = 0.5)),
                                                     col = fpColors(text="black"),
                                                     title="",
                                                     zero=0.5,
                                                     graphwidth=unit(2, "inches"),
                                                     align=c("l"),
                                                     hrzl_lines = list("2" = gpar(lty=2),
                                                                       "43" = gpar(lty=2, col = "#000044")),
                                                     vertices= TRUE,
                                                     fn.ci_norm = fn,
                                                     fn.ci_sum = fn1))
        }
    return(plot)
}

#' @importFrom scales scientific
#' @importFrom forestplot forestplot fpTxtGp fpColors fpDrawNormalCI fpDrawSummaryCI
#' @importFrom grid unit grid.grabExpr grid.draw gpar
#' @importFrom ggplot2 ggsave
.forestPlotDindex4 <- function(labelText, plotData, isSummary, ...) {
    if (!missing(...)) {
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     plotData[, c("mean", "lower", "upper")],
                                                     is.summary=isSummary,
                                                     ...))
    } else {
        # Set plot colouring functions
        # ##TODO:: Determine if there is a more readable way to write this?
        fn <- local({
            i = 0
            l_clrs = c(c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#E7298A"))
            b_clrs = c(c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#1B9E77","#E7298A"),
                       c("#666666","#1B9E77","#E7298A"), c("#666666","#E7298A"))
            function(..., clr.line, clr.marker){
                i <<- i + 1
                fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
            }
        })
        fn1 <- local({
            i = 0

            s_clrs =c("#666666","#1B9E77","#E7298A","#666666","#1B9E77","#E7298A","#666666","#1B9E77","#E7298A")
            function(..., col){
                i <<- i + 1
                fpDrawSummaryCI(...,col=s_clrs[i])
            }
        })
        # Make the plot
        plot <-grid::grid.grabExpr(forestplot::forestplot(labelText,
                                                     plotData[, c("mean", "lower", "upper")],
                                                     xlab="Log2 D-index",
                                                     is.summary=isSummary,
                                                    new_page=FALSE,
                                                    clip=c(-1, 2.5),
                                                    txt_gp = fpTxtGp(label = gpar(fontfamily = "Helvetica"),
                                                                     ticks = gpar(cex=0.8),
                                                                     xlab  = gpar(fontfamily = "Helvetica")),
                                                    col = fpColors(summary ="blue",
                                                                   text="black"),
                                                    title=" ",
                                                    zero=0,
                                                    graphwidth=unit(2, "inches"),
                                                    align=c("l"),
                                                    hrzl_lines = list("2" = gpar(lty=2),
                                                                      "43" = gpar(lty=2, col = "#000044")),
                                                    vertices= TRUE,
                                                    fn.ci_norm = fn,
                                                    fn.ci_sum = fn1))
    }
    return(plot)
}
