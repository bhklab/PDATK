#' Density plot AUC distribution of a PCOSP model
#'
#' Create a density plot of the distrubtion of AUCs between the overall,
#'    sequencing and microarray cohorts in a PCOSP model.
#'
#' @param formattedValCohorts A \code{list} of formatted validation cohorts
#' @param selectedModels A \code{list} of kTSP classificaiton models.
#' @param seqCohorts A \code{character} vector of names of cohorts
#'     containing sequencing data.
#' @param vlines A \code{numeric} vector of length 3 indicating
#'     where to draw a vertical line for the combine, sequencing
#'     and microarray density plots, respectively.
#' @param nthread A \code{numeric} vector indicating the integer
#'     number of threads to parallelize over.
#' @param title A \code{character} vector with the title for the
#'     plot
#' @param filePath A \code{character} vector with the path to
#'     save the plot to. If absent, the plot is just returned,
#'     not saved to disk. Passed to `ggsave` function.
#' @param fileName A \code{character} vector with the file name
#'     and image file extension to save the plot under. Passed
#'     to `ggsave` function.
#'
#' @return A \code{grob} object from \code{ggplot2}. Also saves to disk if `filePath`
#'     and `fileName` are specified.
#'
#' @importFrom ggplot2 geom_density geom_vline labs ggtitle facet_wrap theme ggsave
#' @export
densityPlotModel <- function(formattedValCohorts, selectedModels, seqCohorts, title,
                             vlines, nthread, filePath, fileName) {

    dStats <- .densityStats(formattedValCohorts, selectedModels, seqCohorts,
                            nthread)

    densityDF <- data.frame(
        "platforms"=rep(c("1. Array-based", "2. Sequencing", "3. Overall"), each=length(selectedModels)),
        "BAC"=unlist(c(dStats$arrayAUCs, dStats$seqAUCs, dStats$metaAUCs))
    )
    vlineDF <- data.frame(
        "platforms"=unique(densityDF$platforms), v1=vlines
    )

    plt <- ggplot(densityDF, aes(x=BAC, color=platforms)) +
        geom_density(aes(fill=platforms), alpha=0.5) +
        geom_vline(data=vlineDF,
                   aes(xintercept=v1,
                       color=platforms),
                   linetype="dashed",
                   size=1.5) +
        labs(x = "Balanced Accuracy", y="Density")+
        ggtitle(title) +
        facet_wrap(~ platforms, ncol=1, strip.position="right") +
        theme(plot.title=element_text(hjust = 0.5),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "grey"),
              legend.position="none")

    if (!missing(filePath) && !missing(fileName)) {
        ggsave(filename=fileName, path=filePath, plot=plt)
    }
    plt
}


#' Calcualte the density plot statistics for a list of validation
#'    cohorts and trainedd kTSP classificaiton models
#'
#' @param formattedValCohorts A \code{list} of formatted validation cohorts
#' @param selectedModels A \code{list} of kTSP classificaiton models.
#' @param seqCohorts A \code{character} vector of names of cohorts
#'     containing sequencing data.
#' @param nthread A \code{numeric} vector indicating the integer
#'
#' @keywords interal
.densityStats <- function(formattedValCohorts, selectedModels, seqCohorts,
                          nthread) {
    KTSPs <- predictKTSPs(formattedValCohorts, selectedModels, nthread)

    isSeq <- names(KTSPs) %in% seqCohorts

    nModels <- length(selectedModels)

    metaAUCs <- .summarizeKTSPs(KTSPs, nModels, nthread)
    seqAUCs <- .summarizeKTSPs(KTSPs[isSeq], nModels, nthread)
    arrayAUCs <- .summarizeKTSPs(KTSPs[!isSeq], nModels, nthread)

    return(list(
        "metaAUCs"=metaAUCs,
        "seqAUCs"=seqAUCs,
        "arrayAUCs"=arrayAUCs
    ))

}

#' Predict kTSP classification for a list of validation cohorts and
#'    classifier models
#'
#' @param formattedValCohorts A \code{list} of formatted validation cohorts
#' @param selectedModels A \code{list} of kTSP classificaiton models.
#' @param nthread A \code{numeric} vector with the integer number of
#'     thread to parallelize over.
#'
#' @return A \code{list} predicted kTSP classifications
#'
#' @export
predictKTSPs <- function(formattedValCohorts, selectedModels, nthread) {
    lapply(formattedValCohorts,
           function(cohort, models, nthread) .predictKTSP(cohort, models, nthread),
           models=selectedModels,
           nthread=nthread)
}

#' Calculate the `combine.est` for a list of kTSP model predictions
#'    repeatedly
#'
## FIXME:: Is this documentation correct?
#' @param KTSPs A \code{list} kTSP classificaiton predictions
#'    for AUC and AUC standard error.
#' @param nModels A \code{numeric} vector indicating the integer
#'    number of KTSPs
#'
#' @importFrom BiocParallel bplapply
#' @importFrom survcomp combine.est
.summarizeKTSPs <- function(KTSPs, nModels, nthread) {
    # Temporily change number of cores to parallelize over
    opts <- options()
    options("mc.cores"=nthread)
    on.exit(options(opts))

    AUCs <- do.call(rbind, lapply(KTSPs, function(KTSP) KTSP$AUCs))
    aucSEs <- do.call(rbind, lapply(KTSPs, function(KTSP) KTSP$aucSEs))

    bplapply(seq_len(nModels),
             function(i, AUCs, aucSEs)
                 combine.est(
                     AUCs[, i],
                     aucSEs[, i],
                     na.rm=TRUE
                 )$estimate,
             AUCs=AUCs,
             aucSEs=aucSEs)
}