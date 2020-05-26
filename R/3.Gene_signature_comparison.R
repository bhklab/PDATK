########################################################################################################
### Required libraries
#########################################################################################################
library(ggpubr)
library(gridExtra)

#' Get the per cohort sample metaclusters and return as a data.table
#'
#' @param annotatedClusters
#'
#' @export
extractSampleMetaClasses <- function(annotatedClusters) {
  metaClasses <- lapply(annotatedClusters, `[[`, "metaClasses")
  sampleNames <- unlist(lapply(metaClasses, names))
  classLengths <- lapply(metaClasses, length)
  cohorts <- mapply(rep, names(annotatedClusters), classLengths, SIMPLIFY=FALSE)
  data.table("cohorts"=unlist(cohorts), "samples"=sampleNames, 
             "metaClasses"=unlist(metaClasses))
}

#'
#' @param signatureGenes A \code{character} vector of gene names in the
#'    expression signature.
#' @param signatureName A \code{character} vector containing the name of the
#'    gene signature.
#' @param comparison A \code{list} of character vectors containing the names
#'    of cohorts to compare. All comparisons must be pairwise (i.e., each
#'    character vector can have only two names).
#'
#' @export
computeSignatureScoreDT <- function(cohortsDataL, sampleMetaClassDT, signatureGenes) {
    cohortsSigGenes <- lapply(cohortsDataL, `[`, i=signatureGenes, j=TRUE)
    cohortsSigGenes <- lapply(cohortsSigGenes, na.omit)
    .normalize <- function(cohort) t(scale(t(cohort)))
    normalizedCohorts <- lapply(cohortsSigGenes, .normalize)
    geneScoreList <- lapply(normalizedCohorts, colMeans)
    sampleNames <- unlist(lapply(geneScoreList, names))
    geneScoreDT <- data.table("samples"=sampleNames, 
                              "sigScores"=unlist(geneScoreList))
    sampleClassDT <- annotateSampleMetaClassDT(sampleMetaClassDT, 
                                         c('Basal', 'Exocrine', 'Classical'))
    sigScoreDT <- merge(geneScoreDT[!duplicated(samples)], 
                        sampleClassDT[!duplicated(samples)], 
                        on="samples")
    
}

#'
#'
#'
#'
annotateSampleMetaClassDT <- function(sampleMetaClassDT, clusterLabels) {
  DT <- sampleMetaClassDT
  DT[, metaClasses := as.character(DT$metaClasses)]
  i <- 1
  for (val in na.omit(unique(DT$metaClasses))) {
    DT[metaClasses == val, metaClasses := clusterLabels[i]]
    i <- i + 1
  }
  return(na.omit(DT))
}


#' Draw a boxplot of the signature scores per sample per meta-class
#'
#' @param sigScore
#'
#' @exprot
plotSigScores <- function(cohortSignatureScoreDT, comparisons, signatureScoreName, saveDir, fileName) {
  plot <- ggboxplot(cohortSignatureScoreDT, x = "metaClasses", y = "sigScores",
                    color="metaClasses",
                    palette="jco",
                    ylab=signatureScoreName,
                    add="jitter", xlab="", legend="none") + 
          stat_compare_means(comparisons=comparisons)
  if (!missing(saveDir) && !missing(fileName)) {
    ggsave(plot, file=file.path(saveDir, fileName))
    return(plot)
  } else {
    return(plot)
  }
}
