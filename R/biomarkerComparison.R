#'
#'
#'
#'
#'
#' @export
computeGeneBiomarkerScores <- function(cohortsDataL, sampleMetaClassDT, geneName) {
  cohortsDataL <- lapply(cohortsDataL, na.omit)
  subsetExpressionL <- lapply(cohortsDataL, `[`, i=geneName, j=TRUE)
  subsetExpressionL <- normalizeCohortsList(subsetExpressionL)
  sampleMetaClassDT <- na.omit(sampleMetaClassDT)[!duplicated(samples), ]
  subsetCohortExpression <- unlist(subsetExpressionL)
  names(subsetCohortExpression) <- unlist(lapply(cohortsDataL, colnames))
  keep <- intersect(names(subsetCohortExpression), sampleMetaClassDT$samples)
  keepDT <- sampleMetaClassDT[samples %in% keep, ]
  keepDT$expression <- subsetCohortExpression[keep]
  annotateSampleMetaClassDT(keepDT, c("Basal", "Exocrine", "Classical"))
}

#'
#'
#'
boxplotBiomarkerScores <- function(geneBiomarkerScores, comparisons, geneName, saveDir, fileName) {
  plot <- ggboxplot(geneBiomarkerScores, x="metaClasses", y = "expression",
                    color = "metaClasses", palette = "jco",
                    add = "jitter", xlab = "", ylab=geneName, 
                    legend =NULL) +
          stat_compare_means(label="p.format", comparisons=comparisons)
  
  if (!missing(saveDir) && !missing(fileName)) {
    ggsave(file.path(saveDir, fileName), plot)
    return(plot)
  } else {
    return(plot)
  }
}

#'
#'
#'
#' 
#'
boxplotBiomarkerScoresList <- function(geneBiomarkerScoreL, biomarkersOfInterest, comparisons) {
  plotL <- mapply(boxplotBiomarkerScores, 
                  geneBiomarkerScores=geneBiomarkerScoreL,
                  geneName=biomarkersOfInterest,
                  MoreArgs=list(comparisons=comparisons),
                  SIMPLIFY=FALSE)
  names(plotL) <- names(biomarkersOfInterest)
  return(plotL)
}

