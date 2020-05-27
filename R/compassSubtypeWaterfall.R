#'
#'
#' @param
#' @param
#' @param
#'
#' @export
predictSampleMetaClass <- function(compassLogExprMat, binaryRFmodel, topGenesDT, classLabs=c("Basal", "Classical", "Exocrine")) {
  # Subset data
  topFeatures <- unlist(topGenesDT[, .(topGenes1, topGenes2)])
  topGn1 <- topGenesDT$topGenes1
  topGn2 <- topGenesDT$topGenes2
  compassExprFeatsM <- compassLogExprMat[topFeatures, ]

  # Compare log expression between topGn1 and topGn2
  gn1Gt2Matrix <- compassExprFeatsM[topGn1, ] > compassExprFeatsM[topGn2, ] * 1
  gn1Gt2Matrix <- t(gn1Gt2Matrix * 1) # Convert logical to numeric
  colnames(gn1Gt2Matrix) <- topGenesDT$genePairs
  
  # Make predictions
  binaryModelPredProbs <- predict(binaryRFmodel, gn1Gt2Matrix, type="prob")
  colnames(binaryModelPredProbs) <- paste0('p', classLabs)
  binaryModelPredClasses <- as.character(predict(binaryRFmodel, gn1Gt2Matrix))
  
  # Build data.table
  binaryModPredDT <- as.data.table(binaryModelPredProbs)
  binaryModPredDT$sample <- rownames(binaryModelPredProbs)
  binaryModPredDT$predClass <- binaryModelPredClasses
  
  # Name classes by subtypes
  for (i in seq_along(classLabs)) {
    binaryModPredDT[predClass == i, predClass := classLabs[i]]
  }
  
  return(binaryModPredDT[, .(sample, predClass, pBasal, pClassical, pExocrine)])
}

#'
#'
#'
#'
barplotTumorResponse<- function(classSurvCompDT, noXaxis=FALSE, pVal, saveDir, fileName) {
  plot <- ggbarplot(classSurvCompDT[!is.na(tumorResponse), ], 
                    x="studyID", y="tumorResponse",
                    fill="predClass", color="predClass",
                    pallette="jco", xlab="Study ID", ylab="% Change in Tumor Size",) + 
          theme(axis.text.x=element_text(angle=90, size=9, vjust=0))
  if (noXaxis) {
    plot <- plot + theme(axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank())
  }
  if (!missing(pVal)) {
    plot <- plot + annotate("text", x=5, y=120, 
                            label=paste("Significance", pVal))
  } 
  
         
  if (!missing(saveDir) && ! missing(fileName)) {
    ggsave(file.path(saveDir, fileName), plot)
  }
  return(plot)
}

#'
#'
#'
#'
boxplotClassNLR <- function(classSurvCompDT, method="kruskal.test", 
                            saveDir, fileName) {
  plot <- ggboxplot(classSurvCompDT, x="predClass", y="NLR", color="predClass",
                    palette="jco", add="jitter") + 
          stat_compare_means(method=method, label.x=0.6, 
                             label.y=max(classSurvCompDT$NLR) * 1.1)
  
  if(!missing(saveDir) && !missing(fileName)) {
    ggsave(plot, file.path(saveDir, fileName))
    message(paste0("Saved to ", file.path(saveDir, fileName)))
  }
  return(plot)
}

#'
#'
#'
#'
#'
compassPlotBiomarkers <- function(compassLogExprMat, sampClassPredwSurvivalDT, 
                                  biomarkers, saveDir, fileName) {
  keepSamples <- intersect(colnames(compassLogExprMat),
                           sampClassPredwSurvivalDT$sample)
  
  keepBiomarkers <- biomarkers[which(biomarkers %in% rownames(compassLogExprMat))]
  if (!all(biomarkers %in% keepBiomarkers)) {
    message(paste0("Excluded biomarkers not in expresion matrix:", 
                   paste0(setdiff(biomarkers, keepBiomarkers), collapse=", ")))
  }

  
  biomarkersMat <- compassLogExprMat[biomarkers, keepSamples]
  subtypes <- sampClassPredwSurvivalDT$predClass
  plots <- lapply(keepBiomarkers, compassPlotBiomarker, 
                  subtypes=subtypes, exprMat=biomarkersMat)
  
  plot <- plot_grid(plotlist=plots, ncol=ceiling(sqrt(length(plots))))
  if(!missing(saveDir) && !missing(fileName)) {
    ggsave(plot, file.path(saveDir, fileName))
    message(paste0("Saved to ", file.path(saveDir, fileName)))
  }
  plot
}

#'
#'
#'
#'
#'
#'
compassPlotBiomarker <- function(subtypes, biomarker, exprMat) {
  DF <- data.frame(subtypes, as.numeric(exprMat[biomarker, ]))
  colnames(DF) <- c("Subtype", biomarker)
  ggboxplot(DF, x="Subtype", y=colnames(DF)[2], color="Subtype", add="jitter", 
            pallette="jco", legend="none")
}