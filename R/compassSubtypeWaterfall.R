#'
#'
#' @param
#' @param
#' @param
#'
#' @export
predictSampleMetaClass <- function(exprTable, classifModel, topGenesDT, 
                                   trainData, trainLabels,  classLabs=c("Basal", "Classical", "Exocrine")) {
  # Find matching features
  topFeatures <- unlist(topGenesDT[, .(topGenes1, topGenes2)])
  
  # Retrain on subset of original training data if test data is missing feats
  if (!all(rownames(exprTable) %in% topFeatures)) {
    # Subset to shared features
    keepFeats <- intersect(rownames(exprTable), topFeatures)
    topGenesDT <- topGenesDT[topGenes1 %in% keepFeats & topGenes2 %in% keepFeats, ]
    exprTable <- exprTable[keepFeats, ]
    
    # Get the rownames to keep in the training data
    genePairs <- topGenesDT$genePairs
    
    # Retrain the model
    trainGn1Gt2Matrix <- trainData[genePairs, ]
    classifModel <- randomForest(t(trainGn1Gt2Matrix), trainLabels)
  }
  
  # Extract features
  topGn1 <- topGenesDT$topGenes1
  topGn2 <- topGenesDT$topGenes2
  genePairs <- topGenesDT$genePairs
  
  # Compare log expression between topGn1 and topGn2
  gn1Gt2Matrix <- exprTable[topGn1, ] > exprTable[topGn2, ] * 1
  gn1Gt2Matrix <- t(gn1Gt2Matrix * 1) # Convert logical to numeric
  colnames(gn1Gt2Matrix) <- topGenesDT$genePairs
  
  # Make predictions
  binaryModelPredProbs <- predict(classifModel, gn1Gt2Matrix, type="prob")
  colnames(binaryModelPredProbs) <- paste0('p', classLabs)
  binaryModelPredClasses <- as.character(predict(classifModel, gn1Gt2Matrix))
  
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
waterfallPlotTumorResponse<- function(classSurvCompDT, noXaxis=FALSE, pVal, saveDir, fileName) {
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