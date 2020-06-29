#' Predict the metaclass and class probabilities for each sample
#'
#' @param exprTable A gene by sample matrix-like object containing normalized
#'    log2 expression values for each gene in each sample.
#' @param classifModel A trained \code{randomForest} classifier, as returned
#'    by the `randomForest::randomForest` function.
#' @param topGenesDT A \code{data.table} with the columns 'topGenes1'
#'    and 'topGenes2' containing the top scoring pairs, 'genePairs' containing
#'    a label for each pair, and 'classes' indicating the meta-class each gene pair
#'    is associated with.
#' @param trainData A gene-pair by sample binary \code{matrix} where 1 represents
#'    the predicted top scoring pair being correctly ordered, as returned
#'    by the `calcTopGenes1Gt2Matrix` function in this pacakge.
#' @param trainLabels A \code{factor} vector with the meta-class for each
#'    sample in the training data.
#'
#' @return A \code{data.table} with the columns 'sample', 'predClass', 'pBasal',
#'     'pClassical' and 'pExocrine' holding the predictions and probabilities
#'     for each sample's metaclass.
#'
#' @importFrom randomForest randomForest
#' @import data.table
#' @export
predictSampleMetaClass <- function(exprTable, classifModel, topGenesDT,
                                   trainData, trainLabels) {
  # Class labels
  classLabs <- levels(trainLabels)

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
  gn1Gt2Matrix <- exprTable[topGn1, ] > exprTable[topGn2, ]
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

#' Create a waterfall plot of tumour response (measured as change in volume)
#'
#' @param classSurvCompDT A \code{data.table} containing two meta-classes
#'    to compare in the waterfall plot. To tests differences between specific
#'    drug response, also subset to the drugs of interest.
#' @param noXaxis A \code{boolean} indicating wheter to hide x-axis labels and
#'    ticks.
#' @param pVal An optional \code{character} vector containing the p-value for
#'    difference in survival between the classes in the plot.
#' @param saveDir An optional \code{character} vector specifying the path
#'    to the directory where the plot should be saved. If excluded, fileName
#'    will not work.
#' @param fileName An optional \code{character} vector specifying the
#'    name and extension of the file to save the plot it. This is passed to
#'    the `ggplot2::ggsave`.
#'
#' @return A \code{ggplot} object containing the watefall plot
#'
#' @importFrom ggpubr ggbarplot
#' @importFrom ggplot2 ggsave annotate theme
#' @import data.table
#' @export
waterfallPlotTumorResponse<- function(classSurvCompDT, noXaxis=FALSE, pVal,
                                      saveDir, fileName) {
  plot <- ggbarplot(classSurvCompDT[!is.na(tumorResponse), ],
                    x="studyID", y="tumorResponse",
                    fill="predClass", color="predClass",
                    pallette="Set1", xlab="Study ID", ylab="% Change in Tumor Size",) +
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

#' Boxplot difference between in neutrophil-lymphocyte ratio between
#'   meta-classes, comparin means with `method`
#'
#' @param classSurvCompDT A \code{data.table} containing meta-class labels and
#'    per sample survival data.
#' @param method A \code{character} vector specifying the mean comparison method
#'    passed to `ggpubr::stat_compare_means`.
#' @param saveDir An optional \code{character} vector specifying the path
#'    to the directory where the plot should be saved. If excluded, fileName
#'    will not work.
#' @param fileName An optional \code{character} vector specifying the
#'    name and extension of the file to save the plot it. This is passed to
#'    the `ggplot2::ggsave`.
#'
#' @return A \code{gpplot} with predicted class on the x-axis and the
#'    neutrophil-lymphcyte ratio on the y-axis.
#'
#' @importFrom ggpubr ggboxplot stat_compare_means
#' @importFrom ggplot2 ggsave
#' @export
boxplotClassNLR <- function(classSurvCompDT, method="kruskal.test", palette="Set1",
                            saveDir, fileName) {
  plot <- ggboxplot(classSurvCompDT, x="predClass", y="NLR", color="predClass",
                    palette=palette, add="jitter") +
          stat_compare_means(method=method, label.x=0.6,
                             label.y=max(classSurvCompDT$NLR) * 1.1)

  if(!missing(saveDir) && !missing(fileName)) {
    ggsave(file.path(saveDir, fileName), plot)
    message(paste0("Saved to ", file.path(saveDir, fileName)))
  }
  return(plot)
}

#' Boxplot the log2 expression for each biomarker between meta-clases
#'
#' @param compassLogExprMat A gene by sample \code{matrix} of normalized
#'    log2 expression values.
#' @param sampClassPredwSurvivalDT A \code{data.table} of cohort meta-class
#'    predictions and survival data merged by sample.
#' @param biomarkers A \code{character} vector of biomarker gene names, with
#'    each gene named for the meta-cluster it is a biomarker in.
#' @param saveDir An optional \code{character} vector with the path to
#'    the directory in which the plot should be saved.
#' @param fileName An optional \code{character} vector specifying the name
#'    and extension of the file. Only used if `saveDir` is also specified.
#'    Saving is done via the `ggsave` function from `ggplot2`.
#'
#' @return A \code{gglot} object containing a plot grid with a boxplot
#'   for each specified biomarker.
#'
#' @importFrom cowplot plot_grid
#' @export
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

#' Boxplot expression comparisons for the specified biomarker and subtypes
#'
#' @param subtypes A \code{character} vector containing the subtype for each
#'    sample in the expression matrix.
#' @param biomarker A \code{character} vector specifying the biomarker
#'    to compare between subtypes
#' @param exprMat A \code{matrix} with X by Y gene expression data
#' @param palette A \code{character} vector specifying the name of the
#'    RColorBrewer palette to use for the plot. Passed to `ggpubr::ggboxplot`
#'    as the `palette` argument.
#'
#' @return A \code{ggplot} object containing a boxplot of log2 expression
#'    for the specified biomarker in each subtype.
#'
#' @importFrom ggpubr ggboxplot
#' @export
compassPlotBiomarker <- function(subtypes, biomarker, exprMat) {
  DF <- data.frame(subtypes, as.numeric(exprMat[biomarker, ]))
  colnames(DF) <- c("Subtype", biomarker)
  ggboxplot(DF, x="Subtype", y=colnames(DF)[2], color="Subtype", add="jitter",
            pallette="Set1", legend="none")
}