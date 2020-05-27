#'
#'
#'
#'
computeCohortCommunityMatrix <- function(cohortsDataL, metaGenes) {
  # Subset on metaGenes and convert to matrix
  cohortsMetaGenesL <- lapply(cohortsDataL, `[`, i=metaGenes, j=TRUE)
  cohortsMetaGenesML <- lapply(cohortsMetaGenesL, as.matrix)
  # Impute missing values; Suppress cat/print output
  invisible(capture.output(cohortsImputeL <- lapply(cohortsMetaGenesML, impute.knn)))
  cohortsMetaGenesL <- lapply(cohortsImputeL, `[[`, 'data')
  # Merge matrixes colwise
  cohortCommunityMatrix <- do.call(cbind, cohortsMetaGenesL)
  return(cohortCommunityMatrix)
}

#' Train a K-TSP classifier based on cohort sample meta-classes
#'
#' @param cohortComMat
#' @param sampleMetaClassDT
#' @param maxK
#' 
#' @export
trainSingleSampleClassifer <- function(cohortComMat, sampleMetaClassDT, class, maxK) {
  sampleMetaClassDT <- sampleMetaClassDT[!duplicated(samples), ]
  keepSamples <- intersect(colnames(cohortComMat), sampleMetaClassDT$samples)
  cohortComMat <- cohortComMat[, keepSamples]
  metaClasses <- sampleMetaClassDT[samples %in% keepSamples, ]$metaClasses
  phenoGroup <- as.factor(as.numeric(metaClasses == class))
  SWAP.KTSP.Train(cohortComMat, phenoGroup, maxK)
}

#'
#'
#'
#'
calcTopGenes1Gt2Matrix <- function(cohortTopGenesMat, topGenes1, topGenes2, genePairs) {
  topGenes1Gt2 <- apply(cohortTopGenesMat, 2, 
                        function(col, genes1, genes2) 
                          list(as.numeric(col[genes1] > col[genes2])),
                        genes1=unlist(topGenes1),
                        genes2=unlist(topGenes2))
  topGenes1Gt2 <- lapply(topGenes1Gt2, unlist)
  g1Gt2Matrix <- matrix(unlist(topGenes1Gt2), 
                        nrow=ncol(cohortTopGenesMat), 
                        byrow=TRUE)
  colnames(g1Gt2Matrix) <- genePairs
  rownames(g1Gt2Matrix) <- colnames(cohortTopGenesMat)
  return(t(g1Gt2Matrix))
}

#'
#'
#'
#'
#'
predictSingleSampleClasses <- function(g1Gt2Matrix, metaClassesFactor, nthread) {
  opts <- options()
  options("mc.cores"=nthread)
  on.exit(options(opts))
  
  sampleClassPredictions <-
    bplapply(seq_len(ncol(g1Gt2Matrix)),
             function(i, trainMat, trainLabels) {
               print(i)
               model <- randomForest(trainMat[-i, ], trainLabels[-i])
               print(paste("Finished model", i))
               prediction <- predict(model, trainMat[i, ])
               rm(model); gc()
               print(paste('Finished prediction', i))
               return(as.character(prediction))
             },
            trainMat=t(g1Gt2Matrix),
            trainLabels=metaClassesFactor)
  
  names(sampleClassPredictions) <- colnames(g1Gt2Matrix)
  return(sampleClassPredictions)
}