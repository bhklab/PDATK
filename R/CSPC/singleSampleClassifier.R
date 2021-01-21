#' Subset all cohort expression data to `metaGenes`, cbind them and return a gene by
#'    sample expression matrix with missing values imputed using
#'    `impute::impute.knn`.
#'
#' @param cohortsDataL A \code{list} per cohort gene by sample gene expression
#'   matrixes.
#' @param metaGenes A \code{character} vector of meta genes to subset on, as
#'    returned by `getMetaGenes`.
#'
#' @return A \code{matrix} of gene by sample expression data, with the rows
#'    containing `metaGenes` and the columns containing sample from all cohorts.
#'
#' @importFrom impute impute.knn
#' @export
computeCohortCommMat <- function(cohortsDataL, metaGenes) {
  # Subset on metaGenes and convert to matrix
  cohortsMetaGenesL <- lapply(cohortsDataL, `[`, i=metaGenes, j=TRUE)
  cohortsMetaGenesML <- lapply(cohortsMetaGenesL, as.matrix)
  # Impute missing values; Suppress cat/print output
  invisible(capture.output(cohortsImputeL <- lapply(cohortsMetaGenesML, impute.knn)))
  cohortsMetaGenesL <- lapply(cohortsImputeL, `[[`, '../../data')
  # Merge matrixes colwise
  cohortCommunityMatrix <- do.call(cbind, cohortsMetaGenesL)
  return(cohortCommunityMatrix)
}

#' Train a K-TSP classifier based on cohort sample meta-classes
#'
#' @param cohortComMat A gene by sample \code{matrix} containing all cohort
#'   samples. As returned by the `computeCohortCommMat` function from this package.
#' @param sampleMetaClassDT A \code{data.table} containing the per cohort per
#'   sample meta-class predictions.
#' @param class A \code{character} vector specifying the meta-class to train
#'     the single sample classifer for.
#' @param maxK A \code{numeric} vector with the integer maximum K to use in
#'      `switchBox::SWAP.KTSP.Train`.
#'
#' @return A \code{list} of KTSP classifiers, one for each meta-class.
#'
#' @importFrom switchBox SWAP.KTSP.Train
#' @export
trainSingleSampleClassifer <- function(cohortComMat, sampleMetaClassDT, class, maxK) {
  sampleMetaClassDT <- na.omit(sampleMetaClassDT[!duplicated(samples), ])
  keepSamples <- intersect(colnames(cohortComMat), sampleMetaClassDT$samples)
  cohortComMat <- cohortComMat[, keepSamples]
  metaClasses <- sampleMetaClassDT[samples %in% keepSamples, ]$metaClasses
  phenoGroup <- as.factor(as.numeric(metaClasses == class))
  SWAP.KTSP.Train(cohortComMat, phenoGroup, maxK)
}

#' Calculate top two genes for each biomarker
#'
#' @param cohortTopGenesMat A \code{matrix} of normalized gene expression
#'    values for the genes in `topGenes1` and `topGenes2`. Subset from the
#'    matrix returned by `computeCohortCommMat` function in this pacakge.
#' @param topGenes1 A \code{character} vector of genes names for the top scoring
#'    pairs item predicted to have higher expression on average.
#' @param topGenes2 A \code{character} vector of gene names for the top scoring
#'    gene pairs item predicted to have lower expression on avareage.
#' @param genePairs A \code{character} vector with `topGenes1` and `topGenes`
#'    pasted together with '>'. This is used for label the gene pairs and can
#'    be any character vector of length equal to `topGenes`.`
#'
#' @return A gene-pair by sample binary \code{matrix} where 1 represents gene1
#'     is greater than gene2 (accurate prediction) and 0 represents gene2 less
#'     than gene1 (an inaccurate prediction).
#'
#' @export
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

#' Predict sample metaclasses using a trained single sample classifier
#'
#' @param g1Gt2Matrix A gene-pair by sample binary \code{matrix} where 1 represents gene1
#'     is greater than gene2 (accurate prediction) and 0 represents gene2 less
#'     than gene1 (an inaccurate prediction).
#' @param metaClassesFactor A \code{factor} vector with the meta-classes for
#'     each sample in g1Gt2Matrix.
#' @param nthread A \code{numeric} vector with the integer number of threads to
#'     parallelize over.
#'
#' @return A \code{list} of per sample metaclass predictions
#'
#' @importFrom randomForest randomForest
#' @export
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