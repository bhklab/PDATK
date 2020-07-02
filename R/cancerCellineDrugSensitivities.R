#' Take a list of drug by cell-line drug sensitivity data.frames and merge
#'    them into a single `data.table``
#'
#' @param CCLsensitivityL A \code{list} of per dataset drug by cell-line
#'   drug sensitivity `data.frame`s
#'
#' @return A \code{list} of `data.table`s annotated with drugs as columns
#'   and celline names in the `sample` column.
#'
#' @import data.table
#' @export
formatCCLsensitivityLtoDT <- function(CCLsensitivityL) {
  CCLsensitivityDtL <- lapply(CCLsensitivityL, t)
  .addCellNames <- function(mat) {DF <- as.data.frame(mat)
                                  DF$sample <- row.names(mat)
                                  return(DF)}
  CCLsensitivityDtL <- lapply(CCLsensitivityDtL, .addCellNames)
  CCLsensitivityDtL <- lapply(CCLsensitivityDtL, as.data.table)
  return(CCLsensitivityDtL)
}


## FIXME:: Are the drug sensitivity AUC or AAC? Both are interchaned
## in the original script.
#' Merge per dataset sensitivity data with per dataset predicted cell-line
#'    subtype/meta-class.
#'
#' @param classPredL A \code{list} of per cohort cell-line meta-class predictions,
#'    such as a list of output from the `predictSampleMetaClass` function in
#'    this package.
#' @param drugSensL A \code{list} of per dataset cell-line drug sensitivity
#'    `data.table`s as returned by the `formatCCLsensitvityLtoDT` function
#'    in this package.
#'
#' @return A \code{data.table} containing the columns 'celline', 'dataset',
#'    'predClass', 'pBasal', 'pClassical', 'pExocrine', 'drug' and 'AUC'.
#'
#' @import data.table
#' @export
mergeClassAndDrugs <- function(classPredL, drugSensL) {
  if (!all(names(classPredL) == names(drugSensL))) {
    drugSensL <- drugSensL[names(classPredL)]
    if (any(is.na(names(drugSensL)))) stop("Tried to reorder class and
                                           drug lists to same dataset names,
                                           but something went wrong.\nDo both
                                           lists have the same datasets?\n
                                           Are the list names equal?")
  }
  # Format sample names between class and drugs
  fClassPredL <- lapply(classPredL,
                        function(DT)
                           copy(DT)[, sample := gsub('\\.| ', '-', sample)])
  fDrugSensL <- lapply(drugSensL,
                       function(DT)
                           copy(DT)[, sample := gsub("\\.| ", "-", sample)])

  datasetNames <- names(fClassPredL)

  # Define utility functions
  .addDSnames <- function(DT, dsName) {DT$dataset <- rep(dsName, nrow(DT)); DT}

  fClassPredL <- mapply(.addDSnames,
                       DT=fClassPredL,
                       dsName=datasetNames,
                       SIMPLIFY=FALSE)
  fDrugSensL <- mapply(.addDSnames,
                       DT=fDrugSensL,
                       dsName=datasetNames,
                       SIMPLIFY=FALSE)
  # Get drugs names for metling
  drugNames <- lapply(fDrugSensL,
                      function(DT)
                          colnames(DT[, .SD, .SDcols=-c("sample", "dataset")]))

  # Merge the DTs on cellline and dataset
  mergedByDataset <- mapply(function(DT1, DT2) merge(DT1, DT2, by=c("sample", "dataset")),
                            DT1=fClassPredL,
                            DT2=fDrugSensL,
                            SIMPLIFY=FALSE)

  # Make drugs and AUCs into their own columns
  meltMergedByDataset <- mapply(melt,
                                data=mergedByDataset, measure.vars=drugNames,
                                MoreArgs=list(variable.name="drug", value.name="AUC",
                                              variable.factor=FALSE),
                                SIMPLIFY=FALSE)

  # Unify into a single DT
  DT <- rbindlist(meltMergedByDataset)

  # Rename sample column to cell line
  colnames(DT) <- replace(colnames(DT), which(colnames(DT) == "sample"), "cellline")

  return(DT)
}

#' Calculate the two-sided wilcox p-values for a drug shared between two
#'    meta-classes/subtypes for all drugs common to all included datasets.
#'
#' @param mergedDT A \code{data.table} with merged cell-line class predictions
#'   and cell-line drug sensitivity, as returned by the `mergeClassandDrugs`
#'   function in this paackage. Note: You may want to subset out datasets
#'   which do not have many overlapping drugs between them, as only common drugs
#'   will be included in the test.
#' @param class1 A \code{character} vector with the class/subtype label
#'   to compare against.
#' @param class0 A \code{character} vector with the class/subtype label
#'   to compare with.
#'
#' @return A \code{lost} of p-values for a two-sided wilcox test between
#'   cell-lines for drugs common to all included datasets.
#'
#' @import data.table
#' @export
computeWilcoxOnSharedDrugs <- function(mergedDT, class1, class0) {
  mergedDT <- mergedDT[predClass %in% c(class1, class0)]
  splitOnDataset <- split(mergedDT, by="dataset")
  datasetDrugs <- lapply(splitOnDataset, `[[`, "drug")
  allDSdrugs <- Reduce(intersect, datasetDrugs)

  subsetSharedDrugs <- lapply(splitOnDataset,
                              function(DT, allDSdrugs)
                                  subset(DT, drug %in% allDSdrugs),
                              allDSdrugs=allDSdrugs)
  rm(splitOnDataset); gc()

  datasetL <- list()
  i <- 1
  for (DT in subsetSharedDrugs) {
    message(paste0("Computing for: ", names(subsetSharedDrugs)[i]))
    datasetPvals <- list()
    for (drg in allDSdrugs) {
      DT1 <- DT[drug == drg, ]
      tryCatch({datasetPvals[[drg]] <- wilcox.test(x=DT1[predClass==class1]$AUC, y=DT1[predClass==class0]$AUC, paired=FALSE)},
               error=function(e) {
                 message(paste0('   Caught error:\n     ', e, "   Returning NA for drug: ", drg, "!\n"))
                 datasetPvals[[drg]] <- NA
               })
    }
    if (!all(is.na(datasetPvals)))
      datasetL[[i]] <- simplify2array(datasetPvals)
    else
      datasetL[[i]] <- NA
    i <- i + 1
  }
  names(datasetL) <- names(subsetSharedDrugs)
  rm(subsetSharedDrugs); gc()
  return(datasetL)
}

#' Compute concordance index between all datasets for a specific predicted
#'     subtype/meta-class.
#'
#' @param mergedDT A \code{data.table} with merged cell-line class predictions
#'   and cell-line drug sensitivity, as returned by the `mergeClassandDrugs`
#'   function in this paackage.
#' @param class A \code{character} vecotr specifying the meta-class/subtype
#'     to compute concordance index for.
#' @param na.rm A \code{boolean} specifying whether to
#'     remove NA's from the data when calculating concordance index.
#' @param method A \code{character} vector specifying the method to use
#'     for calculating concordance index. Passed to `survcomp::concordance.index`
#'     as the `method` argument.
#'
#' @importFrom survcomp concordance.index
#' @import data.table
#' @export
computeConcIndOnSharedDrugs <- function(mergedDT, class, na.rm=TRUE, method="noether") {
  ## TODO:: Can I make this faster than O(n^2)?
  splitOnDataset <- split(mergedDT, by="dataset")
  datasetDrugs <- lapply(splitOnDataset, `[[`, "drug")
  allDSdrugs <- Reduce(intersect, datasetDrugs)

  subsetSharedDrugs <- lapply(splitOnDataset,
                              function(DT, allDSdrugs)
                                    subset(DT, drug %in% allDSdrugs),
                              allDSdrugs=allDSdrugs)
  rm(splitOnDataset); gc()

  datasetL <- list()
  i <- 1
  for (DT in subsetSharedDrugs) {
    message(paste0("Computing for:", names(subsetSharedDrugs)[i]))
    datasetCI <- list()
    for (drg in allDSdrugs) {
      DT1 <- DT[drug == drg, ]
      datasetCI[[drg]] <-
        concordance.index(DT1$AUC,
                          cl=as.numeric(DT1$predClass == class),
                          na.rm=na.rm,
                          method=method)
    }
    datasetCI <- simplify2array(datasetCI)
    datasetL[[i]] <- datasetCI
    i <- i + 1
  }
  names(datasetL) <- names(subsetSharedDrugs)
  rm(subsetSharedDrugs); gc()
  return(datasetL)
}

#' Calculate meta concordance index and meta p-values from a list of per
#'   dataset concordance index results.
#'
#' @param concordanceInds A \code{list} of per dataset concordance index
#'   results for a specfic meta-class.
#'
#' @return A \code{data.table} with the columns 'drugs', 'metaCI', and 'metaPval'
#'
#' @importFrom survcomp combine.est
#' @export
calculateMetaStats <- function(concordanceInds) {
    drugs <- colnames(concordanceInds[[1]])

    metaStats <- lapply(drugs, function(drg, ci)
                            combine.est(
                              x=unlist(lapply(ci, `[`, i='c.index', j=drg)),
                              x.se=unlist(lapply(ci, `[`, i='se', j=drg)),
                              hetero=TRUE,
                              na.rm=TRUE
                            ),
                        ci=concordanceInds)

    metaCI <- lapply(metaStats, `[[`, "estimate")

    .metaP <- function(s) 2*pnorm(s$estimate - 0.5/s$se, lower.tail=s$c.index < 0.5)
    metaPval <- lapply(metaStats, .metaP)

    data.table("drugs"=drugs, "metaCI"=metaCI, "metaPval"=metaPval)
}

#' Boxplot the distribution of AUC values for each subtype in each dataset
#'
#' @param mergedDT A \code{data.table}
#' @param conIdx A \code{numeric}
#' @param paletee A \code{character} vector specifying the RColorBrewer palette
#'     to use for the plot. Passed to `ggpubr::ggboxplot`s `palette` argument.
#'
#' @return A \code{list} object
#'
#' @importFrom ggpubr ggboxplot
#' @import scales
#' @import data.table
#' @export
boxplotAUCperSubtypePerDataset <- function(mergedDT, conIdx, palette="Set1") {
  splitOnDataset <- split(mergedDT, by="dataset")
  allDSdrugs <- colnames(conIdx[[1]])

  subsetSharedDrugs <- lapply(splitOnDataset,
                              function(DT, allDSdrugs) subset(DT, drug %in% allDSdrugs),
                              allDSdrugs=allDSdrugs)
  rm(splitOnDataset); gc()

  DTnames <- names(subsetSharedDrugs)
  .boxplot <- function(DT, drg, CI, nm) {
    ci <- scales::scientific(CI["c.index", ][[drg]], 2)
    p <- scales::scientific(CI["p.value", ][[drg]], 2)
    title <- paste0(nm, "\n", drg, " CI: ", ci, ", P: ", p)
    p <- ggboxplot(DT[drug == drg], x="predClass", y="drug", color="predClass", add="jitter",
                   xlab="Subtype", ylab="AUC",
                   title=title, pallette=palette, legend="none")
    return(p)
    }

  plots <- list()
  i <- 1
  for (drug in allDSdrugs) {
    plots[[drug]] <- mapply(.boxplot, DT=subsetSharedDrugs, CI=conIdx,
                            MoreArgs=list(drg=drug, nm=DTnames[i]),
                            SIMPLIFY=FALSE)
    i <- i + 1
  }
  plots <- Reduce(c, plots)
  names(plots) <- paste(rep(names(subsetSharedDrugs), length(allDSdrugs)), allDSdrugs, sep="_")
  return(plots)
}