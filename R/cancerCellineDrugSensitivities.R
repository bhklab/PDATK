#'
#'
#'
#'
#'
#'
formatCCLsensitivityLtoDT <- function(sensitivityList) {
  CCLsensitivityDtL <- lapply(CCLsensitivityL, t)
  .addSampNames <- function(mat) {DF <- as.data.frame(mat)
                                  DF$sample <- row.names(mat)
                                  return(DF)}
  CCLsensitivityDtL <- lapply(CCLsensitivityDtL, .addSampNames)
  CCLsensitivityDtL <- lapply(CCLsensitivityDtL, as.data.table)
  return(CCLsensitivityDtL)
}


#'
#'
#'
#'
#'
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
  fClassPredL <- lapply(classPredL, function(DT) copy(DT)[, sample := gsub('\\.| ', '-', sample)])
  fDrugSensL <- lapply(drugSensL, function(DT) copy(DT)[, sample := gsub("\\.| ", "-", sample)])
  
  datasetNames <- names(fClassPredL)
  
  # Define utilit functions
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
                      function(DT) colnames(DT[, .SD, .SDcols=-c("sample", "dataset")]))
  
  # Merge the DTs on cellline and datset
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
  rbindlist(meltMergedByDataset)
}

#'
#'
#'
#'
#'
computeWilcoxOnSharedDrugs <- function(mergedDT, class1, class0) {
  ## TODO:: Can I make this faster than O(n^2)?
  mergedDT <- mergedDT[predClass %in% c(class1, class0)]
  splitOnDataset <- split(mergedDT, by="dataset")
  datasetDrugs <- lapply(splitOnDataset, `[[`, "drug")
  allDSdrugs <- Reduce(intersect, datasetDrugs)
  
  subsetSharedDrugs <- lapply(splitOnDataset, 
                              function(DT, allDSdrugs) subset(DT, drug %in% allDSdrugs),
                              allDSdrugs=allDSdrugs)
  rm(splitOnDataset); gc()
  
  datasetL <- list()
  i <- 1
  for (DT in subsetSharedDrugs) {
    message(paste0("Computing for:", names(subsetSharedDrugs)[i]))
    datasetPvals <- list()
    for (drg in allDSdrugs) {
      DT1 <- DT[drug == drg, ]
      datasetPvals[[drg]] <- wilcox.test(DT1$AUC ~ as.numeric(DT1$predClass == class1))
    }
    datasetL[[i]] <- simplify2array(datasetPvals)
    i <- i + 1
  }
  names(datasetL) <- names(subsetSharedDrugs)
  rm(subsetSharedDrugs); gc()
  return(datasetL)
}

#'
#'
#'
#'
#'
computeConcIndOnSharedDrugs <- function(mergedDT, class, na.rm=TRUE, method="noether") {
  ## TODO:: Can I make this faster than O(n^2)?
  splitOnDataset <- split(mergedDT, by="dataset")
  datasetDrugs <- lapply(splitOnDataset, `[[`, "drug")
  allDSdrugs <- Reduce(intersect, datasetDrugs)
  
  subsetSharedDrugs <- lapply(splitOnDataset, 
                              function(DT, allDSdrugs) subset(DT, drug %in% allDSdrugs),
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

#'
#'
#'
#'
#'
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

#'
#'
#'
#'
#'
#'
boxplotAUCperSubtypePerDataset <- function(mergedDT, conIdx) {
  splitOnDataset <- split(mergedDT, by="dataset")
  allDSdrugs <- colnames(conIdx[[1]])
  
  subsetSharedDrugs <- lapply(splitOnDataset, 
                              function(DT, allDSdrugs) subset(DT, drug %in% allDSdrugs),
                              allDSdrugs=allDSdrugs)
  rm(splitOnDataset); gc()
  
  DTnames <- names(subsetSharedDrugs)
  .boxplot <- function(DT, drg, CI, nm) {
    ci <- scientific(CI["c.index", ][[drg]], 2)
    p <- scientific(CI["p.value", ][[drg]], 2)
    title <- paste0(nm, "\n", drg, " CI: ", ci, ", P: ", p)
    p <- ggboxplot(DT[drug == drg], x="predClass", y="drug", color="predClass", add="jitter",
                   xlab="Subtype", ylab="AUC",
                   title=title, pallette="jco", legend="none")
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
