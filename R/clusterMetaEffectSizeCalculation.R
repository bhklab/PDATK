#'
#'
#'
#'
#'
#'
#'
makeGenewiseCohortsDT <- function(preprocCohortDataL, annotSampMetaClassDT) {
  cohortNames <- names(preprocCohorts)
  cohortT <- lapply(preprocCohorts, t)
  cohortRNs <- lapply(cohortT, rownames)
  cohortLengths <- lapply(cohortRNs, length)
  cohDT <- lapply(cohortT, as.data.table)
  cohortNames <- mapply(rep, x=cohortNames, times=cohortLengths)
  mergeCohortDT <- rbindlist(cohDT)
  mergeCohortDT[, `:=`(cohorts=unlist(..cohortNames), samples=unlist(..cohortRNs))]
  meltCohortDT <- melt(mergeCohortDT, id.vars=c('cohorts', 'samples'),
                       value.name="expr", variable.name="gene")
  genewiseCohortDT <- merge(meltCohortDT, annotSampMetaClassDT,
                            by=c('cohorts', 'samples'))
  return(genewiseCohortDT)
}

#'
#'
#'
#'
#'
#'
#'
calcClustMetaEstStats <- function(gwCohortsDT) {
  # Split into list of all cohort-gene combinations
  splitDT <- split(gwCohortsDT, by=c('cohorts', 'gene'))

  message("Preprocessing the data...")
  for (DT in splitDT) {
    DT[, `:=`(respBasal=as.numeric(metaClasses == "Basal"),
              respClassical=as.numeric(metaClasses == "Classical"),
              respExocrine=as.numeric(metaClasses == "Exocrine"))]
  }

  message("Calculating cohen.d estimates")
  # For each class, calculate the gene-wise cohen.d the predicted metaclass
  # (i.e., all samples for that gene and class)
  cohenEstimateL <- list()
  for (class in c("Basal", "Classical", "Exocrine")) {
    cohenEstimateL[[class]] <-
      lapply(splitDT, function(DT, class) {
        colName <- paste0("resp",  class)
        cohen.d(DT$expr, DT[[colName]], hedge.correct=TRUE)$estimate
      }, class=class)
  }

  message("Calculating gene-wise standard deviation...")
  # For each class, calculate the gene-wise standard deviation
  # (i.e., all samples of for that gene and class)
  geneSD <- list()
  for (class in c("Basal", "Classical", "Exocrine")) {
    geneSD[[class]] <-
      lapply(splitDT, function(DT, class) {
        sd(DT[metaClasses == class, ]$expr)
      }, class=class)
  }

  message("Reorganizing results...")
  # Flattent he lists
  cohL <- unlist(cohenEstimateL)
  wtL <- 1 / unlist(geneSD)

  # Reform the lists grouped on class(i.e, L[[class]])
  cohLClass <- list()
  wtClass <- list()
  for (class in c("Basal", "Classical", "Exocrine")) {
    cohLClass[[class]] <- cohL[grepl(class, names(cohL))]
    wtClass[[class]] <- wtL[grepl(class, names(geneSD))]
  }

  # Reform the class lists grouped on gene (i.e, L[[class]][[gene]][allCohortExpr])
  clGnCoh <- list()
  clGenWt <- list()
  for (class in c("Basal", "Classical", "Exocrine")) {
    for (gene in unique(gwCohortsDT$gene)) {
      clGnCoh[[class]][[gene]] <-
        cohLClass[[class]][grepl(gene, names(cohLClass[[class]]))]
      clGenWt[[class]][[gene]] <-
        wtClass[[class]][grepl(gene, names(wtClass[[class]]))]
    }
  }


  message("Computing weighted means...")
  # Add na.rm for
  weighted.mean2 <- function (x, w, ..., na.rm = FALSE)
  {
    if (missing(w)) {
      if (na.rm)
        x <- x[!is.na(x)]
      return(sum(x)/length(x))
    }
    if (length(w) != length(x))
      stop("'x' and 'w' must have the same length")
    if (na.rm) {
      i <- !is.na(x) & !is.na(w)
      w <- w[i]
      x <- x[i]
    }
    sum((x * w)[w != 0])/sum(w)
  }

  weightedMean <- list()
  for (class in c("Basal", "Classical", "Exocrine")) {
    weightedMean[[class]] <- mapply(weighted.mean2,
                                    x=clGnCoh[[class]], clGenWt[[class]],
                                    MoreArgs=list(na.rm=TRUE))
  }

  weightedMeans <- as.data.table(do.call(cbind, weightedMean), keep.rownames="gene")
  return(weightedMeans)
}