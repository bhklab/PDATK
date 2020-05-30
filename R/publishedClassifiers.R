#'
#'
#'
#'
#'
subtypeWithClassifier <- function(exprData, centroid, seed=NULL, repeats=100) {
    if (!(all(rownames(exprData) %in% rownames(centroid)))) {
        gnNotMapped <- setdiff(rownames(centroid), rownames(exprData))
        message(paste0(
            nrow(centroid) - length(gnNotMapped), " out of ", nrow(centroid),
            " genes were mapped.\n\nNot mapped: ", paste0(gnNotMapped, collapse=", ")))
    } else {
        message("All genes mapped.")
    }
    keepGenes <- rownames(exprData) %in% rownames(centroid)

    exprData <- .scaleGenewise(limma::avereps(exprData[keepGenes, ]),
                                     scale=TRUE, center=TRUE)

    exprDataClust <- consensusClusterCohort(cohort=exprData,
                                            maxK=4,
                                            distance="pearson",
                                            method="hc",
                                            reps=repeats,
                                            seed=seed)

    keepCentGenes <- intersect(rownames(exprDataClust$centroidClusters),
                           rownames(centroid))

    exprCentroid <- as.data.table(exprDataClust$centroidClusters,
                                  keep.rownames="gene")[gene %in% keepCentGenes]
    classifCentroid <- as.data.table(centroid,
                                     keep.rownames="gene")[gene %in% keepCentGenes, ]

    clustCors <- correlateCentroids(exprCentroid[, -'gene'],
                                    classifCentroid[, -'gene'])

    clustCors[, metaClusters := cent1idx]

    sampleClustLabs <- data.table("sample"=names(exprDataClust$classes),
                                  "metaClusters"=exprDataClust$classes)

    sampClustLabs <- merge(sampleClustLabs, clustCors, by="metaClusters",
                           allow.cartesian=TRUE)

    annotSampClustLabs <- .annotateClusters(sampClustLabs,
                                            gsub('.*-', '', clustCors$pair))

    return(annotSampClustLabs[, .(sample, metaClusters, pair, pairCor)])
}


#'
#'
#'
#'
subtypeDataLwClassifCentroidL <- function(rawDataL, classifCentroidL,
                                          predMetaClassL, seed=NULL, nthread) {
    ## Temporily change number of cores to parallelize over
    if (!missing(nthread)) {
        opts <- options()
        options("mc.cores"=nthread)
        on.exit(options(opts))
    }

    ## Match classification using consenus method with cluaster from different
    ## classifiers
    DtL <- lapply(rawDataL,
                  function(exprData, classifs, seed) {
                      classifSubtypes <-
                          lapply(classifs,
                                 function(centroid, exprData, seed) {
                                     DT <- subtypeWithClassifier(exprData=exprData,
                                                                centroid=centroid,
                                                                seed=seed)
                                 },
                                 exprData=exprData, seed=seed)

                      .nameByClassif <- function(x, y) { x[['classif']] <- rep(y, nrow(x)); x }
                      subList <- mapply(.nameByClassif,
                                        classifSubtypes,
                                        names(classifSubtypes),
                                        SIMPLIFY=FALSE)

                      DT <- rbindlist(subList)
                      DT
                  }, seed=seed, classifs=classifCentroidL)

    # Annotate with cohort name
    for (i in seq_along(DtL)) {
        DtL[[i]][, cohort := rep(names(DtL)[i], nrow(DtL[[i]]))]
    }

    # Make into datatable
    subtypeDT <- rbindlist(DtL)

    # Add metaClass predictions
    if (!missing(predMetaClassL)) {
        cohNames <- names(predMetaClassL)
        for (i in seq_along(predMetaClassL)) {
            nrow <- nrow(predMetaClassL[[i]])
            predMetaClassL[[i]][, `:=`(cohort=rep(cohNames[i], nrow),
                                       classif=rep("metaClass", nrow),
                                       metaClusters=predClass,
                                       pair=paste("meta", predClass, sep='-'),
                                       pairCor=rep(NA, nrow))]
        }
        browser()
        metaDT <- rbindlist(predMetaClassL)[, .SD, .SDcols=colnames(subtypeDT)]
        subtypeDT <- rbind(subtypeDT, metaDT)
    }
    return(subtypeDT)
}


#'
#'
#'
#'
#'
#'
correlateCentroids <- function(centroid1, centroid2) {
    clustCors <- cor(centroid1, centroid2)
    maxRowCors <- apply(clustCors, 1, which.max)
    maxColCors <- apply(clustCors, 2, which.max)[maxRowCors]
    clustMatches <- data.table("pair"=paste(names(maxRowCors), names(maxColCors), sep="-"),
                               "cent1idx"=maxColCors, "cent2idx"=maxRowCors,
                               "pairCor"=apply(clustCors, 1, max))
    return(clustMatches)
}


#'
#'
#'
#'
#'
.annotateClusters <- function(clusterData, clusterLabels) {
    if (!is.data.table(clusterData)) {
        DT <- as.data.table(clusterData, keep.rownames=TRUE)
    } else {
        DT <- copy(clusterData)
    }
    DT[, metaClusters := as.character(metaClusters)]
    i <- 1
    for (val in na.omit(unique(DT$metaClusters))) {
        DT[metaClusters == val, metaClusters := clusterLabels[i]]
        i <- i + 1
    }
    DT
}


#'
#'
#'
#'
#'
.findClusterCombinations <- function(centroid1, centroid2) {
    if (missing(centroid2)) centroid2 <- centroid1

    pairs <- expand.grid(x=seq_len(ncol(centroid1)),
                         y=seq_len(ncol(centroid2)))

    namePairs <- expand.grid(x=colnames(centroid1), y=colnames(centroid2),
                             stringsAsFactors=FALSE)

    # Remove self comparisons
    if (missing(centroid2)) {
        # Remove self comparisons
        pairs <- pairs[-which(pairs[, 1] == pairs[, 2]), ]
        namePairs <- namePairs[-which(namePairs[, 1] == namePairs[, 2]), ]
    }
    allPairs <- pairs
    allNames <- namePairs

    # Paste together pair names
    pairNames <- mapply(paste, allNames[, 1], allNames[, 2], MoreArgs=list(sep="-"))

    # Assign names as rownames to pair data.frame
    rownames(allPairs) <- pairNames

    return(allPairs)
}