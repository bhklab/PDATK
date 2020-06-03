#' Use an existing classifier centroid to classify exprData
#'
#' @param exprData A \code{matrix} of gene by cell-line gene expression values.
#' @param centroid A \code{matrix} of gene by subtype cluster centroids from
#'    an existing classifer.
#' @param seed A \code{numeric} vector with the integer seed to set for
#'    consensus clustering in this function.
#' @param reps A \code{numeric} vector with the integer number of times to
#'    repeat consensus clustering.
#'
#' @return A \code{data.table} with the predicted subtype for each classifier
#'   in each cell-line.
#'
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @import data.table
#' @export
subtypeWithClassifier <- function(exprData, centroid, seed=NULL, reps=100) {
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

    ## FIXME:: Original COINCIDE clustering used Pearson distance?
    exprDataClust <- consensusClusterCohort(cohort=exprData,
                                            maxK=4,
                                            distance="spearman",
                                            method="hc",
                                            reps=reps,
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
                                            gsub('^[^-]*-', '', clustCors$pair))

    return(annotSampClustLabs[, .(sample, metaClusters, pair, pairCor)])
}


#' Predict the meta class for a list of datasets using a list of classifiers
#'
#' @param rawDataL A \code{list} of `data.table`s containing the raw expression
#'    data for a set of sample, such as for an experimental cohort.
#' @param classifCentroidL A \code{list} of classifier centroids to
#' @param predMetaClassL A \code{list} of predicted meta-classes using the COINCIDE
#'    classifier in this pacakge. Each list item should correspond to the dataset
#'    in `rawDataL`.
#' @param seed An optional \code{numeric} vector containing the random seed
#'    to use when consensus clustering data using the supplied classifier list.
#'    If a seed was used for the initial COINCIDE classification, the same seed
#'    should be used here to ensure comparability of the results.
#'
#' @return A \code{data.table} containing the classification results
#'
#' @import data.table
#' @export
subtypeDataLwClassifCentroidL <- function(rawDataL, classifCentroidL,
                                          predMetaClassL, seed=NULL) {

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
        metaDT <- rbindlist(predMetaClassL)[, .SD, .SDcols=colnames(subtypeDT)]
        subtypeDT <- rbind(subtypeDT, metaDT)
    }
    return(subtypeDT)
}


#' Calculate the correlation between the centroids from two classifiers
#'
#' @param centroid1 A \code{data.table}, \code{data.frame} or \code{matrix}
#'     containing the centroids from a classifier
#' @param centroid2 A \code{data.table}, \code{data.frame} or \code{matrix}
#'     containing the centroids from a classifier
#'
#' @return A \code{data.table} with the correlation the cluster centroids
#'     in each clustering result.
#'
#' @import data.table
#' @export
correlateCentroids <- function(centroid1, centroid2) {
    clustCors <- cor(centroid1, centroid2)
    maxRowCors <- apply(clustCors, 1, which.max)
    maxColCors <- apply(clustCors, 2, which.max)[maxRowCors]
    clustMatches <- data.table("pair"=paste(names(maxRowCors), names(maxColCors), sep="-"),
                               "cent1idx"=maxColCors, "cent2idx"=maxRowCors,
                               "pairCor"=apply(clustCors, 1, max))
    return(clustMatches)
}


#' Plot a heatmap showing the classification identify of each cell-line using
#'   each classifier.
#'
#' @param pubClassifSubtypeDT \code{data.table} with the predicted subtype for each classifier
#'   in each cell-line.
#' @param allClassif A \code{boolean} indicating whether plot only cell-line
#'   with sybtypes in all classifiers.
#' @param saveDir An optional \code{character} vector specifying the path
#'    to the directory where the plot should be saved. If excluded, fileName
#'    will not work.
#' @param fileName An optional \code{character} vector specifying the
#'    name and extension of the file to save the plot it. This is passed to
#'    the `ggplot2::ggsave`.
#'
#' @return A \code{ggplot} object showing the classification of each cell-line
#'   in each classifier.
#'
#' @export
plotClassifierComparisons <- function(pubClassifSubtypeDT, allClassif=TRUE, saveDir, fileName) {
    if (allClassif) {
        splitOnClassifL <- split(pubClassifSubtypeDT, by='classif')
        sharedSamples <- Reduce(intersect, lapply(splitOnClassifL, `[[`, "sample"))
        DT <- copy(pubClassifSubtypeDT[sample %in% sharedSamples, ])
    } else {
        DT <- copy(pubClassifSubtypeDT)
    }

    DT[, classifSubtype := mapply(paste, classif, metaClusters, sep=":")]

    plot <- ggplot(DT, aes(x=factor(classif, levels=unique(classif)), y=sample, fill=classifSubtype)) +
        geom_tile(color='black') +
        labs(y="Cell-line", x="Classifier", fill="Subtype")

    if(!missing(saveDir) && !missing(fileName)) {
        ggsave(file.path(saveDir, fileName), plot)
        message(paste0("Saved to ", file.path(saveDir, fileName)))
    }
    return(plot)
}

#' Calculate the correlation betwen
#'
#' @param pubClassifSubtypeDT \code{data.table} with the predicted subtype for each classifier
#'   in each cell-line.
#'
#' @return A \code{data.table} with the cramers V and p.value for each pair-wise
#'   classifier comparison.
#'
#' @importFrom vcd assocstats
#' @import data.table
#' @export
calcAssocStats <- function(pubClassifSubtypeDT) {
    # Get all samples classified with metaClassifier
    samples <- pubClassifSubtypeDT[classif == "metaClass",]$sample

    # Make a copy to prevent modifying source by reference
    DT <- copy(pubClassifSubtypeDT)

    # Get cartesian product of classifiers
    classifComparisons <- as.data.table(DT[, expand.grid(unique(classif), unique(classif),
                                                         stringsAsFactors=FALSE)])

    # Find the best pair for each sample (by pair correlation)
    maxCors <- DT[sample %in% samples, max(pairCor), by=.(classif, sample, cohort)]
    colnames(maxCors)[4] <- "pairCor"
    DT <- merge(DT, maxCors, by=c("classif", "sample", "pairCor", "cohort"))

    # Append the classifier name to the predicted subtypes
    DT[, classifSubtype := mapply(paste, classif, metaClusters, sep=":")]

    # Count the number of classical and metaclassical classes per sample per dataset
    DT[, `:=`(classical=sum(classifSubtype == "metaClass:Classical"),
              basal=sum(classifSubtype == "metaClass:Basal")),
       by=.(cohort, sample)]

    # Count the total times each sample was classified by metalcass for all cohorts
    metaClassCounts <- DT[, list("basal"=sum(basal), "classical"=sum(classical)),
                          by=classifSubtype]

    #
    DT <- dcast(DT, sample + cohort ~ classif, value.var="classifSubtype")

    # Calculate association stats
    assocStats <- vector("list", nrow(classifComparisons))
    for (i in seq_len(nrow(classifComparisons))) {
        assocStats[[i]] <- summary(assocstats(table(DT[[classifComparisons[i, ]$Var1]],
                                                    DT[[classifComparisons[i, ]$Var2]])))
    }

    names(assocStats) <- mapply(paste,
                                classifComparisons$Var1, classifComparisons$Var2,
                                sep='-')

    pvals <- vapply(assocStats, function(stats) stats$summary$p.value,
                    FUN.VALUE=numeric(1))
    cramersVs <- vapply(assocStats, function(stats) stats$object$cramer,
                       FUN.VALUE=numeric(1))

    assocStatsDT <- data.table("comparison"=names(assocStats),
                               "classif1"=classifComparisons$Var1,
                               "classif2"=classifComparisons$Var2,
                               "cramersV"=cramersVs,
                               "pval"=pvals)

    return(assocStatsDT)
}

#' Heatmap the correlation between classifier subtypes
#'
#' @param assocStatsDT
#' @param dendro A \code{boolean} specifying whether or not to add a dendrogram
#'   to the right side of the plot.
#' @param saveDir An optional \code{character} vector specifying the path
#'    to the directory where the plot should be saved. If excluded, fileName
#'    will not work.
#' @param fileName An optional \code{character} vector specifying the
#'    name and extension of the file to save the plot it. This is passed to
#'    the `ggplot2::ggsave`.
#'
#' @return A \code{ggplot} object with a heatmap showing the correlation
#'   between subtypes between each classifier.
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom ggplotify as.ggplot
#' @importFrom ggdendro dendro_data segment
#' @import data.table
#' @export
heatmapClassifCors <- function(assocStatsDT, dendro=TRUE, saveDir, fileName) {

    # Define custom theme
    .noneTheme <- ggplot2::theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(colour=NA),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
    )

    # Set the upper triangle of the cramer stats to NAs
    cramerMat <- as.matrix(dcast(assocStatsDT, classif1 ~ classif2,
                                 value.var="cramersV")[, -'classif1'])

    # Build classifier dendrogram
    classifDendro <- as.dendrogram(hclust(dist(t(cramerMat))))
    sortClassifDendro <- order.dendrogram(classifDendro)
    classifDendroData <- ggdendro::dendro_data(classifDendro)
    dendroPlot <- ggplot(ggdendro::segment(classifDendroData)) +
                    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
                    coord_flip() + .noneTheme

    # Reorder the classifiers in the cramer mat and convert to DT
    cramerMat <- cramerMat[sortClassifDendro, sortClassifDendro]
    cramerMat[upper.tri(cramerMat)] <- 0
    cramerMatDT <- as.data.table(cramerMat)[, classif1 := colnames(cramerMat)]
    cramerDT <- melt(cramerMatDT, id.vars="classif1", variable.name="classif2",
                     value.name="cramersV")

    # Plot the correlation matrix heatmap with the dendrogram
    corHeatmap <-
        ggplot(cramerDT, aes(x=factor(classif1, levels=colnames(cramerMat)),
                                      y=factor(classif2, levels=colnames(cramerMat)),
                             fill=cramersV)) +
            geom_tile(color = "white")+
            scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                 midpoint = 0, limit = c(-1,1), space = "Lab",
                                 name="Cramer's V index") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                             size = 12, hjust = 1)) +
            coord_fixed() +
            theme(
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                panel.grid.major = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.ticks = element_blank(),
                legend.justification = c(1, 0),
                legend.position = c(0.5, 0.8),
                legend.direction = "horizontal") +
            guides(fill = guide_colorbar(barwidth = 10, barheight = 3,
                                         title.position = "top", title.hjust = 0.5))

    if (dendro) {
        spacingPlot1 <- ggplot() + theme_void()
        spacingPlot2 <- ggplot() + theme_void()

        dendro <- as.ggplot(grid.arrange(spacingPlot1, dendroPlot,
                                         spacingPlot2, ncol=1,
                                         heights=c(0.2, 1, 0.2)))
        plot <- as.ggplot(grid.arrange(corHeatmap, dendro, ncol=2,
                                       widths=c(1, 0.2)))
    } else {
        plot <- corHeatmap
    }

    if(!missing(saveDir) && !missing(fileName)) {
        ggsave(file.path(saveDir, fileName), plot)
        message(paste0("Saved to ", file.path(saveDir, fileName)))
    }
    return(plot)
}


##TODO:: Move below to utilities.R

#'
#' @param clusterData A \code{data.table}, \code{data.frame} or \code{matrix}
#'   containing the clustering data.
#' @param clusterLabels A \code{}
#'
#' @import data.table
#' @keywords internal
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
#' @param centroid1 A \code{matrix}
#' @param centroid2 A \code{matrix}
#'
#' @return A \code{character} vector with the names of each pairwise
#'   comparison possible between two classifier centroids.
#'
#' @keywords internal
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