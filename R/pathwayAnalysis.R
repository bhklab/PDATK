#'
#' @param
#'
#'
#' @export
rankMetaClassGenesByEffSize <- function(metaEstimateStatsDT) {
  meltDT <- melt(metaEstimateStatsDT, id.vars='gene', variable.name="class",
                 value.name="wMean")
  topGenes <- meltDT[order(-abs(wMean)), head(.SD, 1000), by=class]
  return(topGenes)
}

#'
#'
#'
#'
#'
loadGnSetGMTs <- function(gmtFldir) {
  pathwayGMTfl <- list.files(gmtFldir, pattern='.gmt', full.names=TRUE)
  pathwayL <- lapply(pathwayGMTfl, loadGSC)
  pathwayNames <- lapply(pathwayL, function(GSC) gsub('_.*$', '', GSC$addInfo[1, 1]))
  names(pathwayL) <- pathwayNames
  return(pathwayL)
}

#'
#'
#'
#'
#'
computePathwayScores <- function(rankedMetaClassGenes, pathwayL, referenceGenes) {
  splitByClass <- split(rankedMetaClassGenes, by='class')

  pathStats <- list()
  pvalsDT <- list()
  for (path in pathwayL) {
    pathwayName <- gsub('_.*$', '', path$addInfo[1, 1])
     pathwayScores <-
        lapply(splitByClass,
           function(DT, refGns, path)
               runGSAhyper(genes=DT$gene, gsc=path, adjMethod="fdr", universe=refGns),
           refGns=referenceGenes,
           path=path)
    pvals <- lapply(pathwayScores, `[[`, "pvalues")
    classScores <- lapply(pathwayScores,
                          function(path) qnorm(1 - (path$pvalues/2)))

    pathwayScoresDT <-
        as.data.table(do.call(cbind, classScores))[, pathway := names(classScores[[1]])]
    pvalsDT <-
        as.data.table(do.call(cbind, pvals))[, pathway := names(classScores[[1]])]

    pathScoreDT <-  melt(pathwayScoresDT, id.vars="pathway",
                          value.name="pathScore", variable.name="class",
                          variable.factor=FALSE
                          )
    pDT <- melt(pvalsDT, id.vars="pathway",
                      value.name="pval", variable.name="class",
                      variable.factor=FALSE)

    pathwayStats <- merge(pathScoreDT, pDT, by=c("class", "pathway"))
    pathStats[[pathwayName]] <- pathwayStats[, geneSet := rep(pathwayName, nrow(pathwayStats))]
  }
  pathStatsDT <- rbindlist(pathStats)[is.infinite(pathScore), pathScore := NA]
  return(pathStatsDT)
}


#'
#'
#'
#'
#'
heatmapPathwayScores <- function(pathStatsDT, exclude, significance=0.05) {

    if (!missing(exclude)) {
        pStatsDT <- pathStatsDT[pval < significance & !(class %in% exclude), ][, -'pval']

    } else {
        pStatsDT <- pathStatsDT[pval < significance, ][, -'pval']
    }

    pStatsDT <- dcast(pStatsDT, geneSet + pathway ~ class , value.var='pathScore')
    pStatsDT <- na.omit(pStatsDT)

    splitDT <- split(pStatsDT, by="geneSet")

    preprocMats <- lapply(splitDT, function(DT) {
        mat <- as.matrix(DT[, .SD, .SDcols=-c("geneSet", "pathway")])
        rownames(mat) <- DT$pathway
        mat
    })

    colorFun <- colorRampPalette(c("gold", "white", "black"))
    plotExpressions <- lapply(preprocMats,
                    function(mat)
                        as.expression(call('heatmap.2', x=mat, Rowv=TRUE, Colv=NA,
                                          col=colorFun(15), lhei=c(0.05, 0.20),
                                          scale="column", margins=c(5, 10),
                                          trace="none", dendrogram="row",
                                          cexRow=0.5, cexCol=1))
                    )

    plots <- lapply(plotExpressions, as.ggplot)
    return(plots)
}


########################################################################################################
### Required libraries
#########################################################################################################
require(piano)
require(GSA)
library(gplots)

#########################################################################################
load("../../data/meta_effectsize.RData")
load('../../data/common_genes_cohorts_new.RData')


metacluster1_genes=abs(mm$Metacluster1) > 0.5
metacluster2_genes=abs(mm$Metacluster2) > 0.5
metacluster3_genes=abs(mm$Metacluster3) > 0.5


c3 <- cbind(metacluster1_genes, metacluster2_genes, metacluster3_genes)

metacluster1_genes=mm[order(abs(mm$Metacluster1), decreasing = TRUE)[1:1000],]
metacluster2_genes=mm[order(abs(mm$Metacluster2), decreasing = TRUE)[1:1000],]
metacluster3_genes=mm[order(abs(mm$Metacluster3), decreasing = TRUE)[1:1000],]


##################### Gene enrichment


reference_genes= rownames(cohorts1$pcsi)

hallmark= loadGSC("../../data/Pathway_files/h.all.v6.1.symbols.gmt")
results_hallmark=runGSAhyper(rownames(metacluster1_genes),gsc= hallmark,  adjMethod="fdr",  universe= reference_genes)
results_hallmark$resTab[which(results_hallmark$resTab[,2]<0.05),]
clust1=qnorm(1 - (results_hallmark$pvalues/2))

results_hallmark=runGSAhyper(rownames(metacluster2_genes),gsc= hallmark,  adjMethod="fdr",  universe= reference_genes)
results_hallmark$resTab[which(results_hallmark$resTab[,2]<0.05),]
clust2=qnorm(1 - (results_hallmark$pvalues/2))

results_hallmark=runGSAhyper(rownames(metacluster3_genes),gsc= hallmark,  adjMethod="fdr",  universe= reference_genes)
results_hallmark$resTab[which(results_hallmark$resTab[,1]<0.05),]
clust3=qnorm(1 - (results_hallmark$pvalues/2))

zz=data.frame(clust1=clust1, clust2=clust2, clust3=clust3)
zz[,1][is.infinite(zz[,1])] = NA
rownames(zz)=sub("HALLMARK_","",rownames(zz))
colnames(zz)= c("Basal","Exocrine","Classical")

colfunc <- colorRampPalette(c("gold", "white", "black"))

pdf("../results/Hallmark_genes_pathway.pdf")

heatmap.2(data.matrix(zz),Rowv=TRUE, Colv=NA, col=colfunc(15), scale="column", margins=c(5,10), trace="none",dendrogram="row",cexRow=0.5,cexCol=1)

####### Conanical
genesets=loadGSC("../../data/Pathway_files/c2.cp.v6.1.symbols.gmt", type="auto")

results=runGSAhyper(rownames(metacluster1_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$p.adj < 0.05),]
clust1=qnorm(1 - (results$pvalues/2))
pval1=results$pvalues

results=runGSAhyper(rownames(metacluster2_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$p.adj < 0.05),]
clust2=qnorm(1 - (results$pvalues/2))
pval2=results$pvalues

results=runGSAhyper(rownames(metacluster3_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$pvalues < 0.01),]
clust3=qnorm(1 - (results$pvalues/2))
pval3=results$pvalues
#
# results=runGSAhyper(rownames(metacluster4_genes),gsc= genesets,  adjMethod="fdr", universe= reference_genes)
# results$resTab[which(results$p.adj < 0.05),]
# clust4=qnorm(1 - (results$pvalues/2))
# pval4=results$pvalues


zz=data.frame(clust1=clust1, clust2=clust2, clust3=clust3)
pp=data.frame(clust1=pval1, clust2=pval2, clust3=pval3)

count=0
ss=vector()
for(i in 1:dim(pp)[1]){
  if( length(which(pp[i,1:3]>0.05))==3){
    ss[count]=i
    count=count+1
  }
}

pp=pp[-ss,]
zz=zz[-ss,]

zz[,1][is.infinite(zz[,1])] = NA
colnames(zz)= c("Basal","Exocrine","Classical")

colfunc <- colorRampPalette(c("gold", "white", "black"))

pdf("../results/Conanical_pathway.pdf")
heatmap.2(data.matrix(zz),Rowv=TRUE, Colv=NA, col=colfunc(15), scale="column", margins=c(5,10),lhei = c(0.05, 0.20), trace="none",dendrogram="row",cexRow=0.5,cexCol=1)


