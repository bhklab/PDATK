########################################################################################################
### Required libraries
#########################################################################################################


#library(survminer)

require(IDPmisc)
library(dplyr)
library(matrixStats)
#library(factoextra)
library(NbClust)
library(cluster)
library(ConsensusClusterPlus)
library(clusterRepro)
library(igraph)
library(gtools)
library(plyr)
library(survival)
library(KMsurv)
library(limma)
library(piano)
require(ggplot2)
library(ggpubr)
library(doParallel)
library(foreach)
library(effsize)
library(vcdExtra)
library(survcomp)
library(RColorBrewer)
library(lattice)

#' Draw a plot showing the network graph of the significant cluster edges
#' 
#'
#' @param clusterEdges A \code{data.table} containing the statistics for significant
#'    inter-cohort cluster comparisons, as returned by `compareClusters`.
#' @param seed
#' @param savePath
#' @param fileName
#'
#' @import igraph
#' @export
plotClusterNetwork <- function(clusterEdges, seed=NULL, savePath, fileName, ...) {
  # Set seed for reproducible results
  if (!is.null(seed)) set.seed(seed)
  
  # Prepare the edge labels
  cohort1 <- vapply(strsplit(clusterEdges$comparison, '-'), `[`,i= 1, character(1))
  cohort2 <- vapply(strsplit(clusterEdges$comparison, '-'), `[`, i=2, character(1))
  edges <- cbind(
    'cluster1'=unlist(mapply(paste0, cohort1, '-', clusterEdges$c1Clust, SIMPLIFY=FALSE)),
    'cluster2'=unlist(mapply(paste0, cohort2, '-', clusterEdges$c2Clust, SIMPLIFY=FALSE))
  )
  
  # Format the network graph
  graph <- graph_from_edgelist(edges)
  coords <- layout_with_fr(graph)
  ugraph <- as.undirected(graph)
  metaClusters <- fastgreedy.community(ugraph, weights=clusterEdges$threshold)
  colours <-  c("blue", brewer.pal(n=8, name="Dark2")[c(4, 5)])
  
  plotNetwork <- call('.plotNetwork', ugraph, metaClusters, coords, colours)
  
  # Plot the network graph
  plot <- base2grob(as.expression(plotNetwork))
  
  if (!missing(savePath) && !missing(fileName))
    ggsave(file.path(savePath, fileName), plot)
  
  grid.draw(plot)
}


#'
#'
#'
#'
#'
.plotNetwork <- function(ugraph, metaClusters, coords, colours) {
  plot(ugraph, 
       vertex.color=colours[membership(metaClusters)],  
       vertex.shape="sphere", 
       vertex.size=6,
       edge.arrow.size=0.5,
       vertex.label.cex=0.8, 
       vertex.label.dist=2, 
       edge.curved=0.1, 
       vertex.color=colours[membership(metaClusters)],
       edge.arrow.size=0.4, 
       layout=layout_with_dh(ugraph), 
       layout=coords)
  legend('topleft',
         legend=c("Basal","Exocrine","Classical"), 
         pt.cex=1.8, pch=21, pt.bg=colours, bty='n', col=colours)
}



















#########################################################################################################
#########################################################################################################

load('../data/edges.RData')
load('../data/clusters.RData')
source('../../code/Functions/mgsub_function.R')


######### FIND EDGES ###############################
edges_defined = edges_define[,1:2]

g = graph_from_edgelist(edges_defined)
coords = layout_with_fr(g)


#### Fast greedy

g1=as.undirected(g)
c5=fastgreedy.community(g1,weights = as.numeric(edges_define[,3]))
table(membership(c5))
plot(c5, g1, layout=coords,vertex.size=6)

l <- layout_with_fr(g1)
#plot(g1, vertex.color=membership(c5),  vertex.size=6,edge.arrow.size=0.5,vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.4, vertex.color="gray",edge.arrow.size=.4,layout=l)


colrs <-  c( "blue", brewer.pal(n = 8, name = "Dark2")[c(4,5)])
shapes <- c("circle", "square", "sphere")
plot(g1, vertex.color=colrs[membership(c5)],  vertex.shape= "sphere", vertex.size=6,edge.arrow.size=0.5,vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.1, vertex.color=colrs[membership(c5)],
     edge.arrow.size=.4, layout=layout.davidson.harel(g1))


zz=sapply(1:length(membership(c5)), function(x) strsplit(names(membership(c5))[x], "-")[[1]][1])
x1= match(zz, unique(zz))

cc= rainbow(22)

pdf("../results/Meta_clusters_Network_community.pdf")

plot.igraph(g1, vertex.color=colrs[membership(c5)],  vertex.shape= "sphere", vertex.size=6,edge.arrow.size=0.5,vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.1, vertex.color=colrs[membership(c5)],
     edge.arrow.size=.4, layout=layout_with_dh(g1), layout=coords)
legend('topleft',legend=c("Basal","Exocrine","Classical"), pt.cex=1.8,pch=21, pt.bg=colrs, bty='n', col=colrs)

################## Samples not clustered
datasets = c("ICGC_seq","PCSI","TCGA","Kirby","OUH","Winter","Collisson","Zhang","Chen","UNC","ICGC_arr","Balagurunathan","Pei","Grutzmann","Badea", "haider","lunardi","yang","hamidi","janky","bauer")

total_clusters=list()
for(i in 1:length(datasets)){
  
  total_clusters[[i]]=paste(datasets[i], 1:clusters[[i]]$optimumK, sep="-" )
  
}

not_found=setdiff(unlist(total_clusters), c5$names)

not_found

########### Samples in classes ###################
datasets = c("ICGC_seq","PCSI","TCGA","Kirby","OUH","Winter","Collisson","Zhang","Chen","UNC","ICGC_arr","Balagurunathan","Pei","Grutzmann","Badea", "haider","lunardi","yang","hamidi","janky","bauer")


################ ICGC-array and ICGC-seq

common_seq = which(names(clusters$ICGC_seq$meta_classes) %in% names(clusters$ICGC_arr$meta_classes))
common_arr = which(names(clusters$ICGC_arr$meta_classes) %in% names(clusters$ICGC_seq$meta_classes))

names(clusters$ICGC_seq$meta_classes[sort(names(clusters$ICGC_seq$meta_classes)[common_seq])])==  names(clusters$ICGC_arr$meta_classes[sort(names(clusters$ICGC_arr$meta_classes)[common_arr])])

seq_classes= clusters$ICGC_seq$meta_classes[sort(names(clusters$ICGC_seq$meta_classes)[common_seq])]
arr_classes= clusters$ICGC_arr$meta_classes[sort(names(clusters$ICGC_arr$meta_classes)[common_arr])]

table(seq_classes, arr_classes)
length(which(seq_classes == arr_classes))
write.table(table(seq_classes, arr_classes), '../results/commmon_ICGC_seq_array.txt')
##################################################################################################
### SURVIVAL 
##########################################################################
##########################################################################
load('../../data/PDAC_Expression_dataset.RData')
yang_survival=read.csv('../../data/yang_survival.txt', sep="\t", header=T)

samples=c(rownames(rs_coh$PCSI_new), rownames(rs_coh$TCGA), rownames(rs_coh$Kirby),
          rownames(rs_coh$ICGC_arr),rownames(rs_coh$ICGC_seq), rownames(rs_coh$UNC),
          rownames(rs_coh$Chen), rownames(rs_coh$Collisson), rownames(rs_coh$Zhang),
          rownames(rs_coh$OUH), rownames(rs_coh$Winter), as.character(yang_survival$ID))


os= c(as.numeric(as.character(rs_coh$PCSI_new$OS)), as.numeric(as.character(rs_coh$TCGA$OS)), 
      as.numeric(as.character(rs_coh$Kirby$OS)),as.numeric(as.character(rs_coh$ICGC_arr$OS)),
      as.numeric(as.character(rs_coh$ICGC_seq$OS)), as.numeric(as.character(rs_coh$UNC$OS)),
      as.numeric(as.character(rs_coh$Chen$OS)), as.numeric(as.character(rs_coh$Collisson$OS)), 
      as.numeric(as.character(rs_coh$Zhang$OS)), as.numeric(as.character(rs_coh$OUH$OS)), 
      as.numeric(as.character(rs_coh$Winter$OS)),as.numeric(as.character(yang_survival$OS)))


os_status= c(as.numeric(as.character(rs_coh$PCSI_new$OS_Status)), as.numeric(as.character(rs_coh$TCGA$OS_Status)), 
             as.numeric(as.character(rs_coh$Kirby$OS_Status)),as.numeric(as.character(rs_coh$ICGC_arr$OS_Status)),
             as.numeric(as.character(rs_coh$ICGC_seq$OS_Status)), as.numeric(as.character(rs_coh$UNC$OS_Status)),
             as.numeric(as.character(rs_coh$Chen$OS_Status)), as.numeric(as.character(rs_coh$Collisson$OS_Status)), 
             as.numeric(as.character(rs_coh$Zhang$OS_Status)), as.numeric(as.character(rs_coh$OUH$OS_Status)), 
             as.numeric(as.character(rs_coh$Winter$OS_Status)), as.numeric(as.character(yang_survival$OS_Status)))



survival= data.frame(samples=samples, os=os, os_status=os_status)


##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################

clusters1=clusters
classes=list()
sample_names=list()
cohorts=list()
for(i in 1:length(clusters1)){
  
  classes[[i]]=clusters1[[i]]$meta_classes
  sample_names[[i]]=names(clusters1[[i]]$meta_classes)
  cohorts[[i]]=rep(names(clusters1)[i],length(names(clusters1[[i]]$meta_classes)))
}

meta_Cluster= data.frame(sample = unlist(sample_names), meta_class= unlist(classes), cohorts=unlist(cohorts))

#save(meta_Cluster, file="/Users/vandanasandhu/Desktop/meta_clusters.RData")

## Removing common samples between ICGC seq and ICGC array
clusters1$ICGC_arr$meta_classes= clusters1$ICGC_arr$meta_classes[- which(names(clusters$ICGC_arr$meta_classes) %in% names(clusters$ICGC_seq$meta_classes))]
clusters1$ICGC_arr$classes= clusters1$ICGC_arr$classes[- which(names(clusters$ICGC_arr$meta_classes) %in% names(clusters$ICGC_seq$meta_classes))]

classes=list()
sample_names=list()
cohorts=list()
for(i in 1:length(clusters1)){
  
  classes[[i]]=clusters1[[i]]$meta_classes
  sample_names[[i]]=names(clusters1[[i]]$meta_classes)
  cohorts[[i]]=rep(names(clusters1)[i],length(names(clusters1[[i]]$meta_classes)))
}

meta_Cluster= data.frame(sample = unlist(sample_names), meta_class= unlist(classes), cohorts=unlist(cohorts))

#save(meta_Cluster, file="/Users/vandanasandhu/Desktop/meta_clusters.RData")

meta_survival_cluster=merge(survival, meta_Cluster, by.x="samples", by.y="sample")
meta_survival_cluster$meta_class= mgsub(c("1","2","3"), c("Basal","Exocrine","Classical"), meta_survival_cluster$meta_class)


#meta_survival_cluster1=meta_survival_cluster[-is.na(meta_survival_cluster$meta_class),]
su=coxph(Surv(os, os_status == 1) ~ meta_class + strata(cohorts), data=meta_survival_cluster)
summary(su)

meta_survival_cluster1=meta_survival_cluster[meta_survival_cluster$meta_class %in% c("Basal","Exocrine","Classical"),]
fit <- survfit(Surv(meta_survival_cluster1$os, meta_survival_cluster1$os_status == 1) ~ meta_survival_cluster1$meta_class + strata(cohorts), data=meta_survival_cluster1)
survdiff(Surv(meta_survival_cluster1$os, meta_survival_cluster1$os_status == 1) ~ meta_survival_cluster1$meta_class + strata(cohorts), data=meta_survival_cluster1)

fit <- survfit(Surv(meta_survival_cluster1$os, meta_survival_cluster1$os_status == 1) ~ meta_survival_cluster1$meta_class , data=meta_survival_cluster1)


#meta_survival_cluster1= meta_survival_cluster1[-which(meta_survival_cluster1$meta_class == "Exocrine"),]
survdiff(Surv(meta_survival_cluster1$os, meta_survival_cluster1$os_status == 1) ~ meta_survival_cluster1$meta_class , data=meta_survival_cluster1)

pdf("../results/survival_metaclusters.pdf")
plot( fit,col=c("tomato","green","purple"), lwd=2, p.val=TRUE, xlab="Days", ylab="Survival Probability")
  legend("bottomleft",paste("P =", round(0.02,2), sep = ""), bty = "n")

#ggsurvplot(fit, data = meta_survival_cluster1, risk.table = TRUE, legend="none", risk.table.ticks.col = TRUE,pval = TRUE,ggtheme #= theme_minimal(), pval.size=5 , 
 #          pval.coord= c(1500, 0.9))




#pairwise_survdiff(Surv(os, os_status == 1) ~ meta_class, data=meta_survival_cluster1)


######################


pdf("../results/cohorts_survival_metaclusters.pdf")

unique_cohorts=unique(meta_survival_cluster1$cohorts)
par(mfrow=c(3,4))

for(i in 1:length(unique(meta_survival_cluster1$cohorts))){
  zz=meta_survival_cluster1[which(meta_survival_cluster1$cohorts == unique_cohorts[i] ),]
  
  fit <- survfit(Surv(zz$os, zz$os_status == 1) ~ zz$meta_class, data=zz)
  
  sd=survdiff(Surv(zz$os, zz$os_status == 1) ~  zz$meta_class, data=zz)
  p.val <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  plot(fit,col=c("tomato","green","purple"), lwd=2, main=unique_cohorts[i], p.val=TRUE)
  legend("bottomleft",paste("P =", round(p.val,2), sep = ""), bty = "n") 
  
}
##########################################################################################
##########################################################################################
