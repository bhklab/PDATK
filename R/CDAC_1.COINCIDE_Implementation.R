########################################################################################################
### Required libraries
#########################################################################################################

require(IDPmisc)
library(dplyr)
library(matrixStats)
library(factoextra)
library(NbClust)
library(cluster)
library(ConsensusClusterPlus)
library(clusterRepro)
library(igraph)
library(gtools)
library(plyr)
library(survival)
library(KMsurv)
library(survminer)
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

########################################################################################################
### Loading Normal, GTEx and other cohorts
#########################################################################################################

load('../../data/common_genes_cohorts_new.RData')       
load('../../data/normal_GTEx.RData')
load('../../data/adjacent_normal.RData')


########################################################################################################
### Sourcing helping functions
#########################################################################################################

source('../../code/Functions/mad_Cal.R')
source('../../code/Functions/cluster_me_cont.R')
source('../../code/Functions/comare_interclusters_cont_pearson.R')
source('../../code/Functions/mgsub_function.R')
source('../../code/Functions/msm_threshold.r')


########################################################################################################
### Ranking the genes based on MAD
#########################################################################################################

ll=vector()

pcsi_mad=mad_cal(cohorts1$pcsi)
tcga_mad=mad_cal(cohorts1$tcga)
kirby_mad=mad_cal(cohorts1$kirby)

icgc_seq_mad=mad_cal(cohorts1$icgc_seq)
ouh_mad=mad_cal(cohorts1$ouh)
winter_mad=mad_cal(cohorts1$winter)

collisson_mad=mad_cal(cohorts1$collisson)
zhang_mad=mad_cal(cohorts1$zhang)
chen_mad=mad_cal(cohorts1$chen)

unc_mad=mad_cal(cohorts1$unc)
icgc_arr_mad=mad_cal(cohorts1$icgc_arr)
badea_mad=mad_cal(cohorts1$badea)

balagurunathan_mad=mad_cal(cohorts1$balagurunathan)
grutzmann_mad=mad_cal(cohorts1$grutzmann)
pei_mad=mad_cal(cohorts1$pei)

haider_mad= mad_cal(cohorts1$haider)
#van_den_broeck_mad= mad_cal(cohorts1$van_den_broeck)
bauer_mad=mad_cal(cohorts1$bauer)

janky_mad=mad_cal(cohorts1$janky)
lunardi_mad=mad_cal(cohorts1$lunardi)
yang_mad=mad_cal(cohorts1$yang)
hamidi_mad=mad_cal(cohorts1$hamidi)

lm=vector()
seq_mad=vector()
arr_mad =vector()
for(i in 1:10331){
  
  
  lm[i]= weightedMad(c(as.numeric(as.character(pcsi_mad$mad_values[i])),
                       as.numeric(as.character(tcga_mad$mad_values[i])),
                       as.numeric(as.character(kirby_mad$mad_values[i])),
                       as.numeric(as.character(icgc_seq_mad$mad_values[i])),
                       as.numeric(as.character(ouh_mad$mad_values[i])),
                       as.numeric(as.character(winter_mad$mad_values[i])),
                       as.numeric(as.character(collisson_mad$mad_values[i])),
                       as.numeric(as.character(zhang_mad$mad_values[i])),
                       as.numeric(as.character(grutzmann_mad$mad_values[i])),
                       as.numeric(as.character(chen_mad$mad_values[i])),
                       as.numeric(as.character(unc_mad$mad_values[i])),
                       as.numeric(as.character(icgc_arr_mad$mad_values[i])),
                       as.numeric(as.character(balagurunathan_mad$mad_values[i])),
                       as.numeric(as.character(badea_mad$mad_values[i])),
                       as.numeric(as.character(pei_mad$mad_values[i])),
                       as.numeric(as.character(haider_mad$mad_values[i])),
                       as.numeric(as.character(lunardi_mad$mad_values[i])),
                       as.numeric(as.character(yang_mad$mad_values[i])),
                       as.numeric(as.character(hamidi_mad$mad_values[i])),
                       as.numeric(as.character(janky_mad$mad_values[i])),
                       
                       as.numeric(as.character(bauer_mad$mad_values[i]))),
                       
                     c(dim(cohorts1$pcsi)[2], dim(cohorts1$tcga)[2], dim(cohorts1$kirby)[2], dim(cohorts1$icgc_seq)[2],
                       dim(cohorts1$ouh)[2], dim(cohorts1$winter)[2], dim(cohorts1$collisson)[2], dim(cohorts1$zhang)[2],dim(cohorts1$grutzmann)[2],
                       dim(cohorts1$chen)[2], dim(cohorts1$unc)[2], dim(cohorts1$icgc_arr)[2], dim(cohorts1$balagurunathan)[2], dim(cohorts1$badea)[2], 
                       dim(cohorts1$pei)[2], dim(cohorts1$haider)[2],dim(cohorts1$lunardi)[2],dim(cohorts1$yang)[2],dim(cohorts1$hamidi)[2],dim(cohorts1$janky)[2],
                       dim(cohorts1$bauer)[2]),
                     
                       na.rm = TRUE
  )
}



zz=data.frame(genes= rownames(cohorts1$pcsi),weighted_mad = lm)
zz$rank = dense_rank(-zz$weighted_mad)
zz=zz[order(zz$weighted_mad, decreasing = TRUE),]


##########################################################################################################
## Top 20 genes in each cohort

top_20_genes = c(as.character(collisson_mad[order(collisson_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(pcsi_mad[order(pcsi_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(tcga_mad[order(tcga_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(kirby_mad[order(kirby_mad$mad_values, decreasing = TRUE)[1:20],]$genes),  
                 as.character(icgc_seq_mad[order(icgc_seq_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(icgc_arr_mad[order(icgc_arr_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(ouh_mad[order(ouh_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(winter_mad[order(winter_mad$mad_values, decreasing = TRUE)[1:20],]$genes),      
                 as.character(zhang_mad[order(zhang_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(grutzmann_mad[order(grutzmann_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(chen_mad[order(chen_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(unc_mad[order(unc_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(balagurunathan_mad[order(balagurunathan_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(badea_mad[order(badea_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(pei_mad[order(pei_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(haider_mad[order(haider_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(lunardi_mad[order(lunardi_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(yang_mad[order(yang_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(hamidi_mad[order(hamidi_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(janky_mad[order(janky_mad$mad_values, decreasing = TRUE)[1:20],]$genes),
                 as.character(bauer_mad[order(bauer_mad$mad_values, decreasing = TRUE)[1:20],]$genes))



Meta_genes_cohorts = as.character(zz$genes[1:2066])
meta_genes= unique(c(top_20_genes, Meta_genes_cohorts))

save(meta_genes,file='../results/meta_genes.RData')

##################### CLUSTER DATA ##########################################################################

pcsi=t(scale(t(data.matrix(cohorts1$pcsi[as.character(meta_genes),])),center = TRUE, scale = FALSE))
tcga=t(scale(t(data.matrix(cohorts1$tcga[as.character(meta_genes),])),center = TRUE, scale = FALSE))
icgc=t(scale(t(data.matrix(cohorts1$icgc_seq[as.character(meta_genes),])),center = TRUE, scale = FALSE))
kirby=t(scale(t(data.matrix(cohorts1$kirby[as.character(meta_genes),])),center = TRUE, scale = FALSE))
ouh=t(scale(t(data.matrix(cohorts1$ouh[as.character(meta_genes),])),center = TRUE, scale = FALSE))
winter=t(scale(t(data.matrix(cohorts1$winter[as.character(meta_genes),])),center = TRUE, scale = FALSE))
collisson=t(scale(t(data.matrix(cohorts1$collisson[as.character(meta_genes),])),center = TRUE, scale = FALSE))
zhang=t(scale(t(data.matrix(cohorts1$zhang[as.character(meta_genes),])),center = TRUE, scale = FALSE))
chen=t(scale(t(data.matrix(cohorts1$chen[as.character(meta_genes),])),center = TRUE, scale = FALSE))
unc=t(scale(t(data.matrix(cohorts1$unc[as.character(meta_genes),])),center = TRUE, scale = FALSE))
icgc_arr=t(scale(t(data.matrix(cohorts1$icgc_arr[as.character(meta_genes),])),center = TRUE, scale = FALSE))
balagurunathan=t(scale(t(data.matrix(cohorts1$balagurunathan[as.character(meta_genes),])),center = TRUE, scale = FALSE))
pei=t(scale(t(data.matrix(cohorts1$pei[as.character(meta_genes),])),center = TRUE, scale = FALSE))
grutzmann=t(scale(t(data.matrix(cohorts1$grutzmann[as.character(meta_genes),])),center = TRUE, scale = FALSE))
badea=t(scale(t(data.matrix(cohorts1$badea[as.character(meta_genes),])),center = TRUE, scale = FALSE))

haider= t(scale(t(data.matrix(cohorts1$haider[as.character(meta_genes),])),center = TRUE, scale = FALSE))
lunardi= t(scale(t(data.matrix(cohorts1$lunardi[as.character(meta_genes),])),center = TRUE, scale = FALSE))
yang= t(scale(t(data.matrix(cohorts1$yang[as.character(meta_genes),])),center = TRUE, scale = FALSE))
hamidi= t(scale(t(data.matrix(cohorts1$hamidi[as.character(meta_genes),])),center = TRUE, scale = FALSE))
janky= t(scale(t(data.matrix(cohorts1$janky[as.character(meta_genes),])),center = TRUE, scale = FALSE))
bauer= t(scale(t(data.matrix(cohorts1$bauer[as.character(meta_genes),])),center = TRUE, scale = FALSE))


z_PCSI=cluster_me(pcsi, 5, "pearson","hc")
z_TCGA=cluster_me(tcga, 5, "pearson","hc")
z_ICGC=cluster_me(icgc, 5, "pearson","hc")
z_Kirby=cluster_me(kirby, 5, "pearson","hc")
z_OUH=cluster_me(ouh, 5, "pearson","hc")
z_winter=cluster_me(winter, 5, "pearson","hc")
z_Collisson=cluster_me(collisson, 5, "pearson","hc")
z_zhang=cluster_me(zhang, 5, "pearson","hc")
z_chen=cluster_me(chen, 5, "pearson","hc")
z_unc=cluster_me(unc, 5, "pearson","hc")
z_icgc_arr=cluster_me(icgc_arr, 5, "pearson","hc")
z_balagurunathan =cluster_me(balagurunathan, 5, "pearson","hc")
z_pei=cluster_me(pei, 5, "pearson","hc")
z_grutzmann =cluster_me(grutzmann, 5, "pearson","hc")
z_badea= cluster_me(badea, 5, "pearson","hc")
z_haider= cluster_me(haider, 5, "pearson","hc")
z_lunardi= cluster_me(lunardi, 5, "pearson","hc")
z_yang= cluster_me(yang, 5, "pearson","hc")
z_hamidi= cluster_me(hamidi, 5, "pearson","hc")
z_janky= cluster_me(janky, 5, "pearson","hc")
z_bauer= cluster_me(bauer, 5, "pearson","hc")

rownames(normal_gtex)[650] =  "ALPPL2"
which(rownames(data.matrix(normal_gtex)) %in% "MICAL2")
#12716
setdiff(as.character(meta_genes), rownames(data.matrix(normal_gtex)))
#MICALCL"
rownames(normal_gtex)[12716] =  "MICALCL"

z_normal= cluster_me(data.matrix(nn2)[as.character(meta_genes),], 5, "pearson","hc")
z_gtex= cluster_me(data.matrix(normal_gtex)[as.character(meta_genes),], 5, "pearson","hc")


datasets = c("ICGC_seq","PCSI","TCGA","Kirby","OUH","Winter","Collisson","Zhang","Chen","UNC","ICGC_arr","Balagurunathan","Pei","Grutzmann","Badea", "haider","lunardi","yang","hamidi","janky","bauer", "Normals","GTEX_Normals")
clusters = list(z_ICGC, z_PCSI, z_TCGA, z_Kirby, z_OUH, z_winter,z_Collisson,z_zhang,z_chen,z_unc,z_icgc_arr,z_balagurunathan, z_pei, z_grutzmann, z_badea, z_haider, z_lunardi, z_yang, z_hamidi, z_janky, z_bauer, z_normal,z_gtex)
zz=expand.grid(x =1:length(datasets), y = 1:length(datasets))
allPairs <-  zz[-which(zz[,1] == zz[,2]),]
data=list(icgc, pcsi, tcga, kirby, ouh, winter, collisson, zhang, chen, unc,icgc_arr, balagurunathan, pei, grutzmann, badea, haider, lunardi,yang,hamidi,janky,bauer, data.matrix(nn2)[as.character(meta_genes),],data.matrix(normal_gtex)[as.character(meta_genes),])


cores=detectCores()
cl <- makeCluster(cores[1]-1, outfile="") #not to overload your computer
registerDoParallel(cl)


threshold_msm = list()
threshold_msm = foreach(i = 1: dim(allPairs)[1], .packages = c( "doParallel","ConsensusClusterPlus", "clusterRepro","foreach")) %dopar% {                                                        
  msm_thresholds(datasets[allPairs[i,1]], datasets[allPairs[i,2]], clusters[[allPairs[i,1]]],clusters[[allPairs[i,2]]],data.frame(data[allPairs[i,2]]),i)
  
}
pdf("../results/densityplot.pdf")

densityplot(unlist(threshold_msm ))

cores=detectCores()
cl <- makeCluster(cores[1]-1, outfile="") #not to overload your computer
registerDoParallel(cl)


edges=list()
edges= foreach(i = 1: dim(allPairs)[1], .packages = c( "doParallel","ConsensusClusterPlus", "clusterRepro","foreach")) %dopar% {                                                        
  compare_interclusters(datasets[allPairs[i,1]], datasets[allPairs[i,2]], clusters[[allPairs[i,1]]],clusters[[allPairs[i,2]]],data.frame(data[allPairs[i,2]]),i)
  
}
edges_define = matrix(unlist(edges), ncol = 5, byrow = TRUE)
save(edges_define, file='../results/edges.RData')

######### FIND EDGES ###############################
edges_defined = edges_define[,1:2]

######### FIND METACLUSTERS ########################
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
datasets = c("ICGC_seq","PCSI","TCGA","Kirby","OUH","Winter","Collisson","Zhang","Chen","UNC","ICGC_arr","Balagurunathan","Pei","Grutzmann","Badea", "haider","lunardi","yang","hamidi","janky","bauer", "Normals","GTEX_Normals")
clusters = list(z_ICGC, z_PCSI, z_TCGA, z_Kirby, z_OUH, z_winter,z_Collisson,z_zhang,z_chen,z_unc,z_icgc_arr,z_balagurunathan, z_pei, z_grutzmann, z_badea, z_haider, z_lunardi, z_yang, z_hamidi, z_janky, z_bauer, z_normal,z_gtex)

total_clusters=list()
for(i in 1:length(datasets)){
  
  total_clusters[[i]]=paste(datasets[i], 1:clusters[[i]]$optimumK, sep="-" )
  
}

not_found=setdiff(unlist(total_clusters), c5$names)

not_found

########### Samples in classes ###################
datasets = c("ICGC_seq","PCSI","TCGA","Kirby","OUH","Winter","Collisson","Zhang","Chen","UNC","ICGC_arr","Balagurunathan","Pei","Grutzmann","Badea", "haider","lunardi","yang","hamidi","janky","bauer", "Normals","GTEX_Normals")
clusters = list(z_ICGC, z_PCSI, z_TCGA, z_Kirby, z_OUH, z_winter,z_Collisson,z_zhang,z_chen,z_unc,z_icgc_arr,z_balagurunathan, z_pei, z_grutzmann, z_badea, z_haider, z_lunardi, z_yang, z_hamidi, z_janky, z_bauer, z_normal,z_gtex)

samples=list()
for(i in 1: length(unique(membership(c5)))){
  
  clust = names(membership(c5))[which(membership(c5)==i)]  
  
  dataset_no = sapply(1:length(clust), function(x) which(strsplit(clust[x],"-")[[1]][1] == datasets))
  cluster_no = as.numeric(sapply(1:length(clust), function(x) strsplit( names(membership(c5))[which(membership(c5)==i)],"-")[[x]][2]))
  
  samples[[i]]=sapply(1:length(dataset_no),function(x) names(clusters[[dataset_no[x]]]$classes)[which(clusters[[dataset_no[x]]]$classes == cluster_no[x])])
  
}

sample_groups<- list()
for(i in 1:length(samples)){
  
  sample_groups[[i]]=unlist(samples[[i]])
  
}

#save(sample_groups, file="/Users/vandanasandhu/Desktop/Subtyping_PDACs/PDAC_meta_Classes.RData")

########### Samples classes cohort wise #################################################
clusters = list(PCSI= z_PCSI, TCGA=z_TCGA, ICGC_seq=z_ICGC, 
                Kirby = z_Kirby, OUH= z_OUH, Winter= z_winter,
                Collisson=z_Collisson,Zhang=z_zhang,Chen=z_chen,
                UNC=z_unc,ICGC_arr=z_icgc_arr,Balagurunathan=z_balagurunathan, 
                Grutzmann=z_grutzmann,Badea= z_badea,Pei=z_pei,
                hamidi=z_hamidi,yang=z_yang, lunardi=z_lunardi,
                janky=z_janky, bauer=z_bauer, haider=z_haider)
#, Normals=z_normal, GTEX_Normals= z_gtex
#)


x=unlist(strsplit(names(membership(c5)),"-"))
mm=matrix(unlist(x), ncol=2, byrow = TRUE)
mm=data.frame(mm)
colnames(mm)= c("dataset","cohort_clusters")

mm$meta_classes=as.numeric( membership(c5))
mm=mm[order(mm$dataset),]

### ADDING NOT FOUND ROWS
dim(mm)

j=1
for(i in (dim(mm)[1] +1 ):(dim(mm)[1]+length(not_found))){
  
  pp = unlist(strsplit(not_found[j],"-"))
  mm[i,]=c(pp[1],pp[2],NA)
  j=j+1
}

for(i in 1:length(clusters)){
  
  z= which(mm[,1] %in% names(clusters)[i])
  clusters[[i]]$meta_classes=mapvalues(clusters[[i]]$classes, mm[,2][z], mm[,3][z]) 
  
}

save(clusters,file= '../results/clusters.RData')

