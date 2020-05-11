########################################################################################################
### Required libraries
#########################################################################################################
library(ggpubr)
library(survcomp)
library(gridExtra)
library(limma)
library(randomForest)
library(limma)

########################################################################################################
### Load data
#########################################################################################################

load('../../data/Single_sample_classifier/gg.RData')
load('../../data/Single_sample_classifier/binary_rf_model.RData')
load('../../data/Single_sample_classifier/RF_response.RData')
load('../../data/Single_sample_classifier/output1_RF.RData')

source('../../code/Functions/meta_subtypes_celllines.R')
source('../../code/Functions/mgsub_function.R')

######################  PanCuRx Organoids and Xenografts #####################################################
xx= read.table('../../data/PanCuRx_rna_O_X.txt', sep="\t", header = TRUE)

mat=xx[, 4:ncol(xx)]
log_mat= as.matrix(log2(mat+1))
rownames(log_mat) = xx[,2]
avg_mat=avereps(log_mat)
classes= meta_subtype_cellines(avg_mat) 

classes

dataset_format <- function(dataset){
  
  dataset1=avereps(dataset)
  dataset1=data.frame(dataset1)
  dataset2= sapply(dataset1, function(x) as.numeric(as.character(x)))
 rownames(dataset2)=rownames(dataset1)
  return(dataset2)
}


####### PDXE clustering
data <- readRDS('../../data/PDXE_PDAC_RNASeq.Rda')
library(Biobase)
dd=exprs(data)


pdxe_subtypes = meta_subtype_cellines(dataset_format(dd)) 
table(pdxe_subtypes$subtypes)


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
######################  PDXE drug comparison #####################################################

################################################################################################################################
#load("/Users/vandanasandhu/Desktop/Subtyping_PDACs/Datasets/PDXe_classes.RData")
load('../../data/PDXe_drug_response.RData')

pdxe_drugs = unique(meta$drug)
pdxe_subtypes$subtypes=mgsub(c("1","3"), c("Basal","Classical"),  pdxe_subtypes$subtypes)

pdf("../results/pdxe_drugs.pdf")

for(i in 1: length(pdxe_drugs)){
  pdxe_mm = meta[which(meta$drug == pdxe_drugs[i]),]
  pdxe_mm$patient.id=mgsub("-",".",pdxe_mm$patient.id)
  
  mm = merge(pdxe_mm, pdxe_subtypes, by.x="patient.id", by.y="id")
  
  p=ggboxplot(mm , x = "subtypes", y = "AAC", fill = "subtypes", add = "jitter",ylab = "AAC", xlab = "Subtype", title= mm$drug[1])+ stat_compare_means(method = "kruskal.test")+theme_bw(base_family = 'Helvetica')
  plot(p)
}


