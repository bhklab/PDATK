########################################################################################################
### Required libraries
#########################################################################################################
library(ggpubr)
library(gridExtra)
########################################################################################################
### Load data
#########################################################################################################

load('../../data/common_genes_cohorts_new.RData')       
load('../data/clusters.RData')
source('../../code/Functions/mgsub_function.R')

#########################################################################################################

meta_class=list()
cohort_names=list()
for(i in 1: 21){
  
  meta_class[[i]] = clusters[[i]]$meta_classes
  cohort_names[[i]]=names(clusters[[i]]$meta_classes)
}
dd=data.frame(samples=unlist(cohort_names), meta_class=unlist(meta_class))

########################################################################################################
### Biomarker plotting function
#########################################################################################################

signature_score <- function(gene_list, List_name){
  my_comparisons <- list( c("Basal","Classical"), c("Classical", "Exocrine"), c("Exocrine", "Basal"))
  
  gene_list_score=list()
  sample_names=list()
  
  for(i in 1:(length(cohorts1))){
    gene_list_score[[i]] = colMeans( t(scale(t(data.matrix(cohorts1[[i]][gene_list,])))), na.rm = TRUE)
    sample_names[[i]] = names(colMeans(cohorts1[[i]][gene_list,], na.rm = TRUE))
  }
  
  zz= data.frame(gene_list_score=unlist(gene_list_score), class=dd$meta_class)
  REM = which(zz$class %in% c(NA))
  zz$class = mgsub(c("1","2","3"),  c("Basal","Exocrine","Classical"), zz$class)
  
  zz=zz[-REM,]
  #zz= zz[-which(zz$class =="4"),]
  
  p <- ggboxplot(zz, x = "class", y = "gene_list_score",
                 color = "class", palette = "jco",ylab=List_name, 
                 add = "jitter", xlab = "", legend="none")
  
 p= p + stat_compare_means(comparisons = my_comparisons)
  
  
  
  return(p)
}

#1.##################### Cell cycle proliferation score 
pp= read.table('../../data/gene_signature_genelists/CCP_datasets.txt', sep="\t")

p1=signature_score(as.character(pp[,1]), "Proliferation signature") 

#2. ##################### Double Strand break repair
pp= read.table('../../data/gene_signature_genelists/dsbr_signature.txt', sep="\t")
p2=signature_score(as.character(pp[,1]), " Double Strand break repair") 

########################### Hypoxia
pp= read.table('../../data/gene_signature_genelists/HYPOXIA_buffa.txt', sep="\t")
p3=signature_score(as.character(pp[,1]), "Hypoxia signature") 

#3. ##################### iMMUNE CYTOSOLIC ACTIVVITY

pp=c("GZMA","PRF1")
p4= signature_score(as.character(pp), "Immune signature") 

pdf("../results/gene_signatures.pdf")
grid.arrange(p1,p2,p3,p4,nrow=2, ncol=2)

