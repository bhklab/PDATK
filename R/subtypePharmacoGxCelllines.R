########################################################################################################
### Required libraries
#########################################################################################################

library(randomForest)
library(limma)

########################################################################################################
### Load data
#########################################################################################################

load('../../data/Single_sample_classifier/gg.RData')
load('../../data/Single_sample_classifier/binary_rf_model.RData')
load('../../data/Single_sample_classifier/RF_response.RData')
load('../../data/Single_sample_classifier/output1_RF.RData')
load('../../data/pharmacogx_pdac_celllines1.RData')


source('../../code/Functions/mgsub_function.R')
source('../../code/Functions/meta_subtypes_celllines.R')

##################################################################################################################
#### CELLLINES

dataset_format <- function(dataset){
  
  dataset1=avereps(dataset)
  dataset1=data.frame(dataset1)
  dataset2= sapply(dataset1, function(x) as.numeric(as.character(x)))
 rownames(dataset2)=rownames(dataset1)
  return(dataset2)
}

gcsi = dataset_format(panc_celllines1$GCSI)
gdsc = dataset_format(panc_celllines1$GDSC)
ccle = dataset_format(panc_celllines1$CCLE)
ctrpv2= dataset_format(panc_celllines1$CTRPV2)
#rownames(gcsi)[9189] ="MIA2"  ### Aliase for CTAGE5 present in metagene
#rownames(gcsi)[4363] ="APOC2"  ### Aliase for APC2 present in metagene



gcsi_subtypes = meta_subtype_cellines(gcsi) 
ccle_subtypes = meta_subtype_cellines(ccle) 
gdsc_subtypes = meta_subtype_cellines(gdsc) 
ctrpv2_subtypes = meta_subtype_cellines(ctrpv2) 

#### COMPARE CELLLINES SUBTYPING CALLS



aa=merge(gcsi_subtypes, ccle_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]
aa[,1][which(as.character(aa[,2]) != as.character(aa[,6]))]
aa[which(as.character(aa[,2]) != as.character(aa[,6])),]

aa=merge( ccle_subtypes, ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]
aa[,1][which(as.character(aa[,2]) != as.character(aa[,6]))]
aa[which(as.character(aa[,2]) != as.character(aa[,6])),]


aa=merge(gcsi_subtypes,  ctrpv2_subtypes, by.x="id", by.y="id")
length(which(as.character(aa[,2]) == as.character(aa[,6])))/dim(aa)[1]
aa[,1][which(as.character(aa[,2]) != as.character(aa[,6]))]
aa[which(as.character(aa[,2]) != as.character(aa[,6])),]

rownames( gcsi_subtypes)=gcsi_subtypes[,1]
rownames( ccle_subtypes)=ccle_subtypes[,1]
rownames( gdsc_subtypes)=gdsc_subtypes[,1]
rownames( ctrpv2_subtypes)=ctrpv2_subtypes[,1]
celline_subtypes=list(gcsi_subtypes=  gcsi_subtypes, ccle_subtypes= ccle_subtypes, gdsc_subtypes=gdsc_subtypes, ctrpv2_subtypes= ctrpv2_subtypes)

save(celline_subtypes, file="../results/pharmacogx_celllines_subtypes.RData")
