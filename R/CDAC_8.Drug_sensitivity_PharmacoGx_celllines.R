########################################################################################################
### Required libraries
#########################################################################################################
library(ggpubr)
library(survcomp)
library(gridExtra)

########################################################################################################
### Load data
#########################################################################################################

load('../../data/pharmacogx_pdac_celllines1.RData')
load('../../data/sensitivity_celllines_auc_recomputed.RData')
source('../../code/Functions/meta_subtypes_celllines.R')
load('../../data/pharmacogx_celllines_subtypes.RData')
source('../../code/Functions/mgsub_function.R')

#########################################################################################################
#########################################################################################################


length(which(rownames(sensitivity_cellines$GCSI) %in% rownames(sensitivity_cellines$CCLE)))
drugs= intersect(intersect(rownames(sensitivity_cellines$GCSI),rownames(sensitivity_cellines$CCLE)),rownames(sensitivity_cellines$CTRPv2))

drugs1=drugs
########### GCSI drugs, subtypes, 17009 genes

rownames(celline_subtypes$gcsi_subtypes) = colnames(panc_celllines1$GCSI)
A1=mgsub("[-]",".", colnames(sensitivity_cellines$GCSI))
A2=mgsub(" ",".",A1)
colnames(sensitivity_cellines$GCSI)=A2

length(which(rownames(celline_subtypes$gcsi_subtypes) %in% colnames(sensitivity_cellines$GCSI)))
celllines= rownames(celline_subtypes$gcsi_subtypes)[which(rownames(celline_subtypes$gcsi_subtypes)  %in% colnames(sensitivity_cellines$GCSI))]

gcsi_subtypes = celline_subtypes$gcsi_subtypes[celllines,]
gcsi_ic50 = sensitivity_cellines$GCSI[drugs, celllines]
a=mgsub(c("1","3"), c("Basal","Classical"),  gcsi_subtypes[,2])
gcsi_subtypes$meta_class= a


############ CCLE  drugs, subtypes, 17009 genes

rownames(celline_subtypes$ccle_subtypes) = colnames(panc_celllines1$CCLE)
length(which(rownames(celline_subtypes$ccle_subtypes) %in% colnames(sensitivity_cellines$CCLE)))


celllines= rownames(celline_subtypes$ccle_subtypes)[which(rownames(celline_subtypes$ccle_subtypes) %in% colnames(sensitivity_cellines$CCLE))]


ccle_subtypes = celline_subtypes$ccle_subtypes[celllines,]
ccle_ic50 = sensitivity_cellines$CCLE[drugs, celllines]
#ctrpv2_ic50= -1 * log10(ctrpv2_ic50)


############ CTRPV2  drugs, subtypes, 17009 genes

rownames(celline_subtypes$ctrpv2_subtypes) = colnames(panc_celllines1$CTRPV2)
A1=mgsub("[-]",".", colnames(sensitivity_cellines$CTRPv2))
A2=mgsub(" ",".",A1)
colnames(sensitivity_cellines$CTRPv2)=A2


length(which(rownames(celline_subtypes$ctrpv2_subtypes) %in% colnames(sensitivity_cellines$CTRPv2)))
celllines= rownames(celline_subtypes$ctrpv2_subtypes)[which(rownames(celline_subtypes$ctrpv2_subtypes) %in% colnames(sensitivity_cellines$CTRPv2))]


ctrpv2_subtypes = celline_subtypes$ctrpv2_subtypes[celllines,]
ctrpv2_ic50 = sensitivity_cellines$CTRPv2[drugs, celllines]
#ctrpv2_ic50= -1 * log10(ctrpv2_ic50)

p1=sapply(1:length(drugs), function(x)wilcox.test(ccle_ic50[x,]~ ccle_subtypes$subtypes)$p.value)
p2=sapply(1:length(drugs), function(x)wilcox.test(ctrpv2_ic50[x,]~ ctrpv2_subtypes$subtypes)$p.value)
p3=sapply(1:length(drugs), function(x)wilcox.test(gcsi_ic50[x,]~ gcsi_subtypes$subtypes)$p.value)


#meta_p=sapply(1:length(p1), function(x) combine.test(p=c(p1[x], p2[x],p3[x]), weight=c(length(colnames(ccle_ic50)),length( colnames(ctrpv2_ic50)), length( colnames(gcsi_ic50))),method="z.transform"))
#data.frame(meta_p, drugs)




name=function(x){
  ss=ifelse(x==1,"Basal","Classical")
return(ss)  
}

abc1 = sapply(1:length(drugs), function(x) concordance.index(as.numeric(ccle_ic50[x,]), cl= ifelse(ccle_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether"))
abc2 = sapply(1:length(drugs), function(x) concordance.index(as.numeric(ctrpv2_ic50[x,]), cl= ifelse(ctrpv2_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether"))
abc3 = sapply(1:length(drugs), function(x) concordance.index(as.numeric(gcsi_ic50[x,]), cl= ifelse(gcsi_subtypes$subtypes == 1, 1, 0),na.rm = TRUE, method="noether"))

meta_ci=vector()
meta_se=vector()
meta_pvalue=vector()
pdf("../results/pharmaocogx_cellllines_drugs.pdf")

for(i in 1:length(rownames(ccle_ic50))){
mat= data.frame(cellines= colnames(ccle_ic50), drug=ccle_ic50[i,], subtype= name(ccle_subtypes$subtypes) )
mat1= mat[ !is.na(mat$drug),]
p1=ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype", title= paste(rownames(ccle_ic50)[i]," " ,"CI: ", round(abc1[,i][[1]],2), " (P: ",  round(abc1[,i][[5]],2),")", sep="") ,  legend = "none" )+  scale_fill_brewer(palette="Set2")

mat= data.frame(cellines= colnames(ctrpv2_ic50), drug=ctrpv2_ic50[i,], subtype= name(ctrpv2_subtypes$subtypes ))
mat1= mat[ !is.na(mat$drug),]
p2=ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype",title= paste("CI: ", round(abc2[,i][[1]],2), " (P: ",  round(abc2[,i][[5]],2),")", sep="") ,  legend = "none" )+  scale_fill_brewer(palette="Set2")

mat= data.frame(cellines= colnames(gcsi_ic50), drug=gcsi_ic50[i,], subtype= name(gcsi_subtypes$subtypes) )
mat1= mat[ !is.na(mat$drug),]
p3=ggboxplot(mat1 , x = "subtype", y = "drug", fill = "subtype", add = "jitter",ylab = "AAC", xlab = "Subtype", title= paste("CI: ", round(abc3[,i][[1]],2), " (P: ",  round(abc3[,i][[5]],2),")", sep="") ,  legend = "none" )+ scale_fill_brewer(palette="Set2")
pdf(paste("../results/PharmacoGx_drugs",i,".pdf", sep=""))

grid.arrange(p1,p2,p3, nrow=1, ncol=3)
dev.off()
meta_ci[i]=combine.est(x= c( abc1[,i][[1]],  abc2[,i][[1]], abc3[,i][[1]]), x.se=  c( abc1[,i][[2]],  abc2[,i][[2]], abc3[,i][[2]]), hetero = TRUE, na.rm = TRUE)$estimate
meta_se[i]=combine.est(x= c( abc1[,i][[1]],  abc2[,i][[1]], abc3[,i][[1]]), x.se=  c( abc1[,i][[2]],  abc2[,i][[2]], abc3[,i][[2]]), hetero = TRUE, na.rm = TRUE)$se

#meta_pvalue[i]= combine.test(c( abc1[,i][[5]],  abc2[,i][[5]], abc3[,i][[5]]), c(length(colnames(ccle_ic50)),length( colnames(ctrpv2_ic50)), length( colnames(gcsi_ic50))),method="z.transform")
meta_pvalue[i] <- pnorm((meta_ci[i] - 0.5)/meta_se[i], lower.tail = meta_ci[i] < 0.5) * 2

}

data.frame(drugs, meta_ci, meta_pvalue)

#### For running single cohort and 2 datasets, code is located on github