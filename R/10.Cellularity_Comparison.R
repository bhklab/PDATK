########################################################################################################
### Required libraries
#########################################################################################################
library(ggpubr)
library(ggplot2)
library(gridExtra)


########################################################################################################
### Load data
#########################################################################################################

source('../../code/Functions/mgsub_function.R')
load('../data/clusters.RData')
#########################################################################################################
#########################################################################################################

pcsi_cellularity = read.table('../data/Cellularity/pcsi_cellularity.txt', header=TRUE, sep="\t")
pp = data.frame(id= names(clusters$PCSI$meta_classes), classes= clusters$PCSI$meta_classes)
mm = merge(pcsi_cellularity, pp, by.x="Tumor", by.y = "id")
mm$cohorts= rep("PCSI", length(mm$classes))
colnames(mm) = c("Id","Tumor percentage", "Subtypes", "Cohorts")


icgc_cellularity = read.table('../data/Cellularity/icgc_seq_cellularity.txt', header=TRUE, sep="\t")
pp = data.frame(id= names(clusters$ICGC_seq$meta_classes), classes= clusters$ICGC_seq$meta_classes)
mm1 = merge(icgc_cellularity, pp, by.x="Icgc_id", by.y = "id")
mm1$cohorts= rep("ICGC", length(mm1$classes))
mm1=mm1[,c(-2,-3,-4,-6)]
colnames(mm1) = c("Id","Tumor percentage", "Subtypes", "Cohorts")



tcga_cellularity = read.table('../data/Cellularity/tcga_cellularity.txt', header=TRUE, sep="\t")
pp = data.frame(id= names(clusters$TCGA$meta_classes), classes= clusters$TCGA$meta_classes)
mm2 = merge(tcga_cellularity, pp, by.x="Sample_ID", by.y = "id")
mm2$cohorts= rep("TCGA", length(mm2$classes))
colnames(mm2) = c("Id","Tumor percentage", "Subtypes", "Cohorts")


ouh_cellularity = read.table('../data/Cellularity/ouh_cellularity.txt', header=TRUE, sep="\t")
pp = data.frame(id= names(clusters$OUH$meta_classes), classes= clusters$OUH$meta_classes)
mm3 = merge(ouh_cellularity, pp, by.x="id", by.y = "id")
mm3$cohorts= rep("OUH", length(mm3$classes))
mm3$tumor_percentage= mm3$tumor_percentage/100
colnames(mm3) = c("Id","Tumor percentage", "Subtypes", "Cohorts")

haider_cellularity = read.table('../data/Cellularity/haider_cellularity.txt', header=TRUE, sep="\t")
pp = data.frame(id= names(clusters$haider$meta_classes), classes= clusters$haider$meta_classes)
mm4 = merge(haider_cellularity, pp, by.x="id", by.y = "id")
mm4$cellularity =mm4$cellularity/100
mm4$cohorts= rep("Haider", length(mm4$classes))
colnames(mm4) = c("Id","Tumor percentage", "Subtypes", "Cohorts")


mat= rbind(mm,mm1,mm2,mm3,mm4)
colnames(mat)= c("Id","Cellularity", "Subtypes", "Cohorts")
mat$Subtypes = mgsub(c("1","2","3"), c("Basal","Exocrine","Classical"), mat$Subtypes)


pdf("../results/PDAC_cellularity.pdf")


ggplot(data = mat, aes(Cohorts, Cellularity)) + geom_boxplot(aes(fill=Subtypes), width=0.5)+
  scale_fill_brewer(palette = "Pastel1") + ylab("Cellularity")+ 
  guides(fill=guide_legend(title="Cellularity"))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
