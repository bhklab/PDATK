########################################################################################################
### Required libraries
#########################################################################################################

library(randomForest)
library(limma)
library(ggpubr)
library(gridExtra)
library(survival)
#library(survminer)

########################################################################################################
### Load data
#########################################################################################################

load('../../data/Single_sample_classifier/gg.RData')
load('../../data/Single_sample_classifier/binary_rf_model.RData')
load('../../data/Single_sample_classifier/RF_response.RData')
load('../../data/Single_sample_classifier/output1_RF.RData')

source('../../code/Functions/mgsub_function.R')
source('../../code/Functions/meta_subtypes_celllines.R')

####### COMPASS SAMPLES 
load('../../data/log_transformed_compass.RData')
com_sur= read.csv('../../data/compass_survival', sep="\t", header = TRUE)

aa=meta_subtype_cellines(compass_mat) 
table(aa$subtypes)


mm1 = merge(aa, com_sur, by.x="id", by.y ="ID")
mm1$tumor_response=as.numeric(as.character(mm1$tumor_response))
mm1=mm1[order(mm1$tumor_response, decreasing = TRUE),]
mm1= mm1[-which(is.na(mm1$tumor_response) == TRUE),]


ggbarplot(mm1 ,x="Study_ID", y = "tumor_response", fill= "subtypes", color = "subtypes",palette = c("#00AFBB", "#E7B800", "#FC4E07","pink"))


comp_id=gsub("COMP-0","",mm1$Study_ID)
mm1$color =mgsub(c("1","2","3"), c( "blue","cyan", "gold"),  mm1$subtypes )

mm1=mm1[-which(mm1$subtypes ==2), ]

#pdf("../results/compass_waterfall_plot_All_sample.pdf")
#ggboxplot(mm1 , x = "subtypes", y = "tumor_response") + stat_compare_means(method = "kruskal.test")
#barplot(mm1$tumor_response, col=mm1$color, ylab= "Percentage Change in target Lesion", border = "white")
#legend("topright", c( "Basal","Classical", "P=0.03"), fill   = c(  "blue", "gold"), bty = "n")

mm2=mm1        
mm2$subtypes=mgsub(c("1","2","3"), c("Basal","Exocrine","Classical"),  mm1$subtypes )


z1= which(mm1$drug %in% c("FFx", "FFx x 1"))
mm2=mm1[z1,]
mm2$color =mgsub(c("1","2","3","4"), c("blue","cyan", "gold","pink"),  mm2$subtypes )

#mm2=mm2[-which(mm2$subtypes ==c(1,4)),]
pdf("../results/compass_waterfall_plot_FFX.pdf")

#ggboxplot(mm2 , x = "subtypes", y = "tumor_response") + stat_compare_means(method = "kruskal.test")

barplot(as.numeric(mm2$tumor_response), col=mm2$color, ylab= "Percentage Change in target Lesion", border = "white")
legend("topright", c( "Basal","Exocrine","Classical","P=7.0e-03"), fill   = c("blue","cyan", "gold"), bty = "n")



#z2= which(mm1$drug %in% c("GA","Gem only", "GA (PA7)","GA-gem","GA  "))
#mm2=mm1[z2,]
#mm2$color =mgsub(c("1","2","3","4"), c("blue","cyan", "gold","pink"),  mm2$subtypes )

#mm2=mm2[-which(mm2$subtypes ==c(1,4)),]
#pdf("../results/compass_waterfall_plot_gem.pdf")

#ggboxplot(mm2 , x = "subtypes", y = "tumor_response") + stat_compare_means(method = "kruskal.test")

#barplot(as.numeric(mm2$tumor_response), col=mm2$color, ylab= "Percentage Change in target Lesion", border = "white")
#legend("topright", c( "Basal","Classical","P=0.4"), fill   = c("blue", "gold","white"), bty = "n")

#########

pdf("../results/compass_survival.pdf")

mm2$subtypes=mgsub(c(1,3),c("Basal","Classical"), mm2$subtypes)
fit <- survfit(Surv(as.numeric(as.character(mm2$OS)), as.numeric(as.character(mm2$OS_Status)) == 1) ~ mm2$subtypes, data=mm2)
#ggsurvplot(fit, data = mm2, risk.table = TRUE, legend="none", risk.table.ticks.col = TRUE,pval = TRUE,ggtheme = theme_minimal(), #pval.size=5  )

plot( fit,col=c("tomato","green","purple"), lwd=2, p.val=TRUE, xlab="Days", ylab="Survival Probability")
  legend("bottomleft",paste("P =", round(0.03,2), sep = ""), bty = "n")

su=coxph(Surv(as.numeric(as.character(mm2$OS)), as.numeric(as.character(mm2$OS_Status)) == 1) ~ mm2$subtypes, data=mm2)
summary(su)



#############################################################################################################################
### COMPASS#######
######## NLR
aa=meta_subtype_cellines(compass_mat) 

mm1 = merge(aa, com_sur, by.x="id", by.y ="ID")

mm1$NLR=as.numeric(as.character(mm1$NLR))
mm1=mm1[order(mm1$NLR, decreasing = TRUE),]

mm1$color =mgsub(c("1","2","3"), c("gold", "mistyrose", "cyan"),  mm1$subtypes )

mm2=mm1        
mm2$subtypes=mgsub(c("1","2","3"), c("Basal","Exocrine","Classical"),  mm1$subtypes )

zz1= which(mm2$subtypes %in% c("Exocrine","2"))
mm2=mm2[-zz1, ]
mm2=mm2[-which(is.na(mm2$NLR)),]
#pdf("../results/compass_NLR.pdf")

#ggboxplot(mm2 , x = "subtypes", y = "NLR", color = "subtypes", add="jitter")+ stat_compare_means(method = "kruskal.test")



##### Biomarkers


compass_plot <- function(gene_name){
gata6_compass = compass_mat[gene_name,]
dd=as.character(aa$subtypes)

dd[which(dd == "1")] = "Basal"
dd[which(dd == "2")] = "Exocrine"
dd[which(dd == "3")] = "Classical"

df= data.frame(subtype=dd, gene_name=as.numeric(gata6_compass))
which(df$subtype %in% c("Exocrine","2"))
    df= df[-which(df$subtype %in% c("Exocrine","2")), ]

p=ggboxplot(df , x = "subtype", y = "gene_name", color = "subtype", add = "jitter", ylab = gene_name)
    #+ stat_compare_means(method = "kruskal.test")
return(p)
}



pdf("../results/compass_BIOMARKER.pdf")

library(gridExtra)
grid.arrange(compass_plot("GATA6"), compass_plot("CYP3A5"), compass_plot("HNF1A"), compass_plot("EGFR"),compass_plot("KRT81"), compass_plot("SNAI2"), nrow=2, ncol=3)
