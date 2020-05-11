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

########################################################################################################
### Biomarker plotting function
#########################################################################################################

biomarker_plot<- function(gene_name){
  cc= c(scale(t(cohorts1$pcsi[gene_name,]), center = TRUE) , scale(t(cohorts1$tcga[gene_name,]), center = TRUE), scale(t(cohorts1$icgc_seq[gene_name,]), center = TRUE),
        scale(t(cohorts1$kirby[gene_name,]), center = TRUE) , scale(t(cohorts1$ouh[gene_name,]), center = TRUE), scale(t(cohorts1$winter[gene_name,]), center = TRUE),
        scale(t(cohorts1$collisson[gene_name,]), center = TRUE) ,scale(t(cohorts1$zhang[gene_name,]), center = TRUE), scale(t(cohorts1$chen[gene_name,]), center = TRUE),
        scale(t(cohorts1$unc[gene_name,]), center = TRUE) , scale(t(cohorts1$icgc_arr[gene_name,]), center = TRUE), scale(t(cohorts1$balagurunathan[gene_name,]), center = TRUE),
        scale(t(cohorts1$pei[gene_name,]), center = TRUE) , scale(t(cohorts1$grutzmann[gene_name,]), center = TRUE), scale(t(cohorts1$badea[gene_name,]), center = TRUE),
        
        scale(t(cohorts1$hamidi[gene_name,]), center = TRUE),  scale(t(cohorts1$yang[gene_name,]), center = TRUE),  scale(t(cohorts1$lunardi[gene_name,]), center = TRUE),
        scale(t(cohorts1$janky[gene_name,]), center = TRUE),  scale(t(cohorts1$bauer[gene_name,]), center = TRUE),  scale(t(cohorts1$haider[gene_name,]), center = TRUE)
  )
  
  
  dd= c(clusters$PCSI$meta_classes, clusters$TCGA$meta_classes, clusters$ICGC_seq$meta_classes, 
        clusters$Kirby$meta_classes, clusters$OUH$meta_classes, clusters$Winter$meta_classes,
        clusters$Collisson$meta_classes, clusters$Zhang$meta_classes, clusters$Chen$meta_classes,
        clusters$UNC$meta_classes, clusters$ICGC_arr$meta_classes,clusters$Balagurunathan$meta_classes,
        clusters$Pei$meta_classes, clusters$Grutzmann$meta_classes, clusters$Badea$meta_classes,
        clusters$hamidi$meta_classes, clusters$yang$meta_classes, clusters$lunardi$meta_classes,
        clusters$janky$meta_classes, clusters$bauer$meta_classes, clusters$haider$meta_classes
        
  )
  
  
  aa=c(rep("01.PCSI", length(cohorts1$pcsi[gene_name,])), 
       rep("02.TCGA", length( cohorts1$tcga[gene_name,])),  
       rep("03.ICGC-seq", length(cohorts1$icgc_seq[gene_name,])),  
       rep("04.Kirby", length( cohorts1$kirby[gene_name,] )),
       rep("05.OUH", length(cohorts1$ouh[gene_name,])),
       rep("06.Winter", length(cohorts1$winter[gene_name,])),
       rep("07.Collisson", length(cohorts1$collisson[gene_name,])),
       rep("08.Zhang", length(cohorts1$zhang[gene_name,])),
       rep("09.Chen", length(cohorts1$chen[gene_name,])),
       rep("10.UNC", length( cohorts1$unc[gene_name,] )),
       rep("11.ICGC-array", length( cohorts1$icgc_arr[gene_name,])),
       rep("12.Balagurnathan", length(cohorts1$balagurunathan[gene_name,])),
       rep("13.Pei", length( cohorts1$pei[gene_name,] )),
       rep("14.Grutzmann", length(cohorts1$grutzmann[gene_name,])),
       rep("15.Badea", length(cohorts1$badea[gene_name,])),
       rep("16.Hamidi", length(cohorts1$hamidi[gene_name,])),
       rep("17.Yang", length(cohorts1$yang[gene_name,])),
       rep("18.Lunardi", length(cohorts1$lunardi[gene_name,])),
       rep("19.Janky", length(cohorts1$janky[gene_name,])),
       rep("20.Bauer", length(cohorts1$bauer[gene_name,])),
       rep("21.Haider", length(cohorts1$haider[gene_name,]))
       
  )  
  
  dd[which(dd == "1")] ="Basal"
  dd[which(dd == "2")] ="Exocrine" 
  dd[which(dd == "3")] = "Classical"
  
  
  
  
  df= data.frame(Label=dd, variable=aa  ,value=cc)
  zz = which(is.na(df$Label) )
  
  
  df= df[-zz,]
  #df= df[-which(df$Label %in% "3"),]
  df1=df
  require(ggplot2)
  
    Q<-ggplot(data = df1, aes(x=variable, y=value)) + geom_boxplot(aes(fill=Label), width=0.5) +  
    scale_fill_manual(values=c("gold", "tan","red","pink","seagreen")) + ylab(gene_name)+ 
    guides(fill=guide_legend(title="Subtypes"))+ 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  my_comparisons <- list( c("Basal","Classical"), c("Classical", "Exocrine"), c("Exocrine", "Basal"))
  p <- ggboxplot(df1, x = "Label", y = "value",
                 color = "Label", palette = "jco",
                 add = "jitter", xlab = "", ylab=gene_name, legend =NULL)
  
  p = p + stat_compare_means(aes(label = ..p.signif..),comparisons = my_comparisons)
 #if want to print p values remove aes(label=..p.signif..)
  A=p
   return(A)}


p1=biomarker_plot("GATA6")  # Classical
p2=biomarker_plot("EGFR")   # Basal
p3=biomarker_plot("SNAI2")   # Basal
p4=biomarker_plot("KRT81")   # Basal
p5=biomarker_plot("CYP3A5")    # Classical
p6=biomarker_plot("HNF1A")  # Classical
p7=biomarker_plot("REG1A")  # Exocrine
p8=biomarker_plot("REG1B")  # Exocrine
p9=biomarker_plot("REG3A")  # Exocrine
p10= biomarker_plot("CEL")    # Exocrine

pdf("../results/biomarkers_basal.pdf")
grid.arrange(p2,p3,p4,nrow=2, ncol=2)

pdf("../results/biomarkers_classical.pdf")
grid.arrange(p1,p5,p6, nrow=2, ncol=2)

pdf("../results/biomarkers_exocrine.pdf")
grid.arrange(p7,p8,p10 ,nrow=2, ncol=2)

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################





