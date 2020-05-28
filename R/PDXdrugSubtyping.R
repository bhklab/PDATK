#'
#'
#'
#'
#'
#'
preprocPDXdata <- function(PDXdata, classifModel, topGenesDT, trainData, trainLabels) {
  # Deal with matrix and df input
  if (is.matrix(PDXdata)) {
    PDXdata <- as.data.table(PDXdata)[, gene_name := rownames(PDXdata)]
  } else if (!is.data.table(PDXdata)) {
    if (!grepl('gene', colnames(PDXdata), ignore.case=TRUE)) 
      PDXdata$gene_name <- rownames(PDXdata)
    PDXdata <- as.data.table(PDXdata)
  }
  # Get remove annotation columns and find gene name column
  whichNumeric <- lapply(PDXdata, class) %in% c("numeric")
  IDcol <- grep('gene', colnames(PDXdata), value=TRUE, ignore.case=TRUE)[1]
  
  # Convert to matrix, log and average replicates  
  dataMat <- as.matrix(PDXdata[, .SD, .SDcols=whichNumeric])
  rownames(dataMat) <- PDXdata[[IDcol]]
  normMat <- avereps(log2(dataMat + 1))
  
  # Predict subtypes
  subtypeDT <- predictSampleMetaClass(normMat, classifModel, topGenesDT,
                                      trainData, trainLabels)
  # normDT <- as.data.table(t(normMat))[, sample := colnames(normMat)]
  # 
  # preprocDT <- merge(subtypeDT, normDT, by="sample")
  # meltPreproc
  return(subtypeDT)
}

#'
#'
#'
#'
#'
boxplotPDXsubtypePerDrug <- function(PDXmergedDT) {
  splitOnDrug <- split(PDXmergedDT, by="drug")
  names(splitOnDrug) <- unlist(lapply(splitOnDrug, function(DT) unique(DT$drug)))
  plots <- lapply(splitOnDrug, 
                  function(DT)
                    ggboxplot(DT, x="predClass", y="AAC", color="predClass", 
                              add="jitter", ylab="AAC", xlab="Subtype", 
                              pallette="jco", title=unique(DT$drug), 
                              legend="none") +
                              stat_compare_means(method="kruskal.test")
                              )
  return(plots)
}

#' Take a long list of plots and chunk it into a list of plot chunck x chunk
#'    plot grids.
#'
#' @param plots A \code{list} of grob or ggplot objects which can be
#'     passed to the `ggpubr::ggarrange` function.
#' @param chunk A \code{numeric} vector indicating the integer number
#'     of plots per plot grid. All returned grids have equal ncol and nrow (
#'     i.e., are ~ square).
#'
#' @importFrom BBmisc chunk
#' @export
ggarrangePlotL <- function(plotL, chunk) {
  nrow <- .ceilSqrt(chunk)
  plotChunks <- chunk(plots, chunk.size=chunk)
  plotGrids <- lapply(plotChunks, 
                      function(plots, ncol, nrow) 
                          ggarrange(plotlist=plots, ncol=ncol, nrow=nrow),
                      ncol=ceiling(chunk/nrow),
                      nrow=nrow)
  return(plotGrids)
}

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


