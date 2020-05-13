########################################################################################################
### Load data
#########################################################################################################

load('../data/Oslo_TCGA_score.RData')
SNPpos <- read.table('../data/SNPpos.txt',header=T,sep="\t",row.names=1, stringsAsFactors=F)
ss=read.table('../data/rs-snp-annot.txt',sep="\t",header=T)
load('../data/clusters.RData')
source('../../code/Functions/mgsub_function.R')

#########################################################################################################
#########################################################################################################
#########################################################################################################



score_file = output_score_sc
rownames(score_file) = score_file$ProbeID
score_file = score_file[,2:ncol(score_file)]

resp_in= data.frame(SampleID= names(clusters$OUH$meta_classes), Resp= as.numeric(clusters$OUH$meta_classes))
names(clusters$OUH$meta_classes)= mgsub( "X","",names(clusters$OUH$meta_classes))

resp_in=  resp_in[which(resp_in$SampleID %in%  colnames(score_file)),]

p1.sampleID = resp_in[which(resp_in$Resp==1), 1]
p2.sampleID = resp_in[which(resp_in$Resp==2), 1]
p3.sampleID = resp_in[which(resp_in$Resp==3), 1]
#p4.sampleID = resp_in[which(resp_in$Resp==4), 1]
#p5.sampleID = resp_in[which(resp_in$Resp==5), 1]

countNonzero <- function(data) {
  countUp = apply(data, 1, function(r) { sum(r == 1) } )
  countDown = apply(data, 1, function(r) { sum(r == -1) } )
  res = data.frame(ProbeID = rownames(data), up=countUp, down=countDown)
  res$up = 100.0 * (res$up / ncol(data))
  res$down = 100.0 * (res$down / ncol(data))
  res = merge(res, SNPpos, by.x="ProbeID", by.y="row.names")
  return(res)
}



p1.count = countNonzero(score_file[, p1.sampleID])
p2.count = countNonzero(score_file[, p2.sampleID])
p3.count = countNonzero(score_file[, p3.sampleID])


p1.count$what=rep("Basal", nrow(p1.count))
p2.count$what=rep("Exocrine", nrow(p2.count))
p3.count$what=rep("Classical", nrow(p3.count))


ptotal = rbind(p1.count, p2.count,p3.count)
ptotal$Chr = factor(as.character(ptotal$Chr), levels=c(as.character(1:22), "X"))

#write.table(ptotal,"ptotal_all_gp.txt")

library(ggplot2)

png(filename="../results/OUH.png", width=4096, height=1024)
ggplot(ptotal) + 
  geom_segment(color="red", aes(x=Position, xend=Position, y=0, yend = up)) +
  geom_segment(color="green", aes(x=Position, xend=Position, y=0, yend = -down)) +
  facet_grid(what~Chr , scales="free_x", space="free_x") +
  theme_bw() +
  theme(axis.text.x = element_blank() ) +
  ggtitle("PDAC subtypes") + xlab("Genomic position") + ylab("Frequency (%)")

####################################
resp_in= data.frame(SampleID= names(clusters$TCGA$meta_classes), Resp= as.numeric(clusters$TCGA$meta_classes))

resp_in=  resp_in[which(resp_in$SampleID %in%  colnames(score_file)),]

p1.sampleID = resp_in[which(resp_in$Resp==1), 1]
p2.sampleID = resp_in[which(resp_in$Resp==2), 1]
p3.sampleID = resp_in[which(resp_in$Resp==3), 1]
#p4.sampleID = resp_in[which(resp_in$Resp==4), 1]
#p5.sampleID = resp_in[which(resp_in$Resp==5), 1]

countNonzero <- function(data) {
  countUp = apply(data, 1, function(r) { sum(r == 1) } )
  countDown = apply(data, 1, function(r) { sum(r == -1) } )
  res = data.frame(ProbeID = rownames(data), up=countUp, down=countDown)
  res$up = 100.0 * (res$up / ncol(data))
  res$down = 100.0 * (res$down / ncol(data))
  res = merge(res, SNPpos, by.x="ProbeID", by.y="row.names")
  return(res)
}



p1.count = countNonzero(score_file[, p1.sampleID])
p2.count = countNonzero(score_file[, p2.sampleID])
p3.count = countNonzero(score_file[, p3.sampleID])


p1.count$what=rep("Basal", nrow(p1.count))
p2.count$what=rep("Exocrine", nrow(p2.count))
p3.count$what=rep("Classical", nrow(p3.count))


ptotal = rbind(p1.count, p2.count,p3.count)
ptotal$Chr = factor(as.character(ptotal$Chr), levels=c(as.character(1:22), "X"))

#write.table(ptotal,"ptotal_all_gp.txt")

library(ggplot2)

png(filename="../results/TCGA.png", width=4096, height=1024)
ggplot(ptotal) + 
  geom_segment(color="red", aes(x=Position, xend=Position, y=0, yend = up)) +
  geom_segment(color="green", aes(x=Position, xend=Position, y=0, yend = -down)) +
  facet_grid(what~Chr , scales="free_x", space="free_x") +
  theme_bw() +
  theme(axis.text.x = element_blank() ) +
  ggtitle("PDAC subtypes") + xlab("Genomic position") + ylab("Frequency (%)")
###################
