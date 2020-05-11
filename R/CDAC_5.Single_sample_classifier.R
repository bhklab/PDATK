########################################################################################################
### Required libraries
#########################################################################################################

library(randomForest)
library(vcdExtra)
library(switchBox)
library(impute)
########################################################################################################
### Load data
#########################################################################################################

load('../../data/common_genes_cohorts_new.RData')       
meta_genes=load('../../data/meta_genes.RData')
load('../data/clusters.RData')

########################################################################################################
#########################################################################################################

pcsi=data.matrix(cohorts1$pcsi[as.character(meta_genes),])
tcga=data.matrix(cohorts1$tcga[as.character(meta_genes),])
icgc=data.matrix(cohorts1$icgc_seq[as.character(meta_genes),])
kirby=data.matrix(cohorts1$kirby[as.character(meta_genes),])
ouh=data.matrix(cohorts1$ouh[as.character(meta_genes),])
winter=data.matrix(cohorts1$winter[as.character(meta_genes),])
collisson=data.matrix(cohorts1$collisson[as.character(meta_genes),])
zhang=data.matrix(cohorts1$zhang[as.character(meta_genes),])
chen=data.matrix(cohorts1$chen[as.character(meta_genes),])
unc=data.matrix(cohorts1$unc[as.character(meta_genes),])
icgc_arr=data.matrix(cohorts1$icgc_arr[as.character(meta_genes),])
balagurunathan=data.matrix(cohorts1$balagurunathan[as.character(meta_genes),])
pei=data.matrix(cohorts1$pei[as.character(meta_genes),])
grutzmann=data.matrix(cohorts1$grutzmann[as.character(meta_genes),])
badea=data.matrix(cohorts1$badea[as.character(meta_genes),])

hamidi=data.matrix(cohorts1$hamidi[as.character(meta_genes),])
haider=data.matrix(cohorts1$haider[as.character(meta_genes),])
bauer=data.matrix(cohorts1$bauer[as.character(meta_genes),])
yang=data.matrix(cohorts1$yang[as.character(meta_genes),])
#van_den_broeck=data.matrix(cohorts1$van_den_broeck[as.character(meta_genes),])
lunardi=data.matrix(cohorts1$lunardi[as.character(meta_genes),])
janky=data.matrix(cohorts1$janky[as.character(meta_genes),])
hamidi=impute.knn(hamidi)[[1]]
lunardi=impute.knn(lunardi)[[1]]


##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################

clusters1=clusters
classes=list()
sample_names=list()
cohorts=list()
for(i in 1:(length(clusters1))){
  
  classes[[i]]=clusters1[[i]]$meta_classes
  sample_names[[i]]=names(clusters1[[i]]$meta_classes)
  cohorts[[i]]=rep(names(clusters1)[i],length(names(clusters1[[i]]$meta_classes)))
}

meta_Cluster= data.frame(sample = unlist(sample_names), meta_class= unlist(classes), cohorts=unlist(cohorts))

##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################

#########


com_matrix= cbind( pcsi, tcga, icgc, kirby, ouh, winter, collisson, zhang, chen, unc, icgc_arr, balagurunathan, grutzmann, badea , pei, hamidi, yang, lunardi, janky, bauer, haider)
length(which(colnames(com_matrix) == meta_Cluster$sample))
#com_matrix= NaRV.omit(com_matrix)
#meta_Cluster$meta_class
#########################################################################


dat <- com_matrix
Groups <- meta_Cluster$meta_class
y <- as.numeric(as.character(Groups))

pheno_grp <- function(y, com){
  
  z=setdiff(y,com)
  class_com = unlist(sapply(z, function(x) which(Groups ==x)))
  na_index= which(y %in% NA)
  y1=as.numeric(y)
  y1[c(class_com, na_index)] =0
  return(as.factor(y1))  
}

n=3
#dat[is.na(dat)] <- 0
dat[is.na(dat)] <- 0

model=lapply(1:n, function(x) SWAP.KTSP.Train(dat, pheno_grp(y, x) , k=20))

gene1 = unlist(sapply(1:n, function(x) model[[x]]$TSPs[,1]))
gene2 = unlist(sapply(1:n, function(x) model[[x]]$TSPs[,2]))
genepairs <- as.character(paste(gene1,gene2,sep=">"))

feat <- c(gene1, gene2)
xx= dat[feat,]


ll=lapply(1: dim(xx)[2], function(x)  ifelse(xx[gene1,x] > xx[gene2,x], 1,0))
output <- matrix(unlist(ll), nrow = ncol(xx), byrow = TRUE)

rownames(output)=paste('Patient',1:ncol(xx),sep='') 
colnames(output)=genepairs

y <- as.factor(y)
y_NA = which(y %in% NA)
y1= y[-y_NA]
output1= output[-y_NA,]
output1=t(output1)


cores=detectCores()
cl <- makeCluster(cores[1]-1, outfile="") #not to overload your computer
registerDoParallel(cl)

pred=list()
pred= foreach(i = 1: dim(output1)[2], .packages = c( "randomForest","foreach")) %dopar% {         
  
  output2 <- t(output1)
  print(i)
  train <- output2[-i,]
  test <- output2[i,]
  ytrain <- y1[-i]
  model <- randomForest(train, ytrain)
  pred[[i]] <- as.character(predict(model, test))
  
}

binary_rf_model <- randomForest( t(output1), y1)
gg=list(gene1=gene1, gene2= gene2, genepairs =genepairs)

save(binary_rf_model, file="../results/binary_rf_model.RData")
save(output1, file="../results/output1_RF.RData")
save(y1, file="../results/RF_response.RData")
save(gg, file="../results/gg.RData")

pred1 <- do.call("c", pred)
table(pred1, as.character(y1))
a <- data.frame(predictions = pred1,group = y1)
a1 <- xtabs(~a[,1] + a[,2], data = a)
summary(assocstats(a1))
length(which(pred1 == as.character(y1)))/length(pred1)

        