#'
#'
#'
#'
#'
#'
countCNAscoresByProbe <- function(CNAscoreDT, SNPposDT, annotSampMetaClassDT) {
  
  message("Getting CNA score by probe...")
  scoresByProbe <- melt(CNAscoreDT, id.vars="ProbeID", 
                        variable.name="sample", value.name="score")
  rm(CNAscoreDT)
  
  message("Annotating with metaclass data...")
  copyNumSampMetaclassDT <- merge(annotSampMetaClassDT, scoresByProbe, 
                                  by.x="samples", by.y="sample")
  rm(scoresByProbe)
  
  nSamples <- length(unique(copyNumSampMetaclassDT$samples))
  
  message("Calculating proportion of CNA scores by probe...")
  scoreByClassCount <- copyNumSampMetaclassDT[, .(proportion=nrow(.SD)/..nSamples), 
                                              by=.(metaClasses, ProbeID, score)]
  
  message("Reformatting probe score DT...")
  rm(copyNumSampMetaclassDT)
  scoreCountByProbe <- dcast(scoreByClassCount, 
                             metaClasses + ProbeID ~ score,
                             value.var="proportion")
  colnames(scoreCountByProbe) <- c("metaClass", "probe", "down", "same", "up")
  rm(scoreByClassCount)
  
  message("Annotating with SNP positions...")
  annotScoreByProbe <- merge(scoreCountByProbe, SNPposDT, by.x="probe", by.y="ID")
  rm(scoreCountByProbe)  
  
  return(annotScoreByProbe)
}

#'
#'
#'
#'
#'
plotPerChromosomeCNAscores <- function(CNAscoreDT) {
  CNAscoreDT <- copy(CNAscoreDT)
  CNAscoreDT[, `:=`(up=100 * up, down=100 * down, same=NULL)]
  
  splitDT <- split(CNAscoreDT, by='Chr')
  
  plots <- lapply(splitDT, 
                  function(DT) {
                    ggplot(DT) +
                      geom_line(aes(x=Position, y=up), colour="darkred") +
                      geom_line(aes(x=Position, y=-down), color="black") +
                      xlab("") + xlab("") +
                      ggtitle(paste0("Chromosome: ", unique(DT$Chr))) +
                      theme_classic() +
                      facet_grid(rows=vars(metaClass))
                  })
  return(plots)
}


