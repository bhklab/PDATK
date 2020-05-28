#'
#'
#'
#'
boxplotCellarityByCohort <- function(mergedDT, saveDir, fileName) {
  plot <- ggplot(mergedDT, aes(cohorts, cellularity)) +
          geom_boxplot(aes(fill=metaClasses), width=0.5) +
          scale_fill_brewer(palette="Pastel1") +
          ylab("Cellularity") +
          guides(fill=guide_legend("Cellularity")) +
          theme(axis.text.x=element_text(angle=45, hjust=1))
  
  ## FIXME:: Refactor to utility function to reduce repetiion
  if(!missing(saveDir) && !missing(fileName)) {
    ggsave(file.path(saveDir, fileName), plot)
    message(paste0("Saved to ", file.path(saveDir, fileName)))
  }
  return(plot)
}