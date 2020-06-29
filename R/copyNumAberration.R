#' Calculate the number of CNA probes which indicate copy number gain or loss
#'   feat each probe, annotated with the predicted meta-class of that sample.
#'
#' @param CNAscoreDT A \code{data.table} of probe by sample copy number calls,
#'    where 1 represents a gain, 0 normal and -1 a loss. Probes are included
#'    as a column since `data.table`s cannot have rownames.
#' @param SNFposDT A \code{data.table} with the 'Chr' holding chromosome number,
#'    'Positon' holding the genomic position of that SNP on the chromosome and
#'    'ID' representing the ID of that SNP probe.
#' @param annotSampleMetaClassDT A \code{data.table} with annotated meta-class
#'    predictions for each sample in each cohort.
#'
#' @return A \code{data.table} containing the sample meta-class predictions and
#'   number of gain, loss or normal probes in all datasets.
#'
#' @import data.table
#' @export
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

#' Plot the distribution of gains/losses per genomic position for each
#'   chromosome.
#'
#' @param CNAscoreDT A \code{data.table}  containing the sample meta-class
#'   predictions and number of gain, loss or normal probes in all datasets.
#'
#' @return A \code{list} of `ggplot`s, one for each chromsome, faceted by
#'   predicted meta-class/subtype showing the percenage of copy gain vs
#'   copy losses for each probe (genomic region).
#'
#' @importFrom ggplot2 geom_line xlab ylab ggtitle theme_classic facet_grid
#' @export
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


