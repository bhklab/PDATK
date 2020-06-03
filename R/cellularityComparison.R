#' Compare the cellularity between meta
#'
#' @param mergedDT A \code{data.table} of per sample cellularity and
#'    predicted meta-classes for several cohorts.
#' @param saveDir An optional \code{character} vector specifying the path
#'    to the directory where the plot should be saved. If excluded, fileName
#'    will not work.
#' @param fileName An optional \code{character} vector specifying the
#'    name and extension of the file to save the plot it. This is passed to
#'    the `ggplot2::ggsave`.
#'
#' @return A \code{ggplot} object with a grouped boxplot comparing the
#'    distribution of sample cellularity between the predicted meta-classes.
#'
#' @importFrom ggplot geom_boxplot scale_fill_brewer ylab guides
#' theme ggsave
#' @export
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