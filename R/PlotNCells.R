#' Plot Number of Cells per Sample from Seurat Object Metadata
#'
#' @param seuratObject A Seurat object with a 'Group' column in the metadata.
#' @param sample Optional; the name of the metadata column to use for grouping.
#' @return A ggplot object showing the number of cells per sample.
#' @examples
#' PlotNCellsPerSample(seuratObject, sample = "Group")
#' @import ggplot2
#' @import dplyr
PlotNCells <- function(seuratObject, sample = "orig.ident") {
  # Extract metadata
  metadata <- seuratObject@meta.data

  # Check if sample column is specified and exists
  if (!is.null(sample) && sample %in% names(metadata)) {
    # Plot separate bars for each sample
    p <- ggplot(metadata, aes_string(x = sample, fill = sample)) +
      geom_bar() +
      theme_classic() +
      xlab(sample) +
      ylab("Count") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      ggtitle("Number of Cells per Sample")
  return(p)
  }
  ## test2
}

