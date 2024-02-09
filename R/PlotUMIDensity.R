#' Plot UMI Counts Density from Seurat Object
#'
#' @param seuratObject A Seurat object.
#' @param sample Optional; the name of the metadata column to use for coloring.
#' @return A ggplot object showing the density of UMI counts per cell.
#' @examples
#' PlotUMIDensity(seuratObject)
#' PlotUMIDensity(seuratObject, sample = "sampleColumn")
#' @import ggplot2
#' @import dplyr
PlotUMIDensity <- function(seuratObject, sample = NULL) {
  # Extract metadata
  metadata <- seuratObject@meta.data

  # Create the base plot
  p <- ggplot(metadata, aes(x = nCount_RNA)) +
    geom_density(alpha = 0.2) +
    scale_x_log10(limits = c(100, NA)) +
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500) +
    ggtitle("UMI counts (transcripts) per cell")

  # If 'sample' column is specified and exists in metadata, use it for coloring
  if (!is.null(sample) && sample %in% names(metadata)) {
    p <- p + aes(color = .data[[sample]], fill = .data[[sample]])
  }

  # Return the plot
  return(p)
}
