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
#' @export
PlotUMIDensity <- function(seuratObject, group = NULL) {
  # Extract metadata
  metadata <- seuratObject@meta.data

  # Define the aesthetic mappings
  if (!is.null(group) && group %in% names(metadata)) {
    # If group column is specified and exists
    aes_mappings <- aes(x = nCount_RNA, fill = .data[[group]],
                        color = .data[[group]])
    x_label <- group
  } else {
    # Default to 'orig.ident' if group is not specified
    aes_mappings <- aes(x = orig.ident, fill = orig.ident, color = orig.ident)
    x_label <- "orig.ident"
  }

  # Create the base plot
  p <- ggplot(metadata, aes_mappings) +
    geom_density(alpha = 0.2) +
    scale_x_log10(limits = c(100, NA)) +
    theme_classic() +
    xlab(x_label) +
    ylab("Cell density") +
    geom_vline(xintercept = 500) +
    ggtitle("UMI counts (transcripts) per cell")


  # Return the plot
  return(p)
}
