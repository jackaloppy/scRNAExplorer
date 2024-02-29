#' Plot Novelty Score Density from a Seurat Object
#'
#' This function generates a density plot of the novelty score (log10(genes) detected per log10(UMI))
#' from a Seurat object, allowing for visualization of cell distribution based on this metric.
#' The plot includes a vertical line as a reference.
#'
#' @param seuratObject A Seurat object containing single-cell RNA-seq data.
#' @param group Optional; the name of the metadata column to use for coloring and filling in the plot.
#'              Defaults to orig.ident if not specified.
#' @param vlineIntercept The x-intercept for the vertical line in the plot. Defaults to 0.8.
#' @return A ggplot object showing the density of the novelty score.
#' @import ggplot2
#' @import dplyr
#' @import Seurat
#' @examples
#' plotNoveltyScoreDensity(seuratObject)
#' plotNoveltyScoreDensity(seuratObject, group = "sampleGroup", vlineIntercept = 0.8)
#' @export
plotNoveltyDensity <- function(seuratObject, group = NULL, vlineIntercept = 0.8) {
  # Ensure the Seurat object is correctly provided
  if (!"Seurat" %in% class(seuratObject)) {
    stop("The provided object is not a Seurat object.")
  }

  # Extract metadata
  metadata <- seuratObject@meta.data

  # Calculate log10GenesPerUMI if not already present
  if (!"log10GenesPerUMI" %in% names(metadata)) {
    metadata$log10GenesPerUMI <- log10(metadata$nFeature_RNA) / log10(metadata$nCount_RNA)
  }

  # Check for group or use a default column
  fill_color_column <- if (!is.null(group) && group %in% names(metadata)) {
    group
  } else if ('orig.ident' %in% names(metadata)) {
    'orig.ident'
  } else {
    NULL
  }

  if (!is.null(fill_color_column)) {
    aes_mappings <- aes(x = log10GenesPerUMI, fill = .data[[fill_color_column]], color = .data[[fill_color_column]])
  } else {
    aes_mappings <- aes(x = log10GenesPerUMI)
  }

  # Create the ggplot density plot
  p <- ggplot(metadata, aes_mappings) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = vlineIntercept) +
    xlim(0.7, NA) +
    ggtitle("Gene detected per UMI (novelty score)")

  # Return the plot
  return(p)
}
