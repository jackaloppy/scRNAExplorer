#' Plot Gene Count Density from a Seurat Object
#'
#' This function generates a density plot of gene counts per cell from a Seurat object,
#' allowing for visualization of cell distribution based on the number of detected genes.
#' The plot uses a logarithmic scale for the x-axis and includes a vertical line as a reference.
#'
#' @param seuratObject A Seurat object containing single-cell RNA-seq data.
#' @param sample Optional; the name of the metadata column to use for coloring and filling in the plot.
#'               Defaults to orig.ident if not specified.
#' @param vlineIntercept The x-intercept for the vertical line in the plot. Defaults to 300.
#' @return A ggplot object showing the density of gene counts per cell.
#' @import ggplot2
#' @import dplyr
#' @import Seurat
#' @examples
#' plotGeneDensity(seuratObject)
#' plotGeneDensity(seuratObject, group = "groups" , vlineIntercept = 300)
#' @export
plotGeneDensity <- function(seuratObject, group = NULL, vlineIntercept = 300) {
  # Ensure the Seurat object is correctly provided
  if (!"Seurat" %in% class(seuratObject)) {
    stop("The provided object is not a Seurat object.")
  }

  # Extract metadata
  metadata <- seuratObject@meta.data

  # Ensure nFeature_RNA is treated as numeric, if it's not, try to convert it
  metadata$nFeature_RNA <- as.numeric(metadata$nFeature_RNA)

  # Check for group or use a default column
  fill_color_column <- if (!is.null(group) && group %in% names(metadata)) {
    group
  } else if ('orig.ident' %in% names(metadata)) {
    'orig.ident'
  } else {
    NULL
  }

  # Define the aesthetic mappings based on the presence of a fill/color column
  if (!is.null(fill_color_column)) {
    aes_mappings <- aes(x = nFeature_RNA, fill = .data[[fill_color_column]], color = .data[[fill_color_column]])
  } else {
    aes_mappings <- aes(x = nFeature_RNA)
  }

  # Create the base plot
  p <- ggplot(metadata, aes_mappings) +
    geom_density(alpha = 0.2) +
    scale_x_log10(limits = c(100, NA)) +
    theme_classic() +
    xlab("Gene Counts") +
    ylab("Cell density") +
    geom_vline(xintercept = vlineIntercept) +
    ggtitle("Gene count per cell")


  # Return the plot
  return(p)
}
