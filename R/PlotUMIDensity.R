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
PlotUMIDensity <- function(seuratObject, group = NULL, vlineIntercept = 500) {
  # Ensure the Seurat object is correctly provided
  if (!"Seurat" %in% class(seuratObject)) {
    stop("The provided object is not a Seurat object.")
  }

  # Extract metadata
  metadata <- seuratObject@meta.data

  # Ensure nCount_RNA is treated as numeric, if it's not, try to convert it
  metadata$nCount_RNA <- as.numeric(metadata$nCount_RNA)

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
    aes_mappings <- aes(x = nCount_RNA, fill = .data[[fill_color_column]], color = .data[[fill_color_column]])
  } else {
    aes_mappings <- aes(x = nCount_RNA)
  }

  # Create the base plot
  p <- ggplot(metadata, aes_mappings) +
    geom_density(alpha = 0.2) +
    scale_x_log10(limits = c(100, NA)) +
    theme_classic() +
    xlab("UMI Counts") +
    ylab("Cell density") +
    geom_vline(xintercept = vlineIntercept) +
    ggtitle("UMI counts (transcripts) per cell")


  # Return the plot
  return(p)
}
