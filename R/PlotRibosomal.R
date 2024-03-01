#' Plot Ribosomal Gene Percentage
#'
#' This function takes a Seurat object and plots the density of ribosomal gene percentages stored in the object's metadata.
#' It requires the `Seurat` and `ggplot2` packages. If the ribosomal gene percentages have not been pre-calculated and stored
#' in the object using the `PercentageFeatureSet` function with a pattern matching ribosomal genes (e.g., "^Rp[sl]"),
#' this function will first calculate these percentages.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA sequencing data.
#' @param group The name of the metadata column in the Seurat object that identifies the sample or condition.
#'              This will be used for coloring and filling the density plot. If `NULL`, no grouping is applied.
#' @return A ggplot object representing the density plot of ribosomal gene percentages.
#' @examples
#' # Assuming `seurat_obj` is your Seurat object and it has a column 'sample' in its metadata:
#' plotRibosomal(seurat_obj, group = "sample")
#' @export
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom ggplot2 ggplot aes geom_density scale_x_log10 theme_classic geom_vline ggtitle
PlotRibosomal <- function(seurat_obj, group = NULL) {

  # Calculate percent.ribo if not already present
  if (!"percent.ribo" %in% names(seurat_obj[[]])) {
    seurat_obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^(?i)Rp[sl]")
  }

  # Extract metadata
  metadata <- seurat_obj@meta.data

  # Generate the plot
  p <- ggplot(metadata, aes_string(x = "percent.ribo")) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    geom_vline(xintercept = 0.2, linetype="dashed") +
    ggtitle("Ribosomal Gene Percentage")

  # Check if grouping is specified and exists in metadata
  if (!is.null(group) && group %in% names(metadata)) {
    p <- p + aes_string(color = group, fill = group)
  } else if (!is.null(group)) {
    warning("Specified group column does not exist in metadata. Plotting without grouping.")
  }

  return(p)
}
